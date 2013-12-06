# setuptools.setup must be used for entry_points (needed on Windows),
# package_data, and others.
# at any rate, distutils.core.setup is deprecated now
try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

import os
from os.path import join # this isn't great, confused with str.join
import sys
import shutil
from gratelpy import get_version

py_v = sys.version_info
py_v_str = '.'.join('%s' % i for i in py_v[:2])
if py_v[:2] < (2, 6):
    print('GraTeLPy requires Python 2.6 (or higher) '
          'and you run Python %s.\n' % py_v_str)
    print('You need to update your Python version.\n'
          'More info here: http://python.org/download/releases/2.7.6/\n')
    sys.exit(1)

dependencies = [
    ('networkx', (1, 7, 0), 'https://pypi.python.org/pypi/networkx/')
]

os_name = os.name
if os_name in ('nt', 'dos'):
    os_name = 'windows'

if 'setuptools' not in sys.modules:
    print('Warning: setuptools was not detected for your Python setup.\n'
          'This means that package dependencies will not get installed '
          'automatically!\n'
          'We will now check which dependencies you miss ...\n')

    if os_name == 'windows':
        print('Since you are on Windows, we strongly recommend that you\n'
              'install a recent version of one of the popular Python\n'
              'distributions: Continuum Anaconda or Enthought Canopy.\n'
              'Installing Python packages on Windows by hand is a '
              'painful process ...\n'
              'Anaconda: http://continuum.io/downloads \n'
              'Enthought Canopy Express https://www.enthought.com/store/\n')

    interrupt = False

    for dep in dependencies:
        v_str = '.'.join('%d' % i for i in dep[1]) # (X,Y.Z) -> 'X.Y.Z'
        try:
            _ = __import__(dep[0])
        except ImportError:
            print('Could not find package %s in your Python setup.\n'
                  'Please install version %s or higher of %s.\n'
                  'Find out more here: %s\n' % (dep[0], v_str, dep[0], dep[2]))
            interrupt = True
        else:
            # 'X.Y.Z' -> (X,Y,Z)
            v_tuple = _.__version__.split('.') # asummes 'X.Y.Z' numbering
            v_tuple = tuple([int(s) for s in v_tuple])
	    if len(v_tuple) < 3:
                # version must be 'X.Y.0'
	        v_tuple = (v_tuple[0], v_tuple[1], 0)

            if v_tuple < dep[1]:
                print('Package %s was found in your Python setup, however\n'
                      'you have version %s installed, while GraTeLPy was\n'
                      'developed with version %s of this package.\n'
                      'Please consider upgrading %s to at least version %s.\n'
                      'More info on %s is here: %s\n' % (dep[0], _.__version__,
                                                          v_str, dep[0],
                                                          v_str, dep[0],
                                                          dep[2]))

        if interrupt:
            print('One ore multiple packages that GraTeLPy depends on are '
                  'not installed on your system.\n'
                  'Please install these packages first and then attempt '
                  'installing GraTeLPy again.\n'
                  'Manual installation of these packages can be avoided by '
                  'getting a recent version of setuptools first.\n'
                  'You can get setuptools here: '
                  'https://pypi.python.org/pypi/setuptools')
            sys.exit(1)

    print('We will now attempt to install GraTeLPy.\n'
          'Ignore any distutils warnings in the '
          'following output ("UserWarning")!')

    raw_input('\nPress Enter to continue ...\n')

version = get_version()

gratelpy_scripts = [join('bin', 'gratelpy_subclient'), 
                    join('bin', 'gratelpy_fragment_server'), 
                    join('bin', 'gratelpy_check_data'), 
                    join('bin', 'gratelpy_benchmark'),
                    join('bin', 'gratelpy_test'),
                    join('bin', 'gratelpy_time')]
                        
def windows_prepare(gratelpy_scripts):
    bat = ['@echo off\n',
          'set PYFILE=%~f0\n',
          'set PYFILE=%PYFILE:~0,-4%-script.py\n'
          '"python.exe" "%PYFILE%" %*\n']

    gratelpy_scripts_windows = []
    windows_bat_files = []
    
    for script in gratelpy_scripts:
        gratelpy_scripts_windows.append(script+'-script.py')
        f = open(script+'.bat', 'w')
        for string in bat:
            f.write(string)
        windows_bat_files.append(f.name)
        f.close()
        print('Wrote .bat file: %s.' % str(script+'.bat'))    
    for script, script_w in zip(gratelpy_scripts, gratelpy_scripts_windows):
        shutil.copy(script, script_w)
        print('Copied %s to %s.' % (script, script_w))
    return gratelpy_scripts_windows + windows_bat_files
                        
if os_name == 'windows':
    gratelpy_scripts = windows_prepare(gratelpy_scripts)
    
setup(
    name='GraTeLPy',
    version=version,
    author='Georg Walther and Matthew Hartley',
    author_email='gratelpy@gmail.com',
    packages=['gratelpy', 'gratelpy.tests'],
    package_dir = {'gratelpy': 'gratelpy',
                   'gratelpy.tests': join('gratelpy', 'tests')},
    package_data = {'gratelpy': [join('mechanisms', '*.txt')],
                    'gratelpy.tests': ['*.dat']},
    include_package_data = True,
    scripts=gratelpy_scripts,
    url='http://pypi.python.org/pypi/GraTeLPy',
    license='GPLv3',
    description='Graph theoretic linear stability analysis',
    long_description='GratTeLPy (Graph Theoretic Analysis of Linear Stability) is a software tool for parameter-free, graph-theoretic linear stability analysis. Given a mechanism file that describes a chemical reaction network (CRN) of mass-action reactions, GraTelPy analyzes the provided mechanism and determines if it meets a necessary condition for multistability.\nPlease find out more at https://github.com/gratelpy/gratelpy',
	classifiers= [
		'Development Status :: 3 - Alpha',
		'Environment :: Console',
		'Intended Audience :: Science/Research',
		'Intended Audience :: Education',
                'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
                'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
		'License :: OSI Approved',
		'Operating System :: POSIX :: Linux',
		'Programming Language :: Python :: 2',
		'Programming Language :: Python :: 2.7',
		'Topic :: Education',
		'Topic :: Scientific/Engineering',
		'Topic :: Scientific/Engineering :: Bio-Informatics',
		'Topic :: Scientific/Engineering :: Chemistry',
		'Topic :: Scientific/Engineering :: Mathematics',
		'Topic :: Scientific/Engineering :: Physics',
	],
    test_suite = 'gratelpy.tests.runtests',
    install_requires=[
        "networkx >= 1.7.0"
    ],
)

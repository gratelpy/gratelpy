# setuptools.setup must be used for entry_points (needed on Windows),
# package_data, and others.
# at any rate, distutils.core.setup is deprecated now
try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

import os
from os.path import join
import sys
import shutil

python_version = sys.version_info
if python_version[:2] < (2, 7):
    print('GraTeLPy requires Python 2.7')
    sys.exit(1)

os_name = os.name
if os_name in ('nt', 'dos'):
    os_name = 'windows'

version = __import__('gratelpy').get_version()

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
        with open(script+'.bat', "w") as f:
            for string in bat:
                f.write(string)
            windows_bat_files.append(f.name)
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
    author_email='Georg.Walther@jic.ac.uk',
    packages=['gratelpy', 'gratelpy.tests'],
    package_dir = {'gratelpy': 'gratelpy',
                   'gratelpy.tests': join('gratelpy', 'tests')},
    package_data = {'gratelpy': [join('mechanisms', '*.txt')],
                    'gratelpy.tests': ['*.dat']},
    include_package_data = True,
    scripts=gratelpy_scripts,
    url='http://pypi.python.org/pypi/GraTeLPy',
    license='BSD',
    description='Graph theoretic linear stability analysis',
    long_description='GratTeLPy (Graph Theoretic Analysis of Linear Stability) is a software tool for parameter-free, graph-theoretic linear stability analysis. Given a mechanism file that describes a chemical reaction network (CRN) of mass-action reactions, GraTelPy analyzes the provided mechanism and determines if it meets a necessary condition for multistability.\nPlease find out more at https://github.com/gratelpy/gratelpy',
	classifiers= [
		'Development Status :: 3 - Alpha',
		'Environment :: Console',
		'Intended Audience :: Science/Research',
		'Intended Audience :: Education',
		'License :: OSI Approved :: BSD License',
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
        "networkx >= 1.7.0",
        "numpy >= 1.6.2",
        "matplotlib >= 1.2.1",
        "setuptools >= 1.4"
    ],
)

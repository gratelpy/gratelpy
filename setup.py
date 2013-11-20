# setuptools.setup must be used for entry_points (needed on Windows),
# package_data, and others.
# at any rate, distutils.core.setup is deprecated now
try:
    from setuptools import setup
except:
    print('Setuptools is required for installation of GraTeLPy.')
    exit()

import os
from os.path import join
import shutil

version = __import__('gratelpy').get_version()

gratelpy_scripts = [join('bin', 'gratelpy_subclient'), 
                    join('bin', 'gratelpy_fragment_server'), 
                    join('bin', 'gratelpy_check_data'), 
                    join('bin', 'gratelpy_benchmark'),
                    join('bin', 'gratelpy_test')]
                        
def extension():
    if os.name == 'nt':
        return '.py'
    else:
        return ''

def windows_prepare(gratelpy_scripts):
    gratelpy_scripts_windows = []
    for script in gratelpy_scripts:
        gratelpy_scripts_windows.append(script+'.py')
    for script, script_w in zip(gratelpy_scripts, gratelpy_scripts_windows):
        shutil.copy(script, script_w)
    return gratelpy_scripts_windows 
                        
if os.name == 'nt':
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

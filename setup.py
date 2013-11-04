# Use preferably setuptools.setup as this will install required packages
# automatically. distutils.core.setup does not do so.
try:
    from setuptools import setup
except:
    from distutils.core import setup

version = __import__('gratelpy').get_version()

setup(
    name='GraTeLPy',
    version=version,
    author='Georg Walther and Matthew Hartley',
    author_email='Georg.Walther@jic.ac.uk',
    packages=['gratelpy', 'tests', 'mechanisms'],
    include_package_data = True,
    scripts=['bin/gratelpy_subclient', 'bin/gratelpy_fragment_server', 
             'bin/gratelpy_check_data', 'bin/gratelpy_benchmark',
             'bin/gratelpy_test'],
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
    install_requires=[
        "networkx >= 1.7.0",
        "numpy >= 1.6.2",
        "matplotlib >= 1.2.1"
    ],
)

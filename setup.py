from distutils.core import setup

version = __import__('gratelpy').get_version()

setup(
    name='GraTeLPy',
    version=version,
    author='Georg Walther and Matthew Hartley',
    author_email='Georg.Walther@jic.ac.uk',
    packages=['gratelpy'],
    scripts=['bin/subclient.py'],
    url='http://pypi.python.org/pypi/GraTeLPy',
    license='BSD',
    description='Graph theoretic linear stability analysis',
    long_description=open('README.md').read(),
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
        "networkx >= 1.7.0"
    ],
)

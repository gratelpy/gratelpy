from distutils.core import setup

setup(
    name='GraTeLPy',
    version='0.1.0',
    author='Georg Walther and Matthew Hartley',
    author_email='Georg.Walther@jic.ac.uk',
    packages=['gratelpy'],
    scripts=['bin/subclient.py'],
    url='http://pypi.python.org/pypi/GraTeLPy',
    license='LICENSE.txt',
    description='Graph theoretic linear stability analysis',
    long_description=open('README.txt').read(),
    install_requires=[
        "networkx >= 1.7.0"
    ],
)

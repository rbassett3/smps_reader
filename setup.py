#demo setup file

from setuptools import setup, find_packages
import sys, os

this_dir = os.path.dirname(os.path.realpath(__file__))

__version__ = '0.0.1'


setup(
    name='smps_reader',
    version=__version__,
    author='Robert Bassett',
    author_email='robert.bassett@nps.edu',
    description='An SMPS (Stochastic Mathematical Programming System) file parser',
    url='https://faculty.nps.edu/rbassett/',
    python_requires='>=3.0',
    install_requires=["numpy", "scipy", "mps_reader @ git+https://github.com/rbassett3/mps_reader"],
    packages=find_packages(exclude=['tests']),
    zip_safe=True,

)


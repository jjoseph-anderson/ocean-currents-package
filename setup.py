from setuptools import setup, find_packages

setup(
    name='ocean_currents',
    version='0.1',
    packages=find_packages(),
    install_requires=[
        'numpy',
        'matplotlib',
        'scipy',
        'h5py'
    ],
)
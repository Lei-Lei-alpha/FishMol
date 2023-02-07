#!/usr/bin/env python

from setuptools import setup, find_packages

setup(
    name = 'fishmol',
    version = '0.0.1',
    description='Molecular Dynamics Data Analysis Package',
    author='Lei Lei',
    author_email='Lei.Lei@durham.ac.uk',
    license='MIT',
    url='https://github.com/Lei-Lei-alpha/FishMol',
    packages=find_packages(),
    install_requires=[
    'numpy>=1.14.5',
    'recordclass>=0.17.2',
    'ase>=3.22.1'
    ],
    extras_require={'plotting': ['matplotlib>=3.5.1', 'jupyter']},
    entry_points={
    'console_scripts': [
        'fishmol = fishmol.__main__:main',
        ],
    },
    include_package_data=True,
    package_data={
        'examples': ['../examples/*.ipynb'],
        'docs': ['../Docs/fishmo_doc.md'],
        'test': ['../test/*.py'],
        'license': ['../LICENCE']
        },
     )
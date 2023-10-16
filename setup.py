#!/usr/bin/env python

from setuptools import setup, find_packages

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

modcontents = open('nemopy/__init__.py').read()

setup(
    name = 'astroNEMO',
    version = open('VERSION').read().strip(),
    author = 'Peter Teuben',
    author_email = 'teuben@astro.umd.edu',
    description = 'Python tools to aid in using NEMO in python',
    license = 'BSD',
    keywords = 'astronomy',
    url = 'https://github.com/teuben/nemo',

    packages = find_packages(),
    # packages = ['astroNEMO'],
    install_requires = ['numpy'],
    # package_data = { '.': ['nemo_start.py'] },

    # long_description=read('README.md'),
 
)

#!/usr/bin/env python

#from distutils.core import setup
from setuptools import setup, find_packages, find_namespace_packages

setup(name='Leavitt',
      version='1.0.2',
      description='Variable star fitting code',
      author='Kyle Matt, David Nidever, Pol Massana',
      author_email='dnidever@montana.edu',
      url='https://github.com/NideverAstroResearch/leavitt',
      packages=find_namespace_packages(where="python"),
      package_dir={"": "python"}, 
      #scripts=['bin/leavitt'],
      install_requires=['numpy','astropy(>=4.0)','scipy','dlnpyutils(>=1.0.3)',
                        'numba','astro-datalab','statsmodels','emcee','healpy'],
      include_package_data=True,
)

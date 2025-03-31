from __future__ import print_function

from distutils.dep_util import newer
import os, os.path
import setuptools
import subprocess
import sysconfig
import sys

with open("README.md", "r") as fh:
    long_description = fh.read()

with open("cytrack/VERSION", "r") as fh:
    version_ = fh.read()



setuptools.setup(
    name="cytrack",
    version=version_,
    distclass=conda_build.bdist_conda.CondaDistribution,
    developer="Albenis Pérez-Alarcón",
    CoDevelopers="Patricia Coll-Hidalgo, Ricardo M. Trigo, Raquel Nieto, and Luis Gimeno ",
    author_email="albenis.pérez.alarcon@uvigo.es",
    description="CyTRACK is an open-source and user-friendly Python toolbox for detecting and tracking cyclones",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/apalarcon/CyTRACK/",
    packages=setuptools.find_packages(),
    classifiers=[
        "Development Status :: Alpha",
        "Programming Language :: Python :: 3+",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    install_requires=["numpy","mpi4py","netCDF4","scipy","matplotlib","xarray","scikit-learn"],

    
    include_package_data=True,
    package_data={"":['VERSION', "cytrack_inputs.cfg","_version.py","cytrack_functions.py", "LAST_UPDATE", "ERA5_terrain_high.nc"]},
)

from pathlib import Path

from setuptools import find_packages, setup
module_dir = Path(__file__).resolve().parent

setup(
    name="Export_traj_to_mongoDB",
    version="0.1",
    url="",
    license="MIT",
    author="JPBergsma",
    author_email="johannes.bergsma@epfl.ch",
    description="A tool to extract data from trajectories and put it into a mongoDB",
    long_description=open(module_dir.joinpath("README.md")).read(),
    long_description_content_type="text/markdown",
    keywords="optimade mongodb trajectory",
    include_package_data=True,
    packages=find_packages(),
    classifiers=[
        "Development Status :: 2 - Pre-Alpha",
        "Programming Language :: Python :: 3 :: Only",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Intended Audience :: Developers",
        "Topic :: Database",
        "Topic :: Database :: Database Engines/Servers",
    ],
    python_requires=">=3.8",
    install_requires=[
        "pymongo>=3.12.1",
        "h5py>=3.7.0",
        "mdanalysis>=2.2.0",
        "numpy>=1.23",
        "pymatgen>=2022.7",
        "optimade @ git+https://github.com/JPBergsma/optimade-python-tools.git@optimade_python_tools_trajectory_0.1"
    ]
)

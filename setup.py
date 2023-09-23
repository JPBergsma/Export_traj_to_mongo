from pathlib import Path

from setuptools import find_packages, setup
module_dir = Path(__file__).resolve().parent

setup(
    name="Export_traj_to_mongoDB",
    version="0.2",
    url="",
    license="MIT",
    author="JPBergsma",
    author_email="kwibus@chemist.com",
    description="A tool to extract data from trajectories and put it into a mongoDB",
    long_description=open(module_dir.joinpath("README.md")).read(),
    long_description_content_type="text/markdown",
    keywords="optimade mongodb trajectory",
    include_package_data=True,
    packages=find_packages(),
    classifiers=[
        "Development Status :: 2 - Pre-Alpha",
        "Programming Language :: Python :: 3 :: Only",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Intended Audience :: Developers",
        "Topic :: Database",
        "Topic :: Database :: Database Engines/Servers",
    ],
    python_requires=">=3.9",
    install_requires=[
        "pymongo>=3.12.1",
        "MDAnalysis>=2.2.0",
        "numpy>=1.23",
        "pymatgen>=2022.7",
        "fastapi>=0.103.0",  # Todo figure out whether this should not be installed already with the installation of the optimade python tools.
        "jsonlines>=3.1.0",
        "optimade@git+https://github.com/JPBergsma/optimade-python-tools.git@JPBergsma/trajectory_endpoint"
    ]
)

import pathlib
from setuptools import setup, find_packages
from avadx import __version__, name

HERE = pathlib.Path(__file__).parent
README = (HERE / "README.md").read_text()

setup(
    name=name,
    version=__version__,
    keywords="AVA,Dx (Analysis of Variation for Association with Disease)",
    description="""a python package implementing the a AVA,Dx pipeline - a computational method for defining the functional role of DNA variation in complex diseases""",
    long_description=README,
    long_description_content_type="text/markdown",
    url="https://bitbucket.org/bromberglab/avadx",
    author="Maximilian Miller",
    author_email="mmiller@bromberglab.com",
    license="NPOSL-3.0",
    python_requires='>=3.8',
    include_package_data=True,
    packages=find_packages(),
    install_requires=[
        'requests==2.24.0',
        'docker==4.2.2',
    ],
    entry_points={
        'console_scripts': ['avadx=avadx.__main__:main'],
    },
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Natural Language :: English",
        "Operating System :: OS Independent"
    ],
    project_urls={
        "Bug Tracker": "https://bitbucket.org/bromberglab/avadx/issues",
        "Documentation": "https://bitbucket.org/bromberglab/avadx/wiki/docs",
        "Source Code": "https://bitbucket.org/bromberglab/avadx",
    }
)

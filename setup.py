#!/usr/bin/env python
from setuptools import setup, find_packages


setup(
    use_scm_version=True,
    setup_requires=["setuptools_scm", "setuptools_scm_git_archive"],
    install_requires=["setuptools_scm", "pandas==2.2", "Bio==1.7.1"],
    tests_require=["pytest", "pytest-cov", "flake8", "pep257", "black"],
    name="VirMet",
    description="Viral metagenomics in clinical applications",
    url="http://github.com/medvir/VirMet",
    author="Osvaldo Zagordi",
    author_email="firstname.lastname@gmail.com",
    packages=find_packages("src"),  # include all packages under src
    package_dir={"": "src"},  # tell setuptools packages are under src
    data_files=[("Rscripts", ["Rscripts/covplot.R"])],
    entry_points={"console_scripts": ["virmet = virmet.cli:main"]},
    license="MIT",
    long_description="""
    Set of tools for the analysis of sequencing data to identify and characterize
    the viral fraction in metagenomic samples, especially in the clinical setting.
    """,
    classifiers=[
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        "Development Status :: 5 - Production/Stable",
        # Indicate who your project is intended for
        "Intended Audience :: Developers",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Medical Science Apps.",
        # Pick your license as you wish (should match "license" above)
        "License :: MIT",
        # Specify the Python versions you support here. In particular, ensure
        # that you indicate whether you support Python 2, Python 3 or both.
        # 'Programming Language :: Python :: 2',
        # 'Programming Language :: Python :: 2.6',
        # 'Programming Language :: Python :: 2.7',
        # 'Programming Language :: Python :: 3',
        # 'Programming Language :: Python :: 3.2',
        # 'Programming Language :: Python :: 3.3',
        "Programming Language :: Python :: 3.5",
        "Programming Language :: Python :: 3.6",
    ],
)

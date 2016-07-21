from setuptools import setup, find_packages
setup(
    name = "VirMet",
    version = "1.0rc1",
    packages = find_packages(),
    scripts = ['bin/virmet'],
    data_files = [('Rscripts', ['Rscripts/covplot.R'])],
    entry_points={
        'console_scripts': [
            'virmet = virmet.__main__:main'
        ]
# 'gui_scripts': [
#     'baz = my_package_gui:start_func',
# ]
    },

    package_data = {
        # If any package contains *.txt or *.rst files, include them:
        '': ['*.txt', '*.rst'],
        # And include any *.msg files found in the 'hello' package, too:
        'hello': ['*.msg'],
    },

    # metadata for upload to PyPI
    author = "Osvaldo Zagordi",
    author_email = "firstname.lastname@gmail.com",
    description = "Viral metagenomics in clinical applications",
    license = "MIT",
    download_url = "http://github.com/ozagordi/VirMet",
    url = "http://virmet.readthedocs.io/en/latest/",
    long_description = '''
    Set of tools for the analysis of sequencing data to identify and characterize
    the viral fraction in metagenomic samples, especially in the clinical setting.
    '''
)

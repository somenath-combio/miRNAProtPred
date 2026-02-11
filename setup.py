from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="mirnaprotpred",
    version="0.1.6",
    author="Sudipta Sardar",
    author_email="sudipta@pusan.ac.kr",
    description="A pip installable CLI for miRNAProtPred tool",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/somenath-combio/mirnaprotpred",
    packages=find_packages(),
    package_data={
        'mirnaprotpred.SeqFinder': ['data/*.xlsx'],
        'mirnaprotpred.validator': ['data/*.xlsx'],
    },
    include_package_data=True,
    entry_points={
        'console_scripts': [
            'SeqFinder=mirnaprotpred.SeqFinder.seqfinder:cli',
            'validator=mirnaprotpred.validator.validator:cli',
        ],
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.13",
    install_requires=[
        "pandas>=2.3.3",
        "openpyxl>=3.1.5",
        "biopython>=1.86",
        "pyfiglet>=1.0.4",
        "requests>=2.32.5",
        "viennarna>=2.7.1",
    ],
)

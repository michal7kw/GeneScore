from setuptools import setup, find_packages
import os

# Read the contents of README.md
this_directory = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name="gene_scoring",
    version="0.1.4",
    packages=find_packages(),
    install_requires=[
        "numpy>=1.19.0",
        "pandas>=1.0.0",
        "scipy>=1.5.0",
    ],
    extras_require={
        "gpu": ["cupy>=9.0.0"],
    },
    author="Michal Kubacki",
    author_email="michal.kubacki@example.com",  # Replace with your email
    description="Gene set scoring package for single-cell RNA sequencing data with GPU acceleration.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/michal7kw/GeneScore",  # Replace with your actual repository URL
    project_urls={
        'Documentation': 'https://github.com/michal7kw/GeneScore#readme',
        'Bug Reports': 'https://github.com/michal7kw/GeneScore/issues',
        'Source Code': 'https://github.com/michal7kw/GeneScore',
    },
    keywords=['bioinformatics', 'single-cell', 'RNA-seq', 'gene scoring', 'GPU acceleration'],
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Natural Language :: English",
    ],
    python_requires=">=3.7",
    include_package_data=True,
    zip_safe=False,
    package_data={
        '': ['*.png', '*.jpg', '*.jpeg'],  # Include all image files
    },
)

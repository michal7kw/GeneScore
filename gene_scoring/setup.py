from setuptools import setup, find_packages

setup(
    name="gene_scoring",
    version="0.1.0",
    packages=find_packages(),
    install_requires=[
        "numpy>=1.19.0",
        "pandas>=1.0.0",
        "scipy>=1.5.0",
    ],
    extras_require={
        "gpu": ["cupy>=9.0.0"],
    },
    author="Original Authors",
    author_email="",
    description="A package for calculating gene scores in single-cell RNA sequencing data with GPU acceleration support",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="",
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
    ],
    python_requires=">=3.7",
)

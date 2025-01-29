"""
Gene Scoring Package

A package for calculating gene scores in single-cell RNA sequencing data with GPU acceleration support.
"""

from .core import calculate_gene_variances, score_genes

__version__ = "0.1.0"
__all__ = ["calculate_gene_variances", "score_genes"]

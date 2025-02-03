# GeneScore: Advanced Gene Set Scoring for Single-Cell RNA Sequencing

GeneScore is a sophisticated Python package that implements a novel gene set scoring methodology for single-cell RNA sequencing data analysis. Building upon techniques from popular toolkits like Seurat and Scanpy, GeneScore introduces several key improvements for more accurate and efficient gene set scoring.

## Key Features

- **Weighted Gene Contributions**: Each gene in the evaluated module is assigned an appropriate weight based on its correlation with the phenotype of interest.
- **Individualized Background Selection**: Background genes are selected independently for each gene of interest, providing more accurate control.
- **Variance-scaled Gene Contribution**: Gene contributions are scaled by their variance under examined conditions.
- **GPU-Accelerated Computation**: Optimized score calculation on GPUs for significantly improved performance.

## Mathematical Foundation

The scoring methodology is defined by the following formula:

```
S_c = Σ (w_l / (σ_c,l + ε)) * (X_c,l - X̄_l,control)
```

Where:
- n: number of genes in the gene list
- w_l: weight for gene l
- X_c,l: mean expression value of gene l in condition c
- X̄_l,control: average expression of control genes in the background
- σ_c,l: standard deviation of gene expression in condition c
- ε: small constant to avoid division by zero

## Installation

```bash
pip install gene_scoring
```

For GPU support:
```bash
pip install gene_scoring[gpu]
```

## Usage

Basic usage example:

```python
from gene_scoring import GeneScorer

# Initialize the scorer
scorer = GeneScorer(
    normalize_weights=True,
    abs_diff=False,
    weighted=True
)

# Calculate scores
scores = scorer.calculate_scores(
    expression_matrix,
    gene_list,
    weights=gene_weights
)
```

## Parameters

- `normalize_weights`: Normalize weights to sum to 1
- `abs_diff`: Use absolute difference between expression and background
- `weighted`: Use weighted scoring (if False, all weights = 1)

## Requirements

- Python ≥ 3.7
- NumPy ≥ 1.19.0
- Pandas ≥ 1.0.0
- SciPy ≥ 1.5.0
- CuPy ≥ 9.0.0 (optional, for GPU support)

## License

MIT License

## Citation

If you use this package in your research, please cite:

[Your paper citation here]

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

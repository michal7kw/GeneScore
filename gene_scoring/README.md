# Gene Scoring Package

This is a Python package implementing a novel gene set scoring methodology for single-cell RNA sequencing data analysis. It builds upon and extends techniques from widely used single-cell analysis toolkits such as Seurat and Scanpy, introducing several key improvements for more accurate gene set scoring.

![alt text](https://raw.githubusercontent.com/michal7kw/GeneScore/main/gene_scoring/image.png)

## Key Features

1. **Weighted Gene Contributions**
   - Each gene in the evaluated module is assigned an appropriate weight based on its correlation with the phenotype of interest
   - Weights can be normalized to sum to 1 (optional)

2. **Individualized Background Selection**
   - Background genes are selected independently for each gene of interest
   - Control features are selected from the same expression bin in the histogram
   - Background can be based on control sample, entire dataset, or a selected pool of genes

3. **Variance-scaled Gene Contribution**
   - Individual gene contributions are scaled by their variance under examined conditions
   - Helps account for gene expression variability in different conditions

4. **Optimized Computation**
   - GPU-accelerated score calculation for improved performance
   - Support for sparse matrices and chunked processing for large datasets
   - Efficient background gene selection algorithm

## Installation

### Basic Installation

```bash
pip install gene_scoring
```

### Installation with GPU Support

```bash
pip install gene_scoring[gpu]
```

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

## Usage

```python
from gene_scoring import GeneScorer

# Initialize the scorer
scorer = GeneScorer(
    normalize_weights=True,  # Normalize weights to sum to 1
    abs_diff=False,         # Use signed difference (not absolute)
    weighted=True           # Use weighted scoring
)

# Calculate scores
scores = scorer.calculate_scores(
    expression_matrix,    # Gene expression matrix (genes x cells)
    gene_list,           # List of genes to score
    weights=gene_weights  # Optional weights for each gene
)
```

## Parameters

The main parameters for score calculation include:

- `normalize_weights`: If True, weights are normalized to sum to 1
- `abs_diff`: If True, use absolute difference between expression and background
- `weighted`: If True, use weighted scoring (otherwise all weights = 1)
- `ctrl_size`: Number of control genes to use per target gene
- `gpu`: Whether to use GPU acceleration
- `chunk_size`: Size of chunks for processing large datasets
- `random_state`: Random seed for reproducibility

For detailed parameter descriptions and additional options, see the function documentation.

## Requirements

- Python ≥ 3.7
- NumPy ≥ 1.19.0
- Pandas ≥ 1.0.0
- SciPy ≥ 1.5.0
- CuPy ≥ 9.0.0 (optional, for GPU support)

## License

This project is licensed under the MIT License - see the LICENSE file for details.

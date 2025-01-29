# Gene Scoring Package

A Python package for calculating gene scores in single-cell RNA sequencing data with optional GPU acceleration support.

## Features

- Calculate gene set scores with control gene comparison
- GPU acceleration support using CuPy
- Flexible gene weighting options
- Support for sparse matrices
- Condition-specific scoring
- Variance-based scaling
- Chunked processing for large datasets

## Installation

### Basic Installation

```bash
pip install gene_scoring
```

### Installation with GPU Support

```bash
pip install gene_scoring[gpu]
```

## Usage

```python
import scanpy as sc
from gene_scoring import score_genes

# Load your data
adata = sc.read_h5ad("your_data.h5ad")

# Define your gene list
genes_of_interest = ["GENE1", "GENE2", "GENE3"]

# Optional: Define gene weights
gene_weights = {
    "GENE1": 1.0,
    "GENE2": 0.8,
    "GENE3": 1.2
}

# Calculate scores
score_genes(
    adata=adata,
    gene_list=genes_of_interest,
    gene_weights=gene_weights,
    score_name="my_score",
    gpu=True  # Set to False if no GPU is available
)

# Access scores
scores = adata.obs["my_score"]
```

## Parameters

The main function `score_genes` accepts the following parameters:

- `adata`: AnnData object containing gene expression data
- `gene_list`: List of genes to score
- `gene_weights`: Optional dictionary of gene weights
- `score_name`: Name for the resulting score (default: "score")
- `ctrl_size`: Number of control genes to use (default: 50)
- `gpu`: Whether to use GPU acceleration (default: True)
- And many more customization options...

For full parameter descriptions, see the function documentation.

## Requirements

- Python ≥ 3.7
- NumPy ≥ 1.19.0
- Pandas ≥ 1.0.0
- SciPy ≥ 1.5.0
- CuPy ≥ 9.0.0 (optional, for GPU support)

## License

This project is licensed under the MIT License - see the LICENSE file for details.

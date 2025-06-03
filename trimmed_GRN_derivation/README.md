# Gene Regulatory Network (GRN) Derivation for Gene Set Scoring

This directory contains a series of Python scripts designed to derive gene sets for gene set scoring based on Gene Regulatory Networks (GRNs). The process involves single-cell RNA sequencing (scRNA-seq) and ATAC-seq data processing, peak calling, co-accessibility analysis, motif scanning, and GRN inference using the CellOracle framework.

## Scripts Description:

## Key Thresholding Parameters

This section outlines important thresholding parameters used throughout the pipeline that can influence the final gene sets.

### P-value and False Positive Rate (FPR) Thresholds:

1.  **Motif Scanning (in [`7_scan_motifs_for_different_cell_types.py`](trimmed_GRN_derivation/7_scan_motifs_for_different_cell_types.py:1)):**
    *   A `fpr_threshold` (False Positive Rate) of `0.05` is used during the motif scanning step ([`tfi.scan(fpr=fpr_threshold, ...)`](trimmed_GRN_derivation/7_scan_motifs_for_different_cell_types.py:246)). This threshold determines the stringency for identifying transcription factor binding motifs in the accessible chromatin regions.

2.  **GRN Link Filtering (in [`9_create_scores_procedural.py`](trimmed_GRN_derivation/9_create_scores_procedural.py:1)):**
    *   A p-value threshold of `0.05` is used to filter regulatory links ([`links.filter_links(p=0.05, ...)`](trimmed_GRN_derivation/9_create_scores_procedural.py:411)). This threshold helps in retaining statistically significant regulatory interactions in the gene regulatory network.

### Other Important Thresholds Limiting Gene Set Size:

**1. scRNA-seq Processing ([`1_herring_scRNA_processing.py`](trimmed_GRN_derivation/1_herring_scRNA_processing.py:1)):**
    *   **Cell Quality Control:** Filters cells based on:
        *   Minimum percentile for genes per cell: `min_genes_percentile = 2` ([`1_herring_scRNA_processing.py:83`](trimmed_GRN_derivation/1_herring_scRNA_processing.py:83))
        *   Maximum percentile for genes per cell: `max_genes_percentile = 98` ([`1_herring_scRNA_processing.py:84`](trimmed_GRN_derivation/1_herring_scRNA_processing.py:84))
        *   Minimum percentile for total counts per cell: `min_counts_percentile = 2` ([`1_herring_scRNA_processing.py:85`](trimmed_GRN_derivation/1_herring_scRNA_processing.py:85))
        *   Maximum percentile for total counts per cell: `max_counts_percentile = 98` ([`1_herring_scRNA_processing.py:86`](trimmed_GRN_derivation/1_herring_scRNA_processing.py:86))
        *   Maximum mitochondrial gene percentage: `max_mito_percent = 30` ([`1_herring_scRNA_processing.py:87`](trimmed_GRN_derivation/1_herring_scRNA_processing.py:87))
    *   **Gene Filtering:** Removes genes expressed in fewer than `min_cells = 1` ([`sc.pp.filter_genes`](trimmed_GRN_derivation/1_herring_scRNA_processing.py:280)).
    *   **Highly Variable Gene Selection:** Selects the top `n_top_genes=3000` ([`sc.pp.highly_variable_genes`](trimmed_GRN_derivation/1_herring_scRNA_processing.py:292)).

**2. ATAC-seq Processing ([`2_herring_atac_processing.py`](trimmed_GRN_derivation/2_herring_atac_processing.py:1)):**
    *   **Peak Calling (MACS2):** Uses a q-value (FDR) of `0.05` ([`peak_calling`](trimmed_GRN_derivation/2_herring_atac_processing.py:313)).
    *   **Consensus Peak Generation:** Defines `peak_half_width = 250` ([`get_consensus_peaks`](trimmed_GRN_derivation/2_herring_atac_processing.py:332)).
    *   **Barcode Quality Control:** Uses automatic thresholds in `pycisTopic` ([`get_barcodes_passing_qc_for_sample`](trimmed_GRN_derivation/2_herring_atac_processing.py:450)).

**3. Cell-Type Consensus Regions ([`3_celltype_consensus_regions.py`](trimmed_GRN_derivation/3_celltype_consensus_regions.py:1)):**
    *   Filters consensus peaks using `score_threshold = 1.0` ([`3_celltype_consensus_regions.py:65`](trimmed_GRN_derivation/3_celltype_consensus_regions.py:65)).

**4. Co-accessibility Analysis ([`4_coaccess.py`](trimmed_GRN_derivation/4_coaccess.py:1)):**
    *   Filters co-accessible peak pairs if their absolute correlation is not greater than `0.2` ([`4_coaccess.py:147`](trimmed_GRN_derivation/4_coaccess.py:147)).

**5. TSS Integration with Co-accessible Regions ([`6_integrate_tss_with_coaccess_regions.py`](trimmed_GRN_derivation/6_integrate_tss_with_coaccess_regions.py:1)):**
    *   Filters integrated TSS-peak connections if the `coaccess` score is less than `0.6` ([`6_integrate_tss_with_coaccess_regions.py:92`](trimmed_GRN_derivation/6_integrate_tss_with_coaccess_regions.py:92)).

**6. Motif Scanning ([`7_scan_motifs_for_different_cell_types.py`](trimmed_GRN_derivation/7_scan_motifs_for_different_cell_types.py:1)):**
    *   Filters scanned TF binding motifs using `motif_score_threshold = 6` ([`7_scan_motifs_for_different_cell_types.py:227`](trimmed_GRN_derivation/7_scan_motifs_for_different_cell_types.py:227)).

**7. GRN Inference and Scoring ([`9_create_scores_procedural.py`](trimmed_GRN_derivation/9_create_scores_procedural.py:1)):**
    *   **Initial Gene Selection:** Analysis is limited to highly variable genes and predefined genes of interest ([`9_create_scores_procedural.py:102`](trimmed_GRN_derivation/9_create_scores_procedural.py:102)).
    *   **GRN Link Filtering:** Keeps up to the top `threshold_number=2000` regulatory links after p-value filtering ([`links.filter_links`](trimmed_GRN_derivation/9_create_scores_procedural.py:411)).

### 1. [`1_herring_scRNA_processing.py`](trimmed_GRN_derivation/1_herring_scRNA_processing.py)
*   **Purpose:** This script processes single-cell RNA sequencing (scRNA-seq) data. It handles data loading, initial setup, filtering for quality control, downsampling, normalization, and log transformation. It also generates pseudobulk profiles and visualizes gene expression across cell types.
*   **Implementation Explanation:**
    *   **Data Loading and Setup:** Loads scRNA-seq data using Scanpy and associated metadata. It sets up configurations for specific neuron types and age groups, handling both full counts and downsampled CPM data.
    *   **Filtering and Quality Control:** Filters cells based on chemistry version (v3), selected cell types, and age groups. It applies stringent quality control metrics, including the number of genes per cell, number of counts per cell, and mitochondrial gene percentage, removing outliers based on percentile thresholds.
    *   **Data Processing:** Downsamples large datasets to a manageable size (e.g., 30,000 cells), normalizes the data, and applies log transformation. Raw data is preserved, and different data layers are created for further analysis.
    *   **Pseudobulk and Visualization:** Groups cells by major cluster to calculate mean expression values, creating a pseudobulk AnnData object. It then visualizes the expression of genes of interest across different cell types using bar plots and histograms.
*   **Biological Significance:** This initial processing step is crucial for preparing high-quality scRNA-seq data, removing technical noise and low-quality cells. This ensures that downstream GRN inference and gene set scoring are based on reliable biological signals, allowing for accurate identification of cell-type-specific gene expression patterns.

### 2. [`2_herring_atac_processing.py`](trimmed_GRN_derivation/2_herring_atac_processing.py)
*   **Purpose:** This script processes single-nucleus ATAC-seq (snATAC-seq) data to identify open chromatin regions (peaks) and prepare them for downstream GRN analysis. It involves data exploration, importing and formatting ATAC data, generating pseudobulk profiles, inferring consensus peaks, and performing quality control.
*   **Implementation Explanation:**
    *   **Data Loading and Processing:** Loads ATAC metadata and fragment files. It processes cell data by extracting age and chemistry information, mapping cell types, and filtering based on selected cell types, ages, and quality control metrics (e.g., removing doublets using Scrublet).
    *   **Pseudobulk and Peak Calling:** Exports pseudobulk profiles for different cell types and performs peak calling using MACS2 to identify accessible chromatin regions.
    *   **Consensus Peak Inference:** Infers consensus peaks across all samples, merging overlapping peaks and filtering against a blacklist to remove spurious signals.
    *   **cisTopic Object Creation and QC:** Creates and merges `cisTopic` objects from fragment data, incorporating cell metadata. It performs extensive quality control on the ATAC-seq data, including plotting sample and barcode statistics, and identifies barcodes passing filters.
*   **Biological Significance:** This script identifies active regulatory regions in the genome (open chromatin) that are crucial for gene regulation. By generating cell-type-specific consensus peaks, it provides the foundational genomic coordinates for identifying transcription factor binding sites and understanding chromatin accessibility in a cell-type-specific manner, which is essential for building accurate GRNs.

### 3. [`3_celltype_consensus_regions.py`](trimmed_GRN_derivation/3_celltype_consensus_regions.py)
*   **Purpose:** This script refines the consensus peaks derived from ATAC-seq data by filtering them to create cell-type-specific consensus regions. These regions represent accessible chromatin elements that are active within particular cell populations.
*   **Implementation Explanation:**
    *   **Loading Consensus Regions:** Reads the `consensus_regions.bed` file generated by the previous ATAC processing step.
    *   **Peak Name Creation:** Generates a unique `peak_name` for each region based on its chromosomal coordinates.
    *   **Cell-Type Specific Filtering:** Iterates through predefined cell types and filters the consensus regions. It selects peaks that are specifically associated with each cell type (based on `peak_id` containing the cell type identifier) and meet a defined score threshold (e.g., `score >= 1.0`).
    *   **Output Generation:** Saves the filtered, cell-type-specific consensus regions into separate `.bed` files, organized by cell type.
*   **Biological Significance:** This step is critical for narrowing down the vast number of open chromatin regions to those most relevant to specific cell types. By focusing on cell-type-specific accessible regions, it enhances the precision of subsequent GRN inference, ensuring that regulatory interactions are considered within their appropriate cellular contexts. This allows for the identification of regulatory elements that drive cell identity and function.

### 4. [`4_coaccess.py`](trimmed_GRN_derivation/4_coaccess.py)
*   **Purpose:** This script calculates co-accessibility scores between pairs of genomic regions (peaks) within specific cell types. Co-accessibility indicates that two regions tend to be open or closed together, suggesting a functional relationship, often indicative of regulatory interactions.
*   **Implementation Explanation:**
    *   **Loading cisTopic Object:** Loads the `cistopic_obj.pkl` file, which contains the processed ATAC-seq fragment matrix and peak names.
    *   **Peak Name Formatting:** Formats peak names for consistency.
    *   **Cell-Type Iteration:** Iterates through each specified cell type.
    *   **Loading Cell-Type Regions:** Loads the cell-type-specific consensus regions generated by `3_celltype_consensus_regions.py`.
    *   **Peak Index Mapping:** Creates a mapping from peak IDs to their corresponding indices in the fragment matrix. It also loads or generates a list of peak indices relevant to the current cell type.
    *   **Co-accessibility Calculation:** For each cell type, it extracts the relevant portion of the fragment matrix. It then calculates the cosine similarity between the binarized (open/closed) fragment counts of peak pairs. This is done in chunks to manage memory usage.
    *   **Filtering and Saving:** Filters co-accessibility scores to retain only strong correlations (e.g., absolute correlation > 0.2) and maps the peak indices back to their original names. The results are saved as a CSV file for each cell type.
*   **Biological Significance:** Co-accessibility analysis helps identify distal regulatory elements (like enhancers) that physically interact with gene promoters or other regulatory regions. These interactions are crucial for understanding how gene expression is controlled in a three-dimensional nuclear space. The output of this script provides a list of functionally linked genomic regions, which are essential inputs for building comprehensive GRNs.

### 5. [`5_prepare_peak_data_for_celloracle.py`](trimmed_GRN_derivation/5_prepare_peak_data_for_celloracle.py)
*   **Purpose:** This script prepares the cell-type-specific consensus peak data into a format suitable for input into CellOracle, a tool used for inferring gene regulatory networks. Specifically, it extracts and formats peak identifiers.
*   **Implementation Explanation:**
    *   **Input File Selection:** Defines a dictionary `dir_input_files` that maps `neurons_set` configurations to lists of cell-type-specific consensus region `.bed` files.
    *   **File Processing Function (`process_file`):** Reads each specified `.bed` file, extracts the `chr_start_end` identifier (which is the `peak_name` column from the previous step), and collects these identifiers into a list.
    *   **Cell Type Extraction Function (`extract_cell_type`):** Parses the filename to extract the cell type name.
    *   **Batch Processing and Output:** Iterates through the selected input files, extracts the cell type, processes the file to get the peak identifiers, and then saves these identifiers into a new `.csv` file (e.g., `cell_type_peaks.csv`). The peaks are saved as a space-separated string on a single line.
*   **Biological Significance:** This script acts as a crucial data preparation step, transforming raw genomic coordinates of accessible regions into a standardized format that CellOracle can directly consume. This ensures compatibility and efficiency for the subsequent motif scanning and GRN inference steps, allowing CellOracle to accurately link regulatory regions to target genes.

### 6. [`6_integrate_tss_with_coaccess_regions.py`](trimmed_GRN_derivation/6_integrate_tss_with_coaccess_regions.py)
*   **Purpose:** This script integrates Transcription Start Site (TSS) information with the co-accessible regions identified in previous steps. This integration is vital for linking regulatory elements (peaks) to their putative target genes, forming the basis of gene regulatory networks.
*   **Implementation Explanation:**
    *   **Loading Co-accessibility and Peak Data:** For each cell type, it loads the co-accessibility data (`_coaccess.csv`) and the prepared peak data (`_peaks.csv`).
    *   **TSS Annotation:** Uses `celloracle.motif_analysis.get_tss_info` to annotate the loaded peaks with information about nearby Transcription Start Sites (TSSs) from the `hg19` reference genome. This step identifies which genes are associated with each peak.
    *   **Integration with Co-accessibility:** Integrates the TSS-annotated peaks with the co-accessibility connections using `celloracle.motif_analysis.integrate_tss_peak_with_cicero`. This links co-accessible peak pairs to their associated genes.
    *   **Filtering and Saving:** Filters the integrated data to retain only strong co-accessibility connections (e.g., `coaccess >= 0.6`). It then selects relevant columns (`peak_id` and `gene_short_name`) and saves the results as a CSV file (e.g., `processed_peak_file_cell_type.csv`) for each cell type.
*   **Biological Significance:** This is a critical step in building GRNs as it establishes the potential regulatory links between distal genomic regions and specific genes. By integrating TSS information with co-accessibility, the script identifies which open chromatin regions are likely to regulate which genes, providing a more direct and biologically relevant set of interactions for GRN construction.

### 7. [`7_scan_motifs_for_different_cell_types.py`](trimmed_GRN_derivation/7_scan_motifs_for_different_cell_types.py)
*   **Purpose:** This script identifies transcription factor (TF) binding motifs within the processed peak regions for different cell types. It leverages multiple motif databases to create a comprehensive set of motifs and then scans the accessible chromatin regions to find potential TF binding sites.
*   **Implementation Explanation:**
    *   **Genome Installation:** Ensures the `hg19` reference genome is installed using `genomepy`, which is necessary for motif scanning.
    *   **Motif Database Preparation:** Downloads and prepares additional motif databases (JASPAR, HOCOMOCO) and combines them with a base CisBP motif set. A custom function `parse_meme_to_motifs` is included to handle MEME format motif files.
    *   **TFinfo Object Creation:** For each cell type, it loads the `processed_peak_file_*.csv` (from script 6) and creates a `celloracle.motif_analysis.TFinfo` object, which is designed to store and manage TF-related information.
    *   **Motif Scanning:** Scans the peak regions for occurrences of the combined transcription factor binding motifs. It uses a specified False Positive Rate (FPR) threshold (e.g., 0.05) to control the stringency of motif detection.
    *   **Filtering and GRN Generation:** After scanning, it filters the motifs by a score threshold (e.g., 6) to retain high-confidence binding sites. It then generates a TF-target gene regulatory network (GRN) dataframe, which represents potential regulatory interactions.
    *   **Saving Results and Statistics:** Saves the `TFinfo` object and the generated GRN dataframe (as a Parquet file) for each cell type. It also calculates and saves summary statistics, including the number of peaks, motifs, and non-zero elements in the GRN, providing insights into the density and complexity of the inferred networks.
*   **Biological Significance:** This script is central to inferring GRNs. By identifying TF binding sites within accessible chromatin, it provides direct evidence for which TFs might regulate which genes. The use of multiple motif databases increases the coverage and accuracy of TF binding site prediction, leading to more robust and comprehensive GRNs that reflect the intricate regulatory landscape of different cell types.

### 8. [`8_celloracle_scan_motifs.py`](trimmed_GRN_derivation/8_celloracle_scan_motifs.py)
*   **Purpose:** This script serves as a verification and comparison step for the motif scanning results generated by `7_scan_motifs_for_different_cell_types.py`. It loads the generated GRN dataframes and compares their characteristics (e.g., shape, number of non-zero elements) with a base GRN, ensuring the quality and consistency of the inferred networks.
*   **Implementation Explanation:**
    *   **Configuration and Data Loading:** Sets up parameters for cell types and ages. It then loads the `.celloracle.parquet` file generated by the motif scanning script for a specific cell type.
    *   **GRN Characteristics Analysis:** Prints the shape of the loaded GRN dataframe and calculates the total number of non-zero elements (representing active regulatory links). It also shows the number of non-zero elements per row (target gene).
    *   **Comparison with Base GRN:** Loads a pre-existing "base GRN" (`2023_11_tfi.celloracle.parquet`) for comparison. It performs the same characteristic analysis on the base GRN and calculates the ratio of non-zero elements between the newly generated GRN and the base GRN.
*   **Biological Significance:** This script is crucial for quality control and validation of the inferred GRNs. By comparing the newly generated networks with a known reference, it helps ensure that the motif scanning and GRN inference processes are yielding biologically plausible and consistent results. This step is important for maintaining the integrity and reliability of the GRNs used in downstream gene set scoring.

### 9. [`9_create_scores_procedural.py`](trimmed_GRN_derivation/9_create_scores_procedural.py)
*   **Purpose:** This script is the core of the gene set scoring derivation. It integrates scRNA-seq data with the inferred GRNs to simulate gene perturbations and calculate gene set scores. It uses CellOracle to build a comprehensive regulatory model and then assesses the impact of perturbing specific genes of interest.
*   **Implementation Explanation:**
    *   **Data Loading and Preparation:** Loads the processed scRNA-seq data (`subseted_rna_andata.h5ad`) and filters it to include highly variable genes and genes of interest. It initializes a `celloracle.Oracle` object with the scRNA-seq data.
    *   **TF-TG Dictionary Enhancement:** Enhances the Transcription Factor-Target Gene (TF-TG) dictionary by incorporating information from multiple external databases (e.g., `2023_11_CellOracleProof.tsv`, `trrust_rawdata.human.tsv`, `Brain_GTEx-regulons.txt`, `Fetal-Brain-regulons.txt`, `TFLink_Homo_sapiens_interactions_*.tsv`). This step enriches the knowledge base of known regulatory interactions.
    *   **Dimensionality Reduction and Imputation:** Performs Principal Component Analysis (PCA) on the scRNA-seq data to reduce dimensionality and then applies k-nearest neighbors (KNN) imputation to fill in missing data, preparing the data for GRN inference.
    *   **GRN Inference and Link Filtering:** For each cell type, it loads the motif scan results (GRN) and imports them into the Oracle object. It then infers regulatory links between TFs and target genes using `oracle.get_links`, filters these links based on statistical significance and strength, and calculates network scores.
    *   **GRN Fitting and Simulation:** Fits the GRN for simulation, allowing the model to predict gene expression changes upon perturbation. It then simulates the perturbation (e.g., knockout) of specific genes of interest (`gois`).
    *   **Score Calculation and Saving:** Calculates simulation scores (e.g., log fold change) and GRN scores (e.g., coefficient means) for each perturbed gene and cell type. These scores are then saved to CSV files (`scores_sim_all_new.csv` and `scores_grn_all_from_comb_run_new.csv`).
    *   **Visualization (Optional):** Includes optional plotting steps to visualize PCA elbow plots, gene expression on UMAP, and simulation results.
*   **Biological Significance:** This is the central script for generating the gene set scores. By integrating multi-omics data (scRNA-seq and ATAC-seq-derived GRNs) and simulating gene perturbations, it provides a mechanistic understanding of how changes in specific genes (e.g., those related to environmental exposures or disease) propagate through the regulatory network to affect the expression of other genes. The resulting scores quantify the impact of these perturbations, which can then be used for gene set enrichment analysis and understanding biological pathways.

### 10. [`10_examin_results.py`](trimmed_GRN_derivation/10_examin_results.py)
*   **Purpose:** This script performs a comprehensive examination and analysis of the simulation and GRN scores generated by `9_create_scores_procedural.py`. It provides insights into the distribution of scores, identifies unique perturbed genes and cell types, and facilitates common analysis tasks like finding common target genes.
*   **Implementation Explanation:**
    *   **Loading Scores:** Loads the `scores_sim_all_new.csv` (simulation scores) and `scores_grn_all_from_comb_run_new.csv` (GRN scores) files.
    *   **Data Renaming and Transformation:** Renames columns for consistency (e.g., `gene` to `target_gene`, `log_fold_change` to `score` for simulation scores; `target` to `target_gene`, `coef_mean` to `score` for GRN scores).
    *   **Statistical Summary and Visualization:** Displays basic statistics (`describe()`) for both simulation and GRN scores. It generates histograms to visualize the distribution of key metrics like `score`, `fold_change`, `p-value`, and `X.logp`.
    *   **Unique Identifier Extraction:** Identifies and prints the unique perturbed genes and cell types present in the datasets.
    *   **Common Analysis Functions:** Utilizes helper functions (from `results_analysis.py`) such as `analyze_gene_perturbation` to extract data for a specific gene and `find_common_targets` to identify target genes that are commonly affected in both simulation and GRN results for a given perturbed gene.
*   **Biological Significance:** This script is essential for interpreting the results of the GRN inference and simulation. By analyzing the distributions of scores and identifying common targets, researchers can gain insights into the robustness of the predicted regulatory interactions and the overall impact of gene perturbations. This helps in validating the GRN model and identifying key genes and pathways that are most significantly affected, guiding further biological investigation.
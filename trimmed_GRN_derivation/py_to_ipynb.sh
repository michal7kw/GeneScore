#!/bin/bash

jupytext --to notebook 1_herring_scRNA_processing.py
jupytext --to notebook 2_herring_atac_processing.py
jupytext --to notebook 3_celltype_consensus_regions.py
jupytext --to notebook 4_coaccess.py
jupytext --to notebook 5_prepare_peak_data_for_celloracle.py
jupytext --to notebook 6_integrate_tss_with_coaccess_regions.py
jupytext --to notebook 7_scan_motifs_for_different_cell_types.py
jupytext --to notebook 8_celloracle_scan_motifs.py
jupytext --to notebook 9_create_scores.py
jupytext --to notebook 9_create_scores_procedural.py
jupytext --to notebook 9_create_scores_procedural_original.py
jupytext --to notebook 10_examin_results.py

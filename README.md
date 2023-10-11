# Melanoma transcriptional phenotype switching meta-analysis manuscript

Code for meta-analysis comparing melanoma cell states identified from bulk and
scRNA-seq studies in terms of cell states, gene co-expression, and predicted 
drug sensitivities.

# Overview
An overview for each of the directory structure. Analysis summarised in below
figure.

    ![title](https://github.com/BradBalderson/Melanoma_manuscript/blob/main/docs/MelanomaMethodsFlow_v2.png)

For information about individual scripts and notebooks, see docs/index.md.

scripts/

    X1_tirosh_trajectory/
        -> Code for the tirosh scRNA-seq patient tumour trajectory analysis.
        
    X2_tsoi_vs_tirosh/
        -> Code comparing bulk RNA-seq low-passage cell line states to the 
            tirosh trajectory.
            
    X3_tcga_vs_tirosh/
        -> Compares TCGA bulk RNA-seq patient samples to tirosh scRNA-seq 
            patient tumour trajectory.
        
    X4_rambow_vs_tirosh/
        -> Compares the rambow PDX scRNA-seq trajectory to the tirosh patient
            scRNA-seq trajectory.
            
    X5_drug_analysis/
        -> Performs drug susceptibility analysis for each dataset and compares
            predictions.







# PhenoID_single_cell_flow_cytometry_analysis
A work flow and set of tools to use flow cytometry to quantify cell types within a complex tissue.
Tools with workbooks.
1. R-preprocess: Read in fsc files and preprocess - select live single cells (removing doublets and debris)
2. R-align: transform, align and reverse transform multiple samples to remove batch effects.
3. R-correlation-assigment: assign cell types by correlating each cell to a reference matrix.  To be used for assisting in the cluster annotation step.
4. py-Clustering and/or R-Clustering: normalize and cluster with different methods, steps for optimization.
5. R-cluster comparision: Calculate statistics for different clustering methods and parameters to identify the best clustering for a given dataset. 
5. py-Visualization and/or R-Visualization: Visualize clusters and explore expression levesl and correlations cell assignement labels to assist in expert cell annotation. 
6. R-cell-annotation: Train and apply random forest classifier to assign cell types based on expert annotation 
7. Cell count analysis tools. 

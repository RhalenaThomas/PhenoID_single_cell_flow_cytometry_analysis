### find_correlation()
Find correlation between the test matrix and the reference matrix. 
Compare pre-processed expression matrix with the expected value for each cell type, return the best and second best correlated cell types for each sample in the expression matrix.

#### input
**test_path:** The test matrix in .csv format.

**reference_path:** The reference matrix in .csv format.

**output_path:** Where the output will be saved

**min_corr:** The minimal correlation that will be accepted. If lower than this value, the cell will be labeled with "unknown" cell type. By default, the value is 0.1


**min_diff:** The minimal difference between the best correlation and the second-best correlation. If smaller than this value, the cell will be labeled with a mixed cell type -- "best-correlated - second-best-correlated" cell type. By default, the value is 0.05


#### output
Generate csv and plots pdfs.


### stats_plot()
Plot 3 statistics that measures the quality of the clustering -- silhouette score, Calinski-Harabasz index and Davies–Bouldin index.  
This function is used by flowSOM clustering or phenograph clustering.

#### input
**stats_ls:** a dataframe generated by the clustering functions (for example: louvain_clustering), the dataframe should contain columns: kn, nc, si, ch, db -- corresponding to k-nearest neighbor, number of clusters, silhouette score, Calinski-Harabasz index and Davies–Bouldin index.

**output_path:** Where the output will be saved

**input_name:** The name of the input, for example: "flowSOMalign". It will be part of the plot names for identifying purpose.

**clust_method:** The clustering method used to generate the stats_ls. for example: "Pheno", "Louvain", "FlowSOM". It will be part of the plot names for identifying purpose.


#### output
Generate 1 pdf per stats method. Within each pdf, there are 2 plots, one with the number of clusters as x-axis, one with the k-nearist neighbour as x-axis. In both cases, the statistical result (si, ch, db) will be the y-axis. 


### louvain_stats_plotting()
Plot 3 statistics that measures the quality of the clustering -- silhouette score, Calinski-Harabasz index and Davies–Bouldin index. This function is used by louvain clustering. Louvain clustering has a different plotting function from the other two clustering methods because louvain clustering contains an extra parameter -- resolution.

#### input
**stats_ls:** A dataframe generated by the clustering functions (for example: louvain_clustering), the dataframe should contain columns: kn, resolution, nc, si, ch, db -- corresponding to k-nearest neighbor, resolution that represents the number of the communities, number of clusters, silhouette score, Calinski-Harabasz index and Davies–Bouldin index

**output_path:** Where the output will be saved

**input_name:** The name of the input, for example: "flowSOMalign". It will be part of the plot names for identifying purpose.

**clust_method:** The clustering method used to generate the stats_ls. for example: "Pheno", "Louvain", "FlowSOM". It will be part of the plot names for identifying purpose.

#### output
Generate 1 pdf per stats method. Within each pdf, there are 2 plots, one with the number of clusters as x-axis, one with the k-nearist neighbour as x-axis. In both cases, the statistical result (si, ch, db) will be the y-axis. 


### plot_comparison()
Generate plots that compare the statistical outcome across 3 clustering methods. 3 types of plots will be generated corresponding to 3 statistical methods -- silhouette score, Calinski-Harabasz index and Davies–Bouldin index. 

#### input
**input_name:** The name of the input, for example: "AlignTrans". It will be part of the plot names for identifying purpose.

**output_path:** Where the output will be saved

**louvain:** A dataframe. The stats list generated by louvain_clustering(). It should contain columns: kn, resolution, nc, si, ch, db -- corresponding to k-nearest neighbor, resolution that represents the number of the communities, number of clusters, silhouette score, Calinski-Harabasz index and Davies–Bouldin index.

**flowsom:** A dataframe. The stats list generated by flowsom_clustering(). It should contain columns: krange, nc, si, ch, db -- corresponding to range of k, number of clusters, silhouette score, Calinski-Harabasz index and Davies–Bouldin index.


**phenograph:** A dataframe. The stats list generated by phenograph_clustering(). It should contain columns: kn, nc, si, ch, db -- corresponding to k-nearest neighbor, number of clusters, silhouette score, Calinski-Harabasz index and Davies–Bouldin index

#### output 
3 pdf in total will be generated, and each pdf accounts for one of the three stats method. Within each pdf, there are 2 plots, one with the number of clusters as x-axis, one with the k-nearist neighbour as x-axis. In both cases, the statistical result (si, ch, db) will be the y-axis. On each plot, statistical outcome for all three types of clustering methods will be displayed and compared.


### louvain_clustering()
Cluster the preprocessed expression matrix with Louvain clustering method. Run 3 statistical methods that measures the quality of the clusterings -- silhouette score, Calinski-Harabasz index and Davies–Bouldin index. Generates feature plots, stats plots, and stats list.

#### input
**input_path:** The csv path for the preprocessed expression matrix

**output_path:** Where the output will be saved

**input_name:** The name of the input, for example: "phenographAlignTrans". It will be part of the plot names for identifying purpose.

**clust_method:** The clustering method used to generate the stats_ls. for example: "Louvain". It will be part of the plot names for identifying purpose.

**kn:** k-nearest neighbor

**resolutions:** The parameter is a collection of resolutions, within which each represents the number of the communities and is unique to the Louvain cluster. By default, the parameter is c(0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,1.0,1.8)

#### output
- A heatmap grouped by batches
- UMAP features plots colored by antibodies for each kn
- UMAP plots grouped by batches for each kn
- UMAP plots for each combination of kn and resolution
- Heatmaps for each each combination of kn and resolution
- Clustree plots for each kn
- Seurat object in .Rds format for each kn
- One stats list in .csv and .Rds format. 
	- Columns: kn, resolution, nc, si, ch, db -- corresponding to k-nearest neighbor, resolution that represents the number of the communities, number of clusters, silhouette score, Calinski-Harabasz index and Davies–Bouldin index. 
	- Rows: stats for all possible combinations of kn and resolution
- Output of louvain_stats_plot()
	- Generate 1 pdf for each of the 3 stats methods. Within each pdf, there are 2 plots, one with the number of clusters as x-axis, one with the k-nearist neighbour as x-axis. In both cases, the statistical result (si, ch, db) will be the y-axis. There are 6 plots in total.


### flowsom_clustering()
Cluster the preprocessed expression matrix with flowSOM clustering method. Run 3 statistical methods that measures the quality of the clusterings -- silhouette score, Calinski-Harabasz index and Davies–Bouldin index. Generates feature plots, stats plots, and stats list.

**krange:** k-nearest neighbor

**input_path:** The csv path for the preprocessed expression matrix

**output_path:** Where the output will be saved

**input_name:** The name of the input, for example: "phenographAlignTrans". It will be part of the plot names for identifying purpose.

**clust_method:** The clustering method used to generate the stats_ls. for example: "FlowSOM". It will be part of the plot names for identifying purpose.

#### output
- A heatmap grouped by batches
- A UMAP feature plot colored by antibodies
- FlowSOM clustering UMAP plot for each k range
- FlowSOM clustering heatmap for each k range
- Heatmap feature plot colored by antibodies for each k range 
- Clustree plots for each kn
- UMAP plot with cell types
- Seurat object in .Rds format for each kn
- One stats list in .csv and .Rds format. 
	- Columns: kn, resolution, nc, si, ch, db -- corresponding to k-nearest neighbor, resolution that represents the number of the communities, number of clusters, silhouette score, Calinski-Harabasz index and Davies–Bouldin index. 
	- Rows: stats for all possible combinations of kn and resolution
- Output of stats_plot()
	- Generate 1 pdf for each of the 3 stats methods. Within each pdf, there are 2 plots, one with the number of clusters as x-axis, one with the k-nearist neighbour as x-axis. In both cases, the statistical result (si, ch, db) will be the y-axis. There are 6 plots in total.


### phenograph_clustering()
Cluster the preprocessed expression matrix with phenograph clustering method. Run 3 statistical methods that measures the quality of the clusterings -- silhouette score, Calinski-Harabasz index and Davies–Bouldin index. Generates feature plots, stats plots, and stats list.

#### input
**kn:** k-nearest neighbor

**input_path:** The csv path for the preprocessed expression matrix

**output_path:** Where the output will be saved

**input_name:** The name of the input, for example: "phenographAlignTrans". It will be part of the plot names for identifying purpose.

**clust_method:** The clustering method used to generate the stats_ls. for example: "FlowSOM". It will be part of the plot names for identifying purpose.

#### output
- A heatmap grouped by batches
- A UMAP feature plot colored by antibodies
- A UMAP plot grouped by batches
- Phenograph clustering UMAP plot for each kn
- Phenograph clustering heatmap for each kn
- Clustree plots for each kn
- Seurat object in .Rds format for each kn
- One stats list in .csv and .Rds format. 
	- Columns: kn, resolution, nc, si, ch, db -- corresponding to k-nearest neighbor, resolution that represents the number of the communities, number of clusters, silhouette score, Calinski-Harabasz index and Davies–Bouldin index. 
	- Rows: stats for all possible combinations of kn and resolution
- Output of stats_plot()
	- Generate 1 pdf for each of the 3 stats methods. Within each pdf, there are 2 plots, one with the number of clusters as x-axis, one with the k-nearist neighbour as x-axis. In both cases, the statistical result (si, ch, db) will be the y-axis. There are 6 plots in total.


### clustering()
The clustering function that includes all three clustering options. 

#### input
**resolution:** The parameter is a collection of resolutions, within which each represents the number of the communities and is unique to the Louvain cluster. The default value is NULL.

**kn:** k-nearest neighbor

**input_path:** The csv path for the preprocessed expression matrix

**output_path:** Where the output will be saved

**input_name:** The name of the input, for example: "phenographAlignTrans". It will be part of the plot names for identifying purpose.

**clust_method:** The clustering method used to generate the stats_ls. for example: "Pheno", "Louvain", "FlowSOM". It will be part of the plot names for identifying purpose.

#### output
Depending on the clust_method, different outputs will be generated. For details, see documentations for louvain_clustering(), flowsom_clustering(), or phenograph_clustering()


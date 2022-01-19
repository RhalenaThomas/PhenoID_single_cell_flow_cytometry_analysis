# make the organoid object with meta data

mat <- readRDS("org_expression_mat.rds")
meta <- read.table("meta.tsv", header=T, sep="\t", as.is=T, row.names=1)
genes = mat[,1][[1]]
genes = gsub(".+[|]", "", genes)
mat2 = data.frame(mat[,-1], row.names=genes)
UCSCorg <- CreateSeuratObject(counts = mat, project = "organoids", meta.data=meta)
saveRDS(UCSCorg,"UCSCorg_metadata.rds")


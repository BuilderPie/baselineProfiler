path = "../../../static/data/CCLE/"
outputPath = "../../../static/data_computed/CCLE/"

#############################################################################
# ======= step 1. CCLE_sample_info preprocess
# ======= find unique sample_info lineages

CCLE_sample_info = file.path(path, "sample_info.csv")
dir.create(file.path(outputPath, "sample_info"), recursive = T, showWarnings = F)

sample_info = read.csv(CCLE_sample_info)
sample_info_lineage = data.frame(lineage = sort(unique(sample_info$lineage)))

write.csv(sample_info_lineage, 
          file = file.path(outputPath, "sample_info", "sample_info_lineage.csv"),
          row.names = F
)

sample_info_trim = sample_info[, c("DepMap_ID","cell_line_name", "stripped_cell_line_name", "lineage")]

write.csv(sample_info_trim, 
          file = file.path(outputPath, "sample_info", "sample_info.csv"),
          row.names = F
)
#############################################################################
# ======= step 2. CCLE_expression data preprocess
CCLE_expression = file.path(path, "CCLE_expression.csv")


df = read.csv(CCLE_expression, row.names = 1)
df = t(df)
# ======= remove the Entrez id part of the rowname
rownames_tmp = rownames(df)
rownames_tmp = sapply(rownames_tmp, FUN = function(x) strsplit(x, "\\.\\.")[[1]][1])
rownames_tmp = unname(rownames_tmp)
# ======= check if there are duplicated rownames
which(duplicated(rownames_tmp))

# ======= update the rownames
rownames(df) = rownames_tmp
# ======= we need to keep dot in gene name without hyphen
# ======= because python can not access labels with hyphen
# dot_index = grep("\\.", rownames_tmp)
# gsub("\\.", "-", rownames_tmp[dot_index])



# rownames(df)[grep("HLA", rownames(df))]
# rownames_tmp[grep("HLA", rownames_tmp)]

# ======= export three objects: df, rownames / gene names, colnames / celllines
dir.create(file.path(outputPath, "expression"), recursive = T, showWarnings = F)

df_rowname = data.frame(gene = rownames(df))
df_colname = data.frame(cell_line = colnames(df))

write.csv(df_rowname, 
          file = file.path(outputPath, "expression", "expression_rownames.csv"), 
          row.names = F)
write.csv(df_colname, 
          file = file.path(outputPath, "expression", "expression_colnames.csv"),
          row.names = F)

write.csv(df, 
          file = file.path(outputPath, "expression", "expression.csv"),
          )

#############################################################################
# ======= step 3. CCLE_CRISPR preprocess

CCLE_CRISPR = file.path(path, "CRISPR_gene_effect.csv")


df = read.csv(CCLE_CRISPR, row.names = 1)
df = t(df)

df[1:3,1:3]
# ======= remove the Entrez id part of the rowname
rownames_tmp = rownames(df)
rownames_tmp = sapply(rownames_tmp, FUN = function(x) strsplit(x, "\\.\\.")[[1]][1])
rownames_tmp = unname(rownames_tmp)
# ======= check if there are duplicated rownames
which(duplicated(rownames_tmp))
# ======= update the rownames
rownames(df) = rownames_tmp

# ======= export three objects: df, rownames / gene names, colnames / celllines
dir.create(file.path(outputPath, "CRISPR", "Broad"), recursive = T, showWarnings = F)

df_rowname = data.frame(gene = rownames(df))
df_colname = data.frame(cell_line = colnames(df))

write.csv(df_rowname, 
          file = file.path(outputPath, "CRISPR", "Broad", "gene_effect_rownames.csv"), 
          row.names = F)
write.csv(df_colname, 
          file = file.path(outputPath, "CRISPR", "Broad", "gene_effect_colnames.csv"),
          row.names = F)

write.csv(df, 
          file = file.path(outputPath, "CRISPR", "Broad", "gene_effect.csv"),
)

#############################################################################
# ======= step 4. CCLE_proteomics preprocess
CCLE_proteomics = file.path(path, "protein_quant_current_normalized.csv")
dir.create(file.path(outputPath, "proteomics"), recursive = T, showWarnings = F)

sampleInfo = file.path(path, "Table_S1_Sample_Information.xlsx")


# ======== load expression data ======================= #
library(readr)

expr <- read_csv(CCLE_proteomics)
expr <- as.data.frame(expr)

colnames(expr) = sapply(colnames(expr), FUN = function(x) strsplit(x, "[..]")[[1]][1])

# ======== subset expr to keep only Gene_Symbol and effecive cell line columns
df = expr[, c(2, 49:dim(expr)[2])]
df_colnames = colnames(df)
df_colnames = unname(sapply(df_colnames, FUN = function(x){
  strsplit(x, split = "_")[[1]][1]
}))
# ====== find duplicated cell lines and take average
dup_col_ind = which(duplicated(df_colnames))
dup_col_ind_all = sapply(dup_col_ind, FUN = function(x){
  which(df_colnames == df_colnames[x])
})

merged_dup = lapply(dup_col_ind, FUN = function(x){
  ind = which(df_colnames == df_colnames[x])
  tmp = df[, ind]
  tmp = data.frame(x = apply(tmp, 1,mean,na.rm=T))
  colnames(tmp) = df_colnames[x]
  tmp
})

merged_dup = (do.call(cbind, merged_dup))

df_rm_dup = df[, -as.vector(dup_col_ind_all)]

df_new = cbind(df_rm_dup, merged_dup)
colnames(df_new) = unname(sapply(colnames(df_new), FUN = function(x){
  strsplit(x, split = "_")[[1]][1]
}))
# ======= find duplicated gene names and rename with number
dup_gene_ind = which(duplicated(df_new$Gene))
dup_gene = unique(df_new$Gene[dup_gene_ind])

raw_Gene = df_new$Gene
sapply(dup_gene, FUN = function(x){
  ind = which(df_new$Gene == x)
  raw_Gene[ind] <<- paste0(x, "_", 1:length(ind))
  print(raw_Gene[ind])
})

df_new$Gene = raw_Gene

# ======= remove Gene names with NA
na_ind = which(is.na(df_new$Gene))
df_new = df_new[-na_ind, ]

# === test if there are any duplicated Gene
which(duplicated(df_new$Gene))

# ===== move Gene to rownames and write to csv file
df_trim = df_new[, -1]
rownames(df_trim) = df_new$Gene

write.csv(df_trim, 
          file = file.path(outputPath, "proteomics", "protein_normalized.csv"),
)

df_rowname = data.frame(gene = rownames(df_trim))
df_colname = data.frame(cell_line = colnames(df_trim))

write.csv(df_rowname, 
          file = file.path(outputPath, "proteomics", "protein_normalized_rownames.csv"), 
          row.names = F)
write.csv(df_colname, 
          file = file.path(outputPath, "proteomics", "protein_normalized_colnames.csv"),
          row.names = F)

# ======== load sample_info.csv ======================= #
library(openxlsx)
sampleInfo = read.xlsx(sampleInfo, sheet = "Sample_Information")

# ======== extract lineage info
sample_info_lineage = data.frame(lineage = sort(unique(sampleInfo$Tissue.of.Origin)))

write.csv(sample_info_lineage,
          file = file.path(outputPath, "sample_info", "sample_info_lineage_proteomics.csv"),
          row.names = F
)

sampleInfo_trim = sampleInfo[, c("Cell.Line", "Tissue.of.Origin")]
colnames(sampleInfo_trim) = c("cell_line", "lineage")
sampleInfo_trim$cell_line_CCLE = sapply(sampleInfo$CCLE.Code, function(x){
  strsplit(x, "_")[[1]][1]
})

write.csv(sampleInfo_trim,
          file = file.path(outputPath, "sample_info", "sample_info_proteomics.csv"),
          row.names = F
)

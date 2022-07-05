path = "../../../static/data/HPA/"
outputPath = "../../../static/data_computed/HPA/"

CCLEmetaPath = "../../../static/data/CCLE/"

dir.create(file.path(outputPath, "Exprsn_cell_line"), recursive = T, showWarnings = F)
dir.create(file.path(outputPath, "Exprsn_blood_cell"), recursive = T, showWarnings = F)


#############################################################################
# ======= step 1. HPA_sample_info preprocess
# ======= find unique sample_info lineages

CCLE_sample_info = file.path(CCLEmetaPath, "sample_info.csv")
dir.create(file.path(outputPath, "sample_info"), recursive = T, showWarnings = F)

sample_info = read.csv(CCLE_sample_info)

# ======= read in HPA rna_celline.tsv file to match cell line with corresponding CCLE lineages
HPA_cellline = file.path(path, "rna_celline.tsv")
df = read.csv(HPA_cellline, sep = "\t")

cell_line = unique(df$Cell.line)

HPA_meta = data.frame(HPA_cell_line = cell_line)
HPA_meta$CCLE_cell_line_name = rep("", length(cell_line))
HPA_meta$CCLE_stripped_cell_line_name = rep("", length(cell_line))
HPA_meta$lineage = rep("", length(cell_line))
print(head(HPA_meta))

sapply(1:dim(HPA_meta)[1], FUN = function(x){
  tmp = HPA_meta$HPA_cell_line[x]
  idx1 = grep(paste0("^",tmp,"$"), sample_info$cell_line_name, ignore.case = T)
  idx2 = grep(paste0("^",tmp,"$"), sample_info$stripped_cell_line_name, ignore.case = T)
  
  if (length(idx1) >0 | length(idx2)>0){
    print(x)
    idx = c(idx1[1], idx2[1])
    idx = idx[which(!is.na(idx))[1]]
    HPA_meta$CCLE_cell_line_name[x] <<- sample_info$cell_line_name[idx]
    # print(sample_info$cell_line_name[idx1[1]])
    HPA_meta$CCLE_stripped_cell_line_name[x] <<- sample_info$stripped_cell_line_name[idx]
    # print(sample_info$stripped_cell_line_name[idx2[1]])
    HPA_meta$lineage[x] <<- sample_info$lineage[idx]
    # print(sample_info$lineage[idx1[1]])
  }
})

write.csv(HPA_meta, file = file.path(outputPath, "sample_info", "HPA_meta_tmp.csv"), row.names = F)

# ================================== #
sample_info_lineage = data.frame(lineage = sort(unique(HPA_meta$lineage)))

write.csv(sample_info_lineage, 
          file = file.path(outputPath, "sample_info", "sample_info_lineage.csv"),
          row.names = F
)
# ================================== #
# load annotated sample_info table
HPA_meta = read.csv(file = outputPath, "sample_info", "HPA_meta.csv")
#############################################################################
# ======= step 2. HPA RNA cell lines preprocess
# pivot long 
library(tidyr)
library(dplyr)

hpa_expData_nTPM <- df[,c("Gene.name","Cell.line","nTPM")] %>%
  # dplyr::distinct(., Cell.line)
  dplyr::group_by(Cell.line, Gene.name) %>%
  # dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  # dplyr::filter(dplyr::n() > 1L) %>%
  slice_sample(n = 1) %>% 
  pivot_wider(.,names_from = Cell.line, values_from = nTPM)



write.csv(x = hpa_expData_nTPM, file = file.path(outputPath, "Exprsn_cell_line", "rna_celline_nTPM.csv"), row.names = F)

# ======== load nTPM expData ========================================== #
hpa_expData_nTPM = read.csv(file.path(outputPath, "Exprsn_cell_line", "rna_celline_nTPM.csv"), check.names = F)
rownames(hpa_expData_nTPM) = hpa_expData_nTPM$Gene.name
hpa_expData_nTPM = hpa_expData_nTPM[, -1]

# ======== convert to log2 nTPM ======================================= #
hpa_expData_log2_nTPM = hpa_expData_nTPM
hpa_expData_log2_nTPM[!is.na(hpa_expData_log2_nTPM)] = log2(1 + hpa_expData_log2_nTPM[!is.na(hpa_expData_log2_nTPM)])

write.csv(x = hpa_expData_log2_nTPM, file = file.path(outputPath, "Exprsn_cell_line", "rna_celline_log2_nTPM.csv"), row.names = T)


# output rowname and colname
df_rowname = data.frame(gene = rownames(hpa_expData_log2_nTPM))
df_colname = data.frame(cell_line = colnames(hpa_expData_log2_nTPM))

write.csv(df_rowname, 
          file = file.path(outputPath, "Exprsn_cell_line", "rna_celline_rownames.csv"), 
          row.names = F)
write.csv(df_colname, 
          file = file.path(outputPath, "Exprsn_cell_line", "rna_celline_colnames.csv"),
          row.names = F)

# ========= process TPM matrix ======================================== #
hpa_expData_TPM <- df[,c("Gene.name","Cell.line","TPM")] %>%
  # dplyr::distinct(., Cell.line)
  dplyr::group_by(Cell.line, Gene.name) %>%
  # dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  # dplyr::filter(dplyr::n() > 1L) %>%
  slice_sample(n = 1) %>% 
  pivot_wider(.,names_from = Cell.line, values_from = TPM)

write.csv(x = hpa_expData_TPM, file = file.path(outputPath, "Exprsn_cell_line", "rna_celline_TPM.csv"), row.names = F)

# ======== load TPM expData ========================================== #
hpa_expData_TPM = read.csv(file.path(outputPath, "Exprsn_cell_line", "rna_celline_TPM.csv"), check.names = F)
rownames(hpa_expData_TPM) = hpa_expData_TPM$Gene.name
hpa_expData_TPM = hpa_expData_TPM[, -1]

# ======== convert to log2 TPM ======================================= #
hpa_expData_log2_TPM = hpa_expData_TPM
hpa_expData_log2_TPM[!is.na(hpa_expData_log2_TPM)] = log2(1 + hpa_expData_log2_TPM[!is.na(hpa_expData_log2_TPM)])

write.csv(x = hpa_expData_log2_TPM, file = file.path(outputPath, "Exprsn_cell_line", "rna_celline_log2_TPM.csv"), row.names = T)

# ======== end of conversion expression data to matrix ================ #
####################################################################################

#############################################################################
# ======= step 3. HPA RNA blood cells preprocess
HPA_blood_cell = file.path(path, "rna_blood_cell.tsv")
HPA_expData_blood_cell = read.csv(HPA_blood_cell, sep = "\t")

unique(HPA_expData_blood_cell$Blood.cell)

# ========== pivot nTPM data table to wider
hpa_blood_nTPM <- HPA_expData_blood_cell[,c("Gene.name","Blood.cell","nTPM")] %>%
  # dplyr::distinct(., Cell.line)
  dplyr::group_by(Blood.cell, Gene.name) %>%
  # dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  # dplyr::filter(dplyr::n() > 1L) %>%
  slice_sample(n = 1) %>% 
  pivot_wider(.,names_from = Blood.cell, values_from = nTPM)

write.csv(x = hpa_blood_nTPM, file = file.path(outputPath, "Exprsn_blood_cell", "rna_blood_nTPM.csv"), row.names = F)

# ======== reload nTPM expData ========================================== #
hpa_blood_nTPM = read.csv(file.path(outputPath, "Exprsn_blood_cell", "rna_blood_nTPM.csv"), check.names = F)
rownames(hpa_blood_nTPM) = hpa_blood_nTPM$Gene.name
hpa_blood_nTPM = hpa_blood_nTPM[, -1]

# ======== convert to log2 nTPM ======================================= #
hpa_blood_log2_nTPM = hpa_blood_nTPM
hpa_blood_log2_nTPM[!is.na(hpa_blood_log2_nTPM)] = log2(1 + hpa_blood_log2_nTPM[!is.na(hpa_blood_log2_nTPM)])

write.csv(x = hpa_blood_log2_nTPM, file = file.path(outputPath, "Exprsn_blood_cell", "hpa_blood_log2_nTPM.csv"), row.names = T)


# output rowname and colname
df_rowname = data.frame(gene = rownames(hpa_blood_log2_nTPM))
df_colname = data.frame(cell_line = colnames(hpa_blood_log2_nTPM))

write.csv(df_rowname, 
          file = file.path(outputPath, "Exprsn_blood_cell", "rna_blood_rownames.csv"), 
          row.names = F)
write.csv(df_colname, 
          file = file.path(outputPath, "Exprsn_blood_cell", "rna_blood_colnames.csv"),
          row.names = F)

# ========== pivot TPM data table to wider
hpa_blood_TPM <- HPA_expData_blood_cell[,c("Gene.name","Blood.cell","TPM")] %>%
  # dplyr::distinct(., Cell.line)
  dplyr::group_by(Blood.cell, Gene.name) %>%
  # dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  # dplyr::filter(dplyr::n() > 1L) %>%
  slice_sample(n = 1) %>% 
  pivot_wider(.,names_from = Blood.cell, values_from = TPM)

write.csv(x = hpa_blood_TPM, file = file.path(outputPath, "Exprsn_blood_cell", "rna_blood_TPM.csv"), row.names = F)

# ======== reload TPM expData ========================================== #
hpa_blood_TPM = read.csv(file.path(outputPath, "Exprsn_blood_cell", "rna_blood_TPM.csv"), check.names = F)
rownames(hpa_blood_TPM) = hpa_blood_TPM$Gene.name
hpa_blood_TPM = hpa_blood_TPM[, -1]

# ======== convert to log2 TPM ======================================= #
hpa_blood_log2_TPM = hpa_blood_TPM
hpa_blood_log2_TPM[!is.na(hpa_blood_log2_TPM)] = log2(1 + hpa_blood_log2_TPM[!is.na(hpa_blood_log2_TPM)])

write.csv(x = hpa_blood_log2_TPM, file = file.path(outputPath, "Exprsn_blood_cell", "hpa_blood_log2_TPM.csv"), row.names = T)

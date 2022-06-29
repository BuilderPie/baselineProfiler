gtexPath = "../../../static/data/GTEx/"
outputPath = "../../../static/data_computed/GTEx/"

#############################################################################
# ======= step 1. GTEx_sample_info preprocess
# ======= find unique sample_info lineages

SampleAttributesDS = read.csv(file.path(gtexPath, "GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"), sep = "\t")

lineage = unique(SampleAttributesDS$SMTS)

dir.create(file.path(outputPath, "sample_info"), recursive = T, showWarnings = F)
write.csv(data.frame(lineage = sort(lineage)), 
          file = file.path(outputPath, "sample_info", "sample_info_lineage.csv"), 
          row.names = F)

sample_info = SampleAttributesDS[, c("SAMPID", "SMTS")]
write.csv(sample_info, 
          file = file.path(outputPath, "sample_info", "sample_info.csv"), 
          row.names = F)

# ==== export both row name / gene name and col name /sample name
df = read.csv(file.path(gtexPath, "data_subset", "Kidney.csv"))

gene_name = df$gene
which(duplicated(gene_name))

write.csv(data.frame(gene = gene_name), 
          file = file.path(outputPath, "sample_info", "gene_name.csv"), 
          row.names = F)
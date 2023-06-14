library(biomaRt)
ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")

gene_ids <- read.csv("../../Dataset/mart_export.txt", header = TRUE, stringsAsFactors = FALSE)

# Extract the Ensembl IDs into a vector
ensembl_ids <- gene_ids$Gene.stable.ID
refseq_ids = getBM(attributes = c("ensembl_gene_id", "refseq_mrna"),
                   filters = "ensembl_gene_id",
                   values = ensembl_ids,
                   mart = ensembl)
print(refseq_ids)
refseq_ids_uniq <- unique(refseq_ids$refseq_mrna)
print(length(refseq_ids_uniq))


library("plyr")
library("dplyr")
library("tidyr")
library("devtools")
library("GenomicFeatures")
load_all("../macrophage-gxe-study/macrophage-gxe-study/seqUtils/")

#Load transcript annotations from disk
txdb = loadDb("annotations/geuvadis/TranscriptDb_GRCh37_050815.db")
exons = exonsBy(txdb, by = "tx", use.names = TRUE)
cds = exonsBy(txdb, by = "tx", use.names = TRUE)

#Load Gencode BASIC transcripts
filtered_metadata = readRDS("annotations/geuvadis/Homo_sapiens.GRCh37.050815.transcript_data.filtered.rds")
gencode_basic_transcript = dplyr::select(filtered_metadata, ensembl_gene_id, ensembl_transcript_id) %>% 
  dplyr::rename(gene_id = ensembl_gene_id, transcript_id = ensembl_transcript_id)

#Convert transcript annotations to 
gene_exon_start_end = extractExonsStartEnd(exons, gencode_basic_transcript)
gene_cdss_start_end = extractExonsStartEnd(cds, gencode_basic_transcript)
write.table(gene_exon_start_end, "annotations/geuvadis/Homo_sapiens.GRCh37.050815.gene_exon_start_end.txt", row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
write.table(gene_cdss_start_end, "annotations/geuvadis/Homo_sapiens.GRCh37.050815.gene_cdss_start_end.txt", row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)

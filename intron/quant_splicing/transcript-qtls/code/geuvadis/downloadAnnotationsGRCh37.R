library("GenomicFeatures")
library("biomaRt")
library("dplyr")

#Make transcript database
txdb_grch37 = makeTranscriptDbFromBiomart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host="grch37.ensembl.org")
saveDb(txdb_grch37, "annotations/TranscriptDb_GRCh37_050815.db")
exons = exonsBy(txdb_grch37, by = "tx", use.names = TRUE)

ensembl_mart = useMart("ENSEMBL_MART_ENSEMBL", host = "grch37.ensembl.org")
ensembl_dataset = useDataset("hsapiens_gene_ensembl",mart=ensembl_mart)

attributes = listAttributes(ensembl_dataset)

biomart_attributes = c("ensembl_gene_id", "external_gene_name", "gene_biotype","status", 
                       "chromosome_name","strand", "start_position", "end_position", 
                       "transcript_start","transcript_end","ensembl_transcript_id", 
                       "transcript_status",
                       "transcript_gencode_basic", "percentage_gc_content","external_transcript_name", 
                       "transcript_length", "transcript_biotype", "ccds")
transcript_data = getBM(attributes = biomart_attributes, mart = ensembl_dataset) %>% tbl_df()
saveRDS(transcript_data, "annotations/geuvadis/Homo_sapiens.GRCh37.050815.transcript_data.rds")

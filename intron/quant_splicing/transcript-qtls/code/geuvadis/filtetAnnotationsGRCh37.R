library("dplyr")

#Load transcript metadata
valid_chromosomes = c("1","10","11","12","13","14","15","16","17","18","19",
                      "2","20","21","22","3","4","5","6","7","8","9","MT","X","Y")
valid_gene_biotypes = c("lincRNA","protein_coding","IG_C_gene","IG_D_gene","IG_J_gene",
                        "IG_V_gene", "TR_C_gene","TR_D_gene","TR_J_gene", "TR_V_gene",
                        "3prime_overlapping_ncrna","known_ncrna", "processed_transcript",
                        "antisense","sense_intronic","sense_overlapping")

#Import annotations
transcript_metadata = readRDS("annotations/Homo_sapiens.GRCh37.050815.transcript_data.rds")

#Filter annotations
filtered_metadata = dplyr::filter(transcript_metadata, transcript_gencode_basic == "GENCODE basic") %>%
  dplyr::filter(chromosome_name %in% valid_chromosomes, gene_biotype %in% valid_gene_biotypes) %>%
  dplyr::filter(transcript_biotype != "retained_intron", transcript_biotype != "nonsense_mediated_decay")
saveRDS(filtered_metadata, "annotations/geuvadis/Homo_sapiens.GRCh37.050815.transcript_data.filtered.rds")

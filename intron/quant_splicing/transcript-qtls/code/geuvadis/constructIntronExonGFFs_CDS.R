library("plyr")
library("dplyr")
library("devtools")
library("rtracklayer")
load_all("../macrophage-gxe-study/macrophage-gxe-study/seqUtils/")

#Import exon start-end coordinates
exon_start_end = read.table("annotations/geuvadis/Homo_sapiens.GRCh37.050815.gene_cdss_start_end.txt", stringsAsFactors = FALSE)
colnames(exon_start_end) = c("gene_id", "exon_starts", "exon_ends")

#Import filtered metadata
filtered_metadata = readRDS("annotations/geuvadis/Homo_sapiens.GRCh37.050815.transcript_data.filtered.rds")
gene_metadata = dplyr::select(filtered_metadata, ensembl_gene_id, gene_biotype, chromosome_name, strand) %>% 
  unique() %>%
  dplyr::rename(gene_id = ensembl_gene_id, chr = chromosome_name)

#Compile gene data into a data frame
gene_data = dplyr::left_join(exon_start_end, gene_metadata, by = "gene_id")

#Make a list of gene ids
gene_id_list = as.list(gene_data$gene_id)
names(gene_id_list) = gene_data$gene_id

#Make gff file of intron coorinates
intron_df_list = lapply(gene_id_list, constructIntronExonDf, gene_data, intron_gap = 10, type = "intron")
intron_df_list = intron_df_list[!unlist(lapply(intron_df_list, function(x){is.null(x)}))]#remove nulls
intron_df = ldply(intron_df_list) %>% dplyr::select(-.id) %>% dplyr::filter(end - start > 0)
saveRDS(intron_df, "annotations/geuvadis/Homo_sapiens.GRCh38.79.cds_intron_df.rds")
intron_gr = dataFrameToGRanges(intron_df)
export.gff3(intron_gr, "annotations/geuvadis/Homo_sapiens.GRCh38.79.cds_introns.gff3")

#Make gff file of exon coordinates
exon_df_list = lapply(gene_id_list, constructIntronExonDf, gene_data, intron_gap = 10, type = "exon")
intron_genes = unique(intron_df$gene_id)
exon_df = ldply(exon_df_list) %>% dplyr::select(-.id) %>% dplyr::filter(gene_id %in% intron_genes)
saveRDS(exon_df, "annotations/geuvadis/Homo_sapiens.GRCh38.79.cds_exon_df.rds")
exon_gr = dataFrameToGRanges(exon_df)
export.gff3(exon_gr, "annotations/geuvadis/Homo_sapiens.GRCh38.79.cds_exons.gff3")


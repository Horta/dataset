library("GenomicFeatures")
library("plyr")
library("dplyr")
library("devtools")
library("ggplot2")
load_all("../reviseAnnotations/")
load_all("../wiggleplotr/")

#Read transcript annotations
annotations = readRDS("annotations/Homo_sapiens.GRCh38.78.marked_data.rds")
#txdb78 = loadDb("annotations/TranscriptDb_GRCh38_78.db", "TranscriptDb", "GenomicFeatures")
txdb78 = loadDb("annotations/TranscriptDb_GRCh38_78.db")
exons = exonsBy(txdb78, by = "tx", use.names = TRUE)
cdss = cdsBy(txdb78, by = "tx", use.names = TRUE)

#Exptract annotations for plotting
plotting_annotations = select(annotations, ensembl_transcript_id, ensembl_gene_id, external_gene_name, strand) %>% 
  dplyr::rename(transcript_id = ensembl_transcript_id, gene_id = ensembl_gene_id, gene_name = external_gene_name) %>% 
  dplyr::mutate(strand = ifelse(strand == 1, "+","-"))

### NCOA7
#Extend transcripts for one gene
current_gene_id = "ENSG00000111912"
gene_data = dplyr::filter(annotations, ensembl_gene_id == current_gene_id)
tx_ids = dplyr::select(gene_data, ensembl_transcript_id) %>% unlist()
gene_exons = exons[tx_ids]
gene_cdss = cdss[tx_ids]
new_tx = extendTranscriptsPerGene(gene_data, gene_exons, gene_cdss)
new_exons = removeMetadata(gene_exons)
new_cdss = removeMetadata(gene_cdss)
new_exons[names(new_tx$exons)] = new_tx$exons
new_cdss[names(new_tx$cdss)] = new_tx$cdss

plotTranscripts(gene_exons, gene_cdss, plotting_annotations, rescale_introns = TRUE)
plotTranscripts(new_exons, new_cdss, plotting_annotations, rescale_introns = TRUE)

#Construct alternative events
events = constructAlternativeEvents("ENST00000392477","ENST00000438495",new_exons)
gene_dat = dplyr::select(plotting_annotations, gene_id, gene_name, strand) %>% 
  dplyr::filter(gene_id == current_gene_id) %>% unique()
tx_dat = data_frame(transcript_id = names(events), gene_id = current_gene_id) %>%
  left_join(gene_dat, by = "gene_id")
plotTranscripts(events, events, tx_dat)


### UTRN ###
gene_data = dplyr::filter(annotations, ensembl_gene_id == "ENSG00000152818")
tx_ids = dplyr::select(gene_data, ensembl_transcript_id) %>% unlist()
gene_exons = exons[tx_ids]
gene_cdss = cdss[tx_ids]
new_tx = extendTranscriptsPerGene(gene_data, gene_exons, gene_cdss)
new_exons = removeMetadata(gene_exons)
new_cdss = removeMetadata(gene_cdss)
new_exons[names(new_tx$exons)] = new_tx$exons
new_cdss[names(new_tx$cdss)] = new_tx$cdss

plotTranscripts(gene_exons, gene_cdss, plotting_annotations, rescale_introns = TRUE)
plotTranscripts(new_exons, new_cdss, plotting_annotations, rescale_introns = TRUE)
ggsave("UTRN_old.pdf", plot = old)
ggsave("UTRN_new.pdf", plot = new)



### SUN2 ###
gene_data = dplyr::filter(annotations, ensembl_gene_id == "ENSG00000100242")
tx_ids = dplyr::select(gene_data, ensembl_transcript_id) %>% unlist()
gene_exons = exons[tx_ids]
gene_cdss = cdss[tx_ids]
new_tx = extendTranscriptsPerGene(gene_data, gene_exons, gene_cdss)
new_exons = removeMetadata(gene_exons)
new_cdss = removeMetadata(gene_cdss)
new_exons[names(new_tx$exons)] = new_tx$exons
new_cdss[names(new_tx$cdss)] = new_tx$cdss

plotTranscripts(gene_exons, gene_cdss, plotting_annotations, rescale_introns = TRUE)
plotTranscripts(new_exons, new_cdss, plotting_annotations, rescale_introns = TRUE)

### FN1 ###
gene_data = dplyr::filter(annotations, ensembl_gene_id == "ENSG00000115414")
tx_ids = dplyr::select(gene_data, ensembl_transcript_id) %>% unlist()
gene_exons = exons[tx_ids]
gene_cdss = cdss[tx_ids]
new_tx = extendTranscriptsPerGene(gene_data, gene_exons, gene_cdss)
new_exons = removeMetadata(gene_exons)
new_cdss = removeMetadata(gene_cdss)
new_exons[names(new_tx$exons)] = new_tx$exons
new_cdss[names(new_tx$cdss)] = new_tx$cdss

plotTranscripts(gene_exons, gene_cdss, plotting_annotations, rescale_introns = TRUE)
plotTranscripts(new_exons, new_cdss, plotting_annotations, rescale_introns = TRUE)

#TMEM39A
gene_data = dplyr::filter(annotations, ensembl_gene_id == "ENSG00000176142")
tx_ids = dplyr::select(gene_data, ensembl_transcript_id) %>% unlist()
gene_exons = exons[tx_ids]
gene_cdss = cdss[tx_ids]
new_tx = extendTranscriptsPerGene(gene_data, gene_exons, gene_cdss)
new_exons = removeMetadata(gene_exons)
new_cdss = removeMetadata(gene_cdss)
new_exons[names(new_tx$exons)] = new_tx$exons
new_cdss[names(new_tx$cdss)] = new_tx$cdss

plotTranscripts(gene_exons, gene_cdss, plotting_annotations, rescale_introns = TRUE)
plotTranscripts(new_exons, new_cdss, plotting_annotations, rescale_introns = TRUE)


#TACC1
gene_data = dplyr::filter(annotations, ensembl_gene_id == "ENSG00000147526")
tx_ids = dplyr::select(gene_data, ensembl_transcript_id) %>% unlist()
gene_exons = exons[tx_ids]
gene_cdss = cdss[tx_ids]
new_tx = extendTranscriptsPerGene(gene_data, gene_exons, gene_cdss)
new_exons = removeMetadata(gene_exons)
new_cdss = removeMetadata(gene_cdss)
new_exons[names(new_tx$exons)] = new_tx$exons
new_cdss[names(new_tx$cdss)] = new_tx$cdss

plotTranscripts(gene_exons, gene_cdss, plotting_annotations, rescale_introns = TRUE)
plotTranscripts(new_exons, new_cdss, plotting_annotations, rescale_introns = TRUE)


#CLSTN1
gene_data = dplyr::filter(annotations, ensembl_gene_id == "ENSG00000171603")
tx_ids = dplyr::select(gene_data, ensembl_transcript_id) %>% unlist()
gene_exons = exons[tx_ids]
gene_cdss = cdss[tx_ids]
new_tx = extendTranscriptsPerGene(gene_data, gene_exons, gene_cdss)
new_exons = removeMetadata(gene_exons)
new_cdss = removeMetadata(gene_cdss)
new_exons[names(new_tx$exons)] = new_tx$exons
new_cdss[names(new_tx$cdss)] = new_tx$cdss

plotTranscripts(gene_exons, gene_cdss, plotting_annotations, rescale_introns = TRUE)
plotTranscripts(new_exons, new_cdss, plotting_annotations, rescale_introns = TRUE)

#KLHL36
gene_data = dplyr::filter(annotations, ensembl_gene_id == "ENSG00000135686")
tx_ids = dplyr::select(gene_data, ensembl_transcript_id) %>% unlist()
gene_exons = exons[tx_ids]
gene_cdss = cdss[tx_ids]
new_tx = extendTranscriptsPerGene(gene_data, gene_exons, gene_cdss)
new_exons = removeMetadata(gene_exons)
new_cdss = removeMetadata(gene_cdss)
new_exons[names(new_tx$exons)] = new_tx$exons
new_cdss[names(new_tx$cdss)] = new_tx$cdss

plotTranscripts(gene_exons, gene_cdss, plotting_annotations, rescale_introns = TRUE)
plotTranscripts(new_exons, new_cdss, plotting_annotations, rescale_introns = TRUE)

#FRMD4B
gene_data = dplyr::filter(annotations, ensembl_gene_id == "ENSG00000114541")
tx_ids = dplyr::select(gene_data, ensembl_transcript_id) %>% unlist()
gene_exons = exons[tx_ids]
gene_cdss = cdss[tx_ids]
new_tx = extendTranscriptsPerGene(gene_data, gene_exons, gene_cdss)
new_exons = removeMetadata(gene_exons)
new_cdss = removeMetadata(gene_cdss)
new_exons[names(new_tx$exons)] = new_tx$exons
new_cdss[names(new_tx$cdss)] = new_tx$cdss

plotTranscripts(gene_exons, gene_cdss, plotting_annotations, rescale_introns = TRUE)
plotTranscripts(new_exons, new_cdss, plotting_annotations, rescale_introns = TRUE)



#OSBPL9
gene_data = dplyr::filter(annotations, ensembl_gene_id == "ENSG00000117859")
tx_ids = dplyr::select(gene_data, ensembl_transcript_id) %>% unlist()
gene_exons = exons[tx_ids]
gene_cdss = cdss[tx_ids]
new_tx = extendTranscriptsPerGene(gene_data, gene_exons, gene_cdss)
new_exons = removeMetadata(gene_exons)
new_cdss = removeMetadata(gene_cdss)
new_exons[names(new_tx$exons)] = new_tx$exons
new_cdss[names(new_tx$cdss)] = new_tx$cdss

plotTranscripts(gene_exons, gene_cdss, plotting_annotations, rescale_introns = TRUE)
plotTranscripts(new_exons, new_cdss, plotting_annotations, rescale_introns = TRUE)


#OSBPL1A
gene_data = dplyr::filter(annotations, ensembl_gene_id == "ENSG00000141447")
tx_ids = dplyr::select(gene_data, ensembl_transcript_id) %>% unlist()
gene_exons = exons[tx_ids]
gene_cdss = cdss[tx_ids]
new_tx = extendTranscriptsPerGene(gene_data, gene_exons, gene_cdss)
new_exons = removeMetadata(gene_exons)
new_cdss = removeMetadata(gene_cdss)
new_exons[names(new_tx$exons)] = new_tx$exons
new_cdss[names(new_tx$cdss)] = new_tx$cdss

plotTranscripts(gene_exons, gene_cdss, plotting_annotations, rescale_introns = TRUE)
plotTranscripts(new_exons, new_cdss, plotting_annotations, rescale_introns = TRUE)

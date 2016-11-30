#Download transcript annotations from biomart
Rscript code/geuvadis/downloadAnnotationsGRCh37.R

#Filter transcript annotations to GENCODE basic transcripts + exclude retained introns and NMD
Rscript code/geuvadis/filterAnnotationsGRCh37.R

#Extract union exon (and cds) start-end coordinates from txdb 
bsub -G team170 -n1 -R "span[hosts=1] select[mem>12000] rusage[mem=12000]" -q normal -M 12000 -o FarmOut/exonCoords.%J.jobout "/software/R-3.1.2/bin/Rscript code/geuvadis/txdbToExonCoords.R"

#Construct gff files for exons and introns
Rscript code/geuvadis/constructIntronExonGFFs.R

#Count reads overlapping exon annotations
cut -f1 geuvadis_sample_list.txt | python ~/software/utils/submitJobs_EBI.py --MEM 1000 --jobname exonCounts --farmout /nfs/nobackup/stegle/users/kaur/out/ --queue research-rh6 --command "python ~/software/utils/bam/bamToExonCounts.py --sampleDir /nfs/research/stegle/projects/1000GenomesRNASeq/data/mRNA/align_genome/20130402-STAR/ --outDir /nfs/nobackup/stegle/users/kaur/counts/ --gtf annotations/geuvadis/Homo_sapiens.GRCh38.79.exons.gff3 --strand 0 --countsSuffix .exon_counts.txt --type exon --bamSuffix .star.bam --execute True --subdir False"

#Count reads overlapping intron annotations
cut -f1 geuvadis_sample_list.txt | head -n1 | python ~/software/utils/submitJobs_EBI.py --MEM 1000 --jobname intronCounts --farmout /nfs/nobackup/stegle/users/kaur/out/ --queue research-rh6 --command "python ~/software/utils/bam/bamToExonCounts.py --sampleDir /nfs/research/stegle/projects/1000GenomesRNASeq/data/mRNA/align_genome/20130402-STAR/ --outDir /nfs/nobackup/stegle/users/kaur/counts/ --gtf annotations/geuvadis/Homo_sapiens.GRCh38.79.introns.gff3 --strand 0 --countsSuffix .intron_counts.txt --type intron --bamSuffix .star.bam --execute True --subdir False"

#Construct intron events from the intron and exon counts
cut -f1 geuvadis_sample_list.txt| python ~/software/utils/bam/constructIntronEvents.py --intronGFF annotations/geuvadis/Homo_sapiens.GRCh38.79.introns.gff3 --exonGFF annotations/geuvadis/Homo_sapiens.GRCh38.79.exons.gff3--sampleDir /nfs/nobackup/stegle/users/kaur/counts/

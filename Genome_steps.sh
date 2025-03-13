## Below is the steps used to assemble the final draft genome of Gordius panarensis with a brief explanation of each step ##
# the code was ran on scripts on NeSI and as such the module calling is related to NeSI's structure

#raw nanopore reads called using guppy 6.2.1
guppy_basecaller --device auto -i raw_nanopore/W136 -s raw_nanopore/W136_called_guppy621 --flowcell FLO-MIN106 --kit SQK-LSK109 --num_callers 4

#cleaned using porechop
module load Porechop/0.2.4-gimkl-2020a-Python-3.8.2
porechop -i W136_guppy621.fastq -o W136_guppy621.porechop.fastq -t 8

#assembled using flye
module load Flye/2.9.1-gimkl-2022a-Python-3.10.5
flye --nano-hq W136_guppy621.porechop.fastq --out-dir /nesi/nobackup/uoo00108/eddy_hairworm/flye_W136_porechop  -t 8 -i 3

#four rounds of racon
module load Racon
module load minimap2
#four rounds:
#minimap2 flye_W136_porechop.fasta \
#W136_guppy621.porechop.fastq  > flye_W136_porechop.minimap.racon1.paf
racon W136_guppy621.porechop.fastq flye_W136_porechop.minimap.racon1.paf flye_W136_porechop.fasta > flye_W136_porechop.racon1.fasta
##round 2
minimap2 flye_W136_porechop.racon1.fasta \
W136_guppy621.porechop.fastq  > flye_W136_porechop.minimap.racon2.paf
racon W136_guppy621.porechop.fastq flye_W136_porechop.minimap.racon2.paf flye_W136_porechop.racon1.fasta > flye_W136_porechop.racon2.fasta
#round 3
minimap2 flye_W136_porechop.racon2.fasta \
W136_guppy621.porechop.fastq  > flye_W136_porechop.minimap.racon3.paf
racon W136_guppy621.porechop.fastq flye_W136_porechop.minimap.racon3.paf flye_W136_porechop.racon2.fasta > flye_W136_porechop.racon3.fasta
#round 4
minimap2 flye_W136_porechop.racon3.fasta W136_guppy621.porechop.fastq  > flye_W136_porechop.minimap.racon4.paf
racon W136_guppy621.porechop.fastq flye_W136_porechop.minimap.racon4.paf flye_W136_porechop.racon3.fasta > flye_W136_porechop.racon4.fasta

#medaka2
module purge
export PYTHONNOUSERSITE=1
module load Miniconda3/
CONDA_PATH=$(conda info | grep -i 'base environment' | awk '{print $4}')
source $CONDA_PATH/etc/profile.d/conda.sh
source activate /nesi/project/uoo00108/bin/miniconda_envs/medaka2

#correct model is: r941_min_hac_g507
medaka_consensus -i W136_guppy621.porechop.fastq -d flye_W136_porechop.racon4.fasta -o flye_W136_porechop_racon_medaka/medaka_consensus -t 10 -m r941_min_hac_g507

#purgedup
module load zlib
module load minimap2
module load Python

/nesi/project/uoo00108/bin/purge_dups/scripts/run_purge_dups.py config.json /nesi/project/uoo00108/bin/purge_dups/bin/ hairworm_nanoporebase
minimap2 -t 4 -xmap-ont flye_W136_porechop_racon_medaka.fasta W136_guppy621.porechop.fastq | gzip -c - > W136_guppy621.porechop.fastq.paf.gz

#bring in 10x linked reads processed through longranger basic
/nesi/project/uoo00108/bin/longranger-2.2.2/longranger basic --id W136_longranger_basic --fastq=/nesi/nobackup/uoo00108/eddy_hairworm/2023_combined/10x_reads_w136

#scaff10x with longrnger data and input.dat (note longread 1 is filtering mapping score on small contigs not to do with input data type)
/nesi/project/uoo00108/bin/Scaff10X/src/scaff10x -nodes 10 -longread 1 -data input.dat flye_W136_porechop_racon_medaka_purgeddup.fasta flye_W136_porechop_racon_medaka_purgeddup_scaff10x.fasta

#tigmint and arcs
cp W136_longranger_basic.fq.qz myreads.fq.gz
module load Miniconda3/
#CONDA_PATH=$(conda info | grep -i 'base environment' | awk '{print $4}')
#source $CONDA_PATH/etc/profile.d/conda.sh
source activate /nesi/project/uoo00108/bin/miniconda_envs/arcs
cp flye_W136_porechop_racon_medaka_purgeddup_scaff10x.fasta myassembly.fa
tigmint-make tigmint t=10 draft=myassembly reads=myreads
cp myassembly.tigmint.fa my_scaffolds.fa

#then arcs seen as they wouldnt run together:
module load Miniconda3/
module load Python
module load minimap2
CONDA_PATH=$(conda info | grep -i 'base environment' | awk '{print $4}')
source $CONDA_PATH/etc/profile.d/conda.sh
####conda create -p /nesi/project/uoo00108/bin/miniconda_envs/arcs
source activate /nesi/project/uoo00108/bin/miniconda_envs/arcs
arcs-make arks draft=my_scaffolds reads=myreads k=60 t=10

#blobtools
#removed contigs less than 50bp
module load minimap2
module load SAMtools

#maping ont reads back for coverage
minimap2 -ax map-ont -t 10 flye_W136_porechop_racon_medaka_purgeddup_scaff10x_tigmint_arcs_min50.fasta \
                      ../W136_guppy621.porechop.fastq.gz | samtools sort -@10 -O BAM -o flye_W136_porechop_racon_medaka_purgeddup_scaff10x_tigmint_arcs_min50.bam -

module load BLAST
export BLASTDB='/nesi/nobackup/uoo00108/eddy_hairworm/2023_combined/nanopore_pipe/blobtools/Aug_2023/nt/nt'

#blast to nt
blastn -db /nesi/nobackup/uoo00108/eddy_hairworm/2023_combined/nanopore_pipe/blobtools/Aug_2023/nt/nt \
        -task blastn \
        -query flye_W136_porechop_racon_medaka_purgeddup_scaff10x_tigmint_arcs_min50.fasta \
        -outfmt "6 qseqid staxids bitscore std" \
        -max_target_seqs 10 \
        -max_hsps 1 \
        -evalue 1e-25 \
        -num_threads 20 \
        -out results/flye_W136_porechop_racon_medaka_purgeddup_scaff10x_tigmint_arcs_min50_blastn_taskblastn_Aug2023_out

#blastx uniprot reference proteomes
module load DIAMOND
diamond blastx \
              --query flye_W136_porechop_racon_medaka_purgeddup_scaff10x_tigmint_arcs_min50.fasta \
              --db /nesi/nobackup/uoo00108/eddy_hairworm/2023_combined/nanopore_pipe/blobtools/uniprot_dini/reference_proteomes.dmnd \
              --outfmt 6 qseqid staxids bitscore qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore \
              --sensitive \
              --max-target-seqs 1 \
              --evalue 1e-25 \
              --threads 20 \
              --memory-limit 19G \
              --out flye_W136_porechop_racon_medaka_purgeddup_scaff10x_tigmint_arcs_min50.diamond.blastx.out

#busco for blobtools
module load BUSCO
busco -m genome -i flye_W136_porechop_racon_medaka_purgeddup_scaff10x_tigmint_arcs_min50.fasta -o flye_W136_porechop_racon_medaka_purgeddup_scaff10x_tigmint_arcs_min50_busco -l eukaryota_odb10 -f

#create blobtools database
export PATH="/nesi/project/uoo00108/bin/miniconda_envs/blobtools2:$PATH"
export PATH="/nesi/nobackup/uoo00108/eddy_hairworm/2023_combined/nanopore_pipe/blobtools:$PATH"
export PATH="/nesi/project/uoo00108/bin/miniconda_envs/blobtools2/lib/python3.9/site-packages:$PATH"

module load Miniconda3/
module load Python
module load minimap2
CONDA_PATH=$(conda info | grep -i 'base environment' | awk '{print $4}')
source $CONDA_PATH/etc/profile.d/conda.sh
conda activate /nesi/project/uoo00108/bin/miniconda_envs/blobtools2

blobtools create \
       --fasta flye_W136_porechop_racon_medaka_purgeddup_scaff10x_tigmint_arcs_min50.fasta \
       --meta hairworm_assembly.yaml \
       --taxid 190565 \
       --taxdump taxdump \
       flye_W136_porechop_racon_medaka_purgeddup_scaff10x_tigmint_arcs_min50_Assembly

# Second, add other info init.

#if had diamond would add second hit file here
blobtools add \
       --hits results/flye_W136_porechop_racon_medaka_purgeddup_scaff10x_tigmint_arcs_min50_blastn_out \
       --hits flye_W136_porechop_racon_medaka_purgeddup_scaff10x_tigmint_arcs_min50.diamond.blastx.out \
       --busco flye_W136_porechop_racon_medaka_purgeddup_scaff10x_tigmint_arcs_min50_busco/run_eukaryota_odb10/full_table.tsv \
       --cov flye_W136_porechop_racon_medaka_purgeddup_scaff10x_tigmint_arcs_min50.bam  \
       --taxrule bestsumorder \
        flye_W136_porechop_racon_medaka_purgeddup_scaff10x_tigmint_arcs_min50_Assembly

#blast to hairworm genomes
#database all available hairworms at time of download (two american and spanish)
module load BLAST
export BLASTDB='/nesi/nobackup/uoo00108/eddy_hairworm/2023_combined/nanopore_pipe/blobtools/hairworm_genomes/hairworm_genomes'

blastn -db /nesi/nobackup/uoo00108/eddy_hairworm/2023_combined/nanopore_pipe/blobtools/hairworm_genomes/hairworm_genomes \
        -task blastn \
        -query flye_W136_porechop_racon_medaka_purgeddup_scaff10x_tigmint_arcs_min50.fasta \
        -outfmt "6 qseqid staxids bitscore std" \
        -max_target_seqs 10 \
        -max_hsps 1 \
        -evalue 1e-25 \
        -num_threads 20 \
        -out results/flye_W136_porechop_racon_medaka_purgeddup_scaff10x_tigmint_arcs_min50_blastn_hairwormgenomes_taskblastn_Aug2023_out

#grab reads out and put into two fasta files one hits to hairworm one no hits
seqtk subseq fasta reads > output

#start reads hitting hairworm
module load minimap2
module load SAMtools

minimap2 -ax map-ont -t 10 flye_W136_porechop_racon_medaka_purgeddup_scaff10x_tigmint_arcs_min50_hithairworm.fasta \
                      ../../W136_guppy621.porechop.fastq.gz | samtools sort -@10 -O BAM -o flye_W136_porechop_racon_medaka_purgeddup_scaff10x_tigmint_arcs_min50_hithairworm.bam -

#pull blast hits for these scaffolds
grep -f flye_W136_porechop_racon_medaka_purgeddup_scaff10x_tigmint_arcs_min50_blastn_hairwormgenomes_taskblastn_Aug2023_out_scaffolds ../results/flye_W136_porechop_racon_medaka_purgeddup_scaff10x_tigmint_arcs_min50_blastn_out > flye_W136_porechop_racon_medaka_purgeddup_scaff10x_tigmint_arcs_min50_blastn_hithairworm_out
#this is all scaffolds so those ones with no hit had no hit to anything

#pull hits out for diamond
grep -f flye_W136_porechop_racon_medaka_purgeddup_scaff10x_tigmint_arcs_min50_blastn_hairwormgenomes_taskblastn_Aug2023_out_scaffolds ../flye_W136_porechop_racon_medaka_purgeddup_scaff10x_tigmint_arcs_min50.diamond.blastx.out > flye_W136_porechop_racon_medaka_purgeddup_scaff10x_tigmint_arcs_min50.diamond.blastx_hithairworm.out

module load BUSCO
busco -m genome -c 10 -i flye_W136_porechop_racon_medaka_purgeddup_scaff10x_tigmint_arcs_min50_hithairworm.fasta -o flye_W136_porechop_racon_medaka_purgeddup_scaff10x_tigmint_arcs_min50_hithairworm_busco -l eukaryota_odb10 -f

#then for scafolds no hit hairworm
minimap2 -ax map-ont -t 10 flye_W136_porechop_racon_medaka_purgeddup_scaff10x_tigmint_arcs_min50_scaffolds_nohits_hairworm.fasta \
                      ../../W136_guppy621.porechop.fastq.gz | samtools sort -@10 -O BAM -o flye_W136_porechop_racon_medaka_purgeddup_scaff10x_tigmint_arcs_min50_nohit_hairworm.bam -

#pull hits for blastn (there is no hits)
grep -f flye_W136_porechop_racon_medaka_purgeddup_scaff10x_tigmint_arcs_min50_scaffoldsnoblasthairworm ../results/flye_W136_porechop_racon_medaka_purgeddup_scaff10x_tigmint_arcs_min50_blastn_out > flye_W136_porechop_racon_medaka_purgeddup_scaff10x_tigmint_arcs_min50_blastn_nohit_hairworm_out
#pull hits out for diamond
grep -f flye_W136_porechop_racon_medaka_purgeddup_scaff10x_tigmint_arcs_min50_scaffoldsnoblasthairworm ../flye_W136_porechop_racon_medaka_purgeddup_scaff10x_tigmint_arcs_min50.diamond.blastx.out > flye_W136_porechop_racon_medaka_purgeddup_scaff10x_tigmint_arcs_min50.diamond.blastx_nohit_hairworm.out

#check is anything shows up on busco
module load BUSCO
busco -m genome -c 10 -i flye_W136_porechop_racon_medaka_purgeddup_scaff10x_tigmint_arcs_min50_scaffolds_nohits_hairworm.fasta -o flye_W136_porechop_racon_medaka_purgeddup_scaff10x_tigmint_arcs_min50_nohit_hairworm_busco -l eukaryota_odb10 -f

#manually check if RNAseq maps to anything not hitting hairworm in a manner consistent with genes (e.g. not a holopolymer map)
#make star genomes for both no hits and hits and all reads (just showing no hits here as example)
module load STAR
STAR --runThreadN 10 --runMode genomeGenerate --genomeDir flye_W136_porechop_racon_medaka_purgeddup_scaff10x_tigmint_arcs_min50_scaffolds_nohits_hairworm_stargenome --genomeFastaFiles flye_W136_porechop_racon_medaka_purgeddup_scaff10x_tigmint_arcs_min50_scaffolds_nohits_hairworm.fasta

#mapping

STAR --runThreadN 10 --genomeDir flye_W136_porechop_racon_medaka_purgeddup_scaff10x_tigmint_arcs_min50_hithairworm_stargenome --readFilesCommand zcat \
 --outFileNamePrefix flye_W136_porechop_racon_medaka_purgeddup_scaff10x_tigmint_arcs_min50_hithairworm_RNA --outSAMtype BAM SortedByCoordinate \
 --outFilterScoreMinOverLread 0.15 --outFilterMatchNminOverLread 0.15 --outReadsUnmapped Fastx --outFilterMismatchNmax 15 \
 --readFilesIn ../../../Weta_RNAseq_fastptrim.fastq.gz

STAR --runThreadN 10 --genomeDir flye_W136_porechop_racon_medaka_purgeddup_scaff10x_tigmint_arcs_min50_scaffolds_nohits_hairworm_stargenome --readFilesCommand zcat \
 --outFileNamePrefix flye_W136_porechop_racon_medaka_purgeddup_scaff10x_tigmint_arcs_min50_nohithairworm_RNA --outSAMtype BAM SortedByCoordinate \
 --outFilterScoreMinOverLread 0.15 --outFilterMatchNminOverLread 0.15  --outReadsUnmapped Fastx --outFilterMismatchNmax 15 \
 --readFilesIn ../../../Weta_RNAseq_fastptrim.fastq.gz

STAR --runThreadN 10 --genomeDir flye_W136_porechop_racon_medaka_purgeddup_scaff10x_tigmint_arcs_min50_stargenome  --readFilesCommand zcat \
 --outFileNamePrefix flye_W136_porechop_racon_medaka_purgeddup_scaff10x_tigmint_arcs_min50_RNA --outSAMtype BAM SortedByCoordinate \
 --outFilterScoreMinOverLread 0.15 --outFilterMatchNminOverLread 0.15  --outReadsUnmapped Fastx --outFilterMismatchNmax 15 \
 --readFilesIn ../../../Weta_RNAseq_fastptrim.fastq.gz

#unmapped reads from hit hairworm map to no hit hairworm map
STAR --runThreadN 10 --genomeDir flye_W136_porechop_racon_medaka_purgeddup_scaff10x_tigmint_arcs_min50_scaffolds_nohits_hairworm_stargenome \
 --outFileNamePrefix flye_W136_porechop_racon_medaka_purgeddup_scaff10x_tigmint_arcs_min50_nohithairworm_readsunmappedhithairworm --outSAMtype BAM SortedByCoordinate \
 --outFilterScoreMinOverLread 0.15 --outFilterMatchNminOverLread 0.15  --outReadsUnmapped Fastx --outFilterMismatchNmax 15 \
 --readFilesIn flye_W136_porechop_racon_medaka_purgeddup_scaff10x_tigmint_arcs_min50_hithairworm_RNAUnmapped.out.mate1.fastq

#manually checked for potential gene regions and contigs added back in
#must have done this adding back in on the command line (have contig list on dropbox or just compare rescued to none rescued intermediatary fasta files)

#rails
module load SAMtools/1.8-gimkl-2018b
module load Perl/5.30.1-GCC-9.2.0
module load minimap2

export PATH="/nesi/project/uoo00108/bin/RAILS/bin:$PATH"
sh /nesi/project/uoo00108/bin/RAILS/bin/runRAILSminimapSTREAM.sh flye_W136_porechop_racon_medaka_purgeddup_scaff10x_tigmint_arcs_min50_hithairworm_rescues.fasta W136_guppy621.porechop.fa 250 0.90 250 4 ont /opt/nesi/CS400_centos7_bdw/SAMtools/1.16.1-GCC-11.3.0/bin/samtools 10

#ragtag
#patching across the two genome builds (the other one used a merged 10x and nanopore genome)
#ragtag is target query with target being patched
ragtag.py patch flye_W136_porechop_racon_medaka_purgeddup_scaff10x_tigmint_arcs_min50_hithairworm_rescues_rail.fasta flye_W136_porechop_racon_medaka_merged350pseudohap_quickmerge_scaff10x_tigmint_arcs_min50_hithairworm_rescues_rails.fasta

#nextpolish
module load Python
/nesi/project/uoo00108/bin/NextPolish/nextPolish run.cfg


############################################annotation#######################################

#earl grey for repeats
module load Miniconda3
source activate earlgrey
earlGrey -g flye_W136_porechop_racon_medaka_purgeddup_scaff10x_tigmint_arcs_min50_hithairworm_rescues_rail_ragtag_nextpolish_renamed.fasta -s gordiusParanensis -o ./earlGreyOutputs -t 20 -d yes

#use soft masked genome output

#organise hints with prothint
#the hint file metazoa_nematomorpha is comprised of Metazoa.fa (as suggested by prothint) and GCA_954871325.2_tfGorSpeb1-WG-v2_protein.faa (protein sequences from available nematomorpha on genbank)
module load BRAKER
module load ProtHint
prothint.py ../../final_mask/flye_W136_porechop_racon_medaka_purgeddup_scaff10x_tigmint_arcs_min50_hithairworm_rescues_rail_ragtag_nextpolish_renamed.fasta.masked \
metazoa_nematomorpha.fa --threads 20

#running braker 
#braker will do everything so:
module load BRAKER
module load sratoolkit
export AUGUSTUS_CONFIG_PATH=/nesi/nobackup/uoo00108/eddy_hairworm/2023_combined/nanopore_pipe/earlgrey/braker/config
export PATH=$PATH:/nesi/project/uoo00108/bin/GeneMark-ETP/bin
braker.pl --genome=flye_W136_porechop_racon_medaka_purgeddup_scaff10x_tigmint_arcs_min50_hithairworm_rescues_rail_ragtag_nextpolish_renamed.fasta.masked --species=Gordius_paranensis \
    --prot_seq=prothint/metazoa_nematomorpha.fa --rnaseq_sets_ids=Weta_RNAseq \
    --rnaseq_sets_dir=/nesi/nobackup/uoo00108/eddy_hairworm/2023_combined/trinity \
    --threads 20  \
    --GENEMARK_PATH=/nesi/project/uoo00108/bin/GeneMark-ETP/bin
#     --GENEMARK_PATH=/opt/nesi/CS400_centos7_bdw/GeneMark-ES/4.71-GCC-11.3.0/

#summarise genes in gtf
#awk '{a[$3]++}END{for(k in a){print k,a[k]}}' braker.gtf


#always ends up being low so rerun based on this:
#https://github.com/Gaius-Augustus/BRAKER/issues/612

#TSEBRA tends to discard all transcripts that have no evidence. So if the
#evidence situation is comparably poor, you lose too many transcripts. In
#the Supplementary of the GALBA preprint, I show a command line that
#enforces keeping one gene set regardless the evidence (

#I think given that we are using RNA seq from one life stage (no egg, larvae, cyst, very small or emerged adult) keeping them is fine.

tsebra.py -g braker.gtf -k Augustus/augustus.hints.gtf -e hintsfile.gff -c default.cfg -o tsebra_rerun.gtf
rename_gtf.py --gtf tsebra_rerun.gtf --out tsebra_rerun_renamed.gtf
#either convert it here to gff
cat  tsebra_rerun_renamed.gtf  | perl -ne 'if(m/\tAUGUSTUS\t/ or m/\tGeneMark\.hmm\t/ or m/\tGeneMark\.hmm3\t/ or m/\tgmst\t/) {print $_;}' | gtf2gff.pl --gff3 --out=tsebra_rerun_renamed.gff
module load AGAT
agat_convert_sp_gxf2gxf.pl --gtf tsebra_rerun_renamed.gff -o tsebra_rerun_renamedgff_agatgff.gff
#convert to proper gtf
agat_convert_sp_gxf2gxf.pl --gtf tsebra_rerun_renamedgff_agatgff.gff -o tsebra_rerun_renamedgff_agatgff.gtf
#regenerate aa file
#create a wrapped fasta
seqtk seq -l 60 flye_W136_porechop_racon_medaka_purgeddup_scaff10x_tigmint_arcs_min50_hithairworm_rescues_rail_ragtag_nextpolish_renamed.fasta > flye_W136_porechop_racon_medaka_purgeddup_scaff10x_tigmint_arcs_min50_hithairworm_rescues_rail_ragtag_nextpolish_renamed_wrapped.fasta
#generate cds
agat_sp_extract_sequences.pl --gff tsebra_rerun_renamedgff_agatgff.gff -f flye_W136_porechop_racon_medaka_purgeddup_scaff10x_tigmint_arcs_min50_hithairworm_rescues_rail_ragtag_nextpolish_renamed_wrapped.fasta -o tsebra_rerun_renamedgff_agatgff_cds.fa
#generate protein
agat_sp_extract_sequences.pl --gff tsebra_rerun_renamedgff_agatgff.gff -f flye_W136_porechop_racon_medaka_purgeddup_scaff10x_tigmint_arcs_min50_hithairworm_rescues_rail_ragtag_nextpolish_renamed_wrapped.fasta -p -o tsebra_rerun_renamedgff_agatgff_protein.fa

#have then had to remove * from the protein sequences to run interproscan as it threw a hissy about the *

#annotating genes with two methods interproscan and blastp
module load InterProScan/5.64-96.0-gimkl-2022aPerl-5.34.1-Python-3.11.3
interproscan.sh -cpu 20 \
-i tsebra_rerun_renamedgff_agatgff_protein_nostarsinterpro.fa -dp -iprlookup --goterms --pathways

#blastp to swissprot
module load BLAST
module load BLASTDB
blastp -num_threads 20 -db uniprot_sprot.fasta -query tsebra_rerun_renamedgff_agatgff_protein_nostarsinterpro.fa -out tsebra_rerun_renamedgff_agatgff_protein_nostarsinterpro_blastp_swissproteinv2.txt -outfmt 6 -max_target_seqs 100
#have also set a blast to trembl and refseq but probably swissprot is the most conservative for annotation

#merge into each other with AGAT seen as gff outputed by interproscan is failing all checks
module load AGAT
agat_sp_manage_functional_annotation.pl -f tsebra_rerun_renamedgff_agatgff.gff -b tsebra_rerun_renamedgff_agatgff_protein_nostarsinterpro_blastp_swissproteinv2.txt --db uniprot_sprot.fasta -i tsebra_rerun_renamedgff_agatgff_protein_nostarsinterpro.fa.tsv --output tsebra_rerun_renamedgff_agatgff_annotatedinterproscan_swissprotblastp.gff

agat_sp_statistics.pl --gff tsebra_rerun_renamedgff_agatgff.gff

#statistics on gene annotation
agat_sp_functional_statistics.pl --gff tsebra_rerun_renamedgff_agatgff.gff

#for profiling the transcriptome (not used in any step other than to run busco stats on transcriptome as braker takes raw reads)
module load Trinity

Trinity --seqType fq --max_memory 100G --CPU 20 \
--single Weta_RNAseq.fastq.gz \
--min_contig_length 200 \
--trimmomatic \
--output trinity_out


#prepare to upload to genbank
#made folder containing:
#flye_W136_porechop_racon_medaka_purgeddup_scaff10x_tigmint_arcs_min50_hithairworm_rescues_rail_ragtag_nextpolish_renamed.fsa - genome
#tsebra_rerun_renamedgff_agatgff.gff - gff
#template.sbt - ncbi template form. Project number is PRJNA1234713
tbl2asn -indir tbl2asn_ncbi_2/ -t template.sbt -M n -Z -locus-tag-prefix PRJNA1234713 -j '[organism=Gordius paranensis]'

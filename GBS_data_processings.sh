# Commands used to process GBS data
#data was run on NeSI and commands reflect that environment
#two enzymes pstI and mspI
#samples sequenced at AgResearch, they call it GBS but its dd-RAD seq really 

#demultiplex in stacks
module load Stacks
process_radtags -p illumina/ -o dual_pst1_msp1/ -b barcodes_stacks.txt -r --renz_1 pstI --renz_2 mspI --threads 10

#mapping to the assembled genome
module load SAMtools
module load BWA
bwa index flye_W136_porechop_racon_medaka_purgeddup_scaff10x_tigmint_arcs_min50_hithairworm_rescues_rail_ragtag_nextpolish_renamed.fasta
for line in $(<file_list_hairworm.txt);
do
bwa mem -t 8 flye_W136_porechop_racon_medaka_purgeddup_scaff10x_tigmint_arcs_min50_hithairworm_rescues_rail_ragtag_nextpolish_renamed.fasta "${line}.fq.gz" | samtools sort -@8 -o "${line}.bam"
done
#where file_list_hairworm.txt is a list of all the demultiplexed files from proccess_radtags

#generate ngsrelate input files in ANGSD (allele frequency file)
module load angsd
/nesi/project/uoo00108/bin/angsd940/angsd/angsd -b bamlist -gl 2 -domajorminor 1 -snp_pval 1e-6 -domaf 1 -minmaf 0.05 -doGlf 3 -minInd 2 -P 8 -out angsd_allfiles

#run NGSrelate
#extract the frequency column from the allele frequency file and remove the header (to make it in the format NgsRelate needs)
zcat angsd_allfiles.mafs.gz | cut -f5 |sed 1d >freq

### run NgsRelate
./ngsrelate  -g angsdput.glf.gz -n 56 -f freq  -O newres

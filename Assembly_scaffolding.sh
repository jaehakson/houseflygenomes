#!/bin/bash
#SBATCH --job-name=hifiasm
#SBATCH --nodes=1                 # Number of nodes you require
#SBATCH --ntasks=1                # Total # of tasks across all nodes
#SBATCH --cpus-per-task=32         # Cores per task (>1 if multithread tasks)
#SBATCH --mem=128000                # Real memory (RAM) required (MB)
#SBATCH --time=128:00:00           # Total run time limit (HH:MM:SS)
#SBATCH --output=/project/meisel/users/JaeHakSon/HouseflyGenomeHiFiHiC/M5/hifiasm_hic/hifiasm_slurm.%N.%j.out  # STDOUT output file
#SBATCH --error=/project/meisel/users/JaeHakSon/HouseflyGenomeHiFiHiC/M5/hifiasm_hic/hifiasm_slurm.%N.%j.err   # STDERR output file (optional)


#### Initial Assembly with 'hifiasm'
cd /project/meisel/users/JaeHakSon/HouseflyGenomeHiFiHiC/M5/hifiasm_hic/

/project/meisel/users/JaeHakSon/software/hifiasm/hifiasm -o M5_male -t 32 \
--h1 /project/meisel/HouseflyGenomeHiFiHiC/35-4_M5/35-4males-2_1469108/Hi-C/35-4males-2_1469109_S3HiC_R1.fastq.gz \
--h2 /project/meisel/HouseflyGenomeHiFiHiC/35-4_M5/35-4males-2_1469108/Hi-C/35-4males-2_1469109_S3HiC_R2.fastq.gz \
/project/meisel/HouseflyGenomeHiFiHiC/M5/SRR27707824.fastq.gz \
/project/meisel/HouseflyGenomeHiFiHiC/M5/SRR27707825.fastq.gz \
/project/meisel/HouseflyGenomeHiFiHiC/M5/SRR27707826.fastq.gz

awk '/^S/{print ">"$2;print $3}' M5_male.hic.p_ctg.gfa > M5_male.hic.p_ctg.fa
awk '/^S/{print ">"$2;print $3}' M5_male.hic.hap1.p_ctg.gfa > M5_male.hic.hap1.p_ctg.fa
awk '/^S/{print ">"$2;print $3}' M5_male.hic.hap2.p_ctg.gfa > M5_male.hic.hap2.p_ctg.fa


#### Purge duplications
module load minimap2

# Step 1. Run minimap2 to align pacbio data and generate paf files, then calculate read depth histogram and base-level read depth
cd /project/meisel/users/JaeHakSon/HouseflyGenomeHiFiHiC/M5/purge_dups
minimap2 -x map-hifi -t 32 ../hifiasm_hic/M5_male.hic.p_ctg.fa /project/meisel/HouseflyGenomeHiFiHiC/M5/SRR*.fastq.gz | gzip -c - > M5_male_hifi_hic.paf.gz

/project/meisel/users/JaeHakSon/software/purge_dups/bin/pbcstat M5_male_hifi_hic.paf.gz	#(produces PB.base.cov and PB.stat files)
/project/meisel/users/JaeHakSon/software/purge_dups/bin/calcuts PB.stat > cutoffs 2>calcults.log


# Step 2. Split an assembly and do a self-self alignment
/project/meisel/users/JaeHakSon/software/purge_dups/bin/split_fa ../hifiasm_hic/M5_male.hic.p_ctg.fa > M5_male.hic.p_ctg.fa.split
minimap2 -x asm5 -DP M5_male.hic.p_ctg.fa.split M5_male.hic.p_ctg.fa.split | gzip -c - > M5_male.hic.p_ctg.fa.split.self.paf.gz
rm M5_male.hic.p_ctg.fa.split

# Step 3. Purge haplotigs and overlaps
/project/meisel/users/JaeHakSon/software/purge_dups/bin/purge_dups -2 -T cutoffs -c PB.base.cov M5_male.hic.p_ctg.fa.split.self.paf.gz > dups.bed 2> purge_dups.log


# Step 4. Get purged primary and haplotig sequences from draft assembly.
/project/meisel/users/JaeHakSon/software/purge_dups/bin/get_seqs -e dups.bed ../hifiasm_hic/M5_male.hic.p_ctg.fa	# this command create 'purged.fa' and 'hap.fa'


#### RUN YaHS to scaffold contigs, assembled by 'hifiasm' and purged by 'purged_dups'
module load BWA
module load SAMtools
module load jdk/22

cd /project/meisel/users/JaeHakSon/HouseflyGenomeHiFiHiC/M5/yahs/
mkdir raw_bam
mkdir filter_bam
mkdir tmp_dir
mkdir pair_dir
mkdir rep_dir

bwa index -a bwtsw -p bwa_idx_M5_hifi_hic ../purge_dups/purged.fa

# Step 1: FASTQ to BAM
bwa mem -t 32 bwa_idx_M5_hifi_hic /project/meisel/HouseflyGenomeHiFiHiC/35-4_M5/35-4males-2_1469108/Hi-C/35-4males-2_1469109_S3HiC_R1.fastq.gz | samtools view -@ 16 -Sb - > ./raw_bam/M5_ctg_hic_1.bam
bwa mem -t 32 bwa_idx_M5_hifi_hic /project/meisel/HouseflyGenomeHiFiHiC/35-4_M5/35-4males-2_1469108/Hi-C/35-4males-2_1469109_S3HiC_R2.fastq.gz | samtools view -@ 16 -Sb - > ./raw_bam/M5_ctg_hic_2.bam

# Step 2: Filter 5' end
samtools view -h ./raw_bam/M5_ctg_hic_1.bam | perl /project/meisel/users/JaeHakSon/software/yahs/mapping_pipeline/filter_five_end.pl | samtools view -Sb - > ./filter_bam/M5_ctg_hic_1.bam
samtools view -h ./raw_bam/M5_ctg_hic_2.bam | perl /project/meisel/users/JaeHakSon/software/yahs/mapping_pipeline/filter_five_end.pl | samtools view -Sb - > ./filter_bam/M5_ctg_hic_2.bam

# Step 3A: Pair reads & mapping quality filter
samtools faidx ../purge_dups/purged.fa
perl /project/meisel/users/JaeHakSon/software/yahs/mapping_pipeline/two_read_bam_combiner.pl ./filter_bam/M5_ctg_hic_1.bam ./filter_bam/M5_ctg_hic_2.bam samtools 10 | samtools view -bS -t ../purge_dups/purged.fa.fai - | samtools sort -@ 16 -o ./tmp_dir/M5_ctg_hic.bam -

# Step 3.B: Add read group
java -Xmx64G -Djava.io.tmpdir=temp/ -jar /project/meisel/users/JaeHakSon/software/yahs/mapping_pipeline/picard.jar AddOrReplaceReadGroups INPUT=./tmp_dir/M5_ctg_hic.bam OUTPUT=./pair_dir/M5_ctg_hic.bam ID=M5_ctg_hic LB=M5_ctg_hic SM=M5_scaff PL=ILLUMINA PU=none

# Step 4: Mark duplicates
java -Xmx96G -XX:-UseGCOverheadLimit -Djava.io.tmpdir=temp/ -jar /project/meisel/users/JaeHakSon/software/yahs/mapping_pipeline/picard.jar MarkDuplicates INPUT=./pair_dir/M5_ctg_hic.bam OUTPUT=./rep_dir/M5_ctg_hic_dedup.bam METRICS_FILE=./rep_dir/metrics.M5_ctg_hic_dedup.txt TMP_DIR=./tmp_dir ASSUME_SORTED=TRUE VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=TRUE
samtools index ./rep_dir/M5_ctg_hic_dedup.bam
perl /project/meisel/users/JaeHakSon/software/yahs/mapping_pipeline/get_stats.pl ./rep_dir/M5_ctg_hic_dedup.bam > ./rep_dir/M5_ctg_hic_dedup.bam.stats

rm -r *_bam
rm -r tmp_dir
rm -r pair_dir

# Step 5: YaHS
mkdir yahs_scaff
/project/meisel/users/JaeHakSon/software/yahs/yahs ../purge_dups/purged.fa ./rep_dir/M5_ctg_hic_dedup.bam -o ./yahs_scaff/yahs_out_M5_hifi_hic

# Step 6: For manual curation with Juicebox (JBAT)
cd /project/meisel/users/JaeHakSon/HouseflyGenomeHiFiHiC/M5/yahs/
mkdir yahs_hic_manual
/project/meisel/users/JaeHakSon/software/yahs/juicer pre -a -o ./yahs_hic_manual/out_JBAT ./yahs_scaff/yahs_out_M5_hifi_hic.bin ./yahs_scaff/yahs_out_M5_hifi_hic_scaffolds_final.agp ../purge_dups/purged.fa.fai >./yahs_hic_manual/out_JBAT.log 2>&1
java -jar -Xmx64G /project/meisel/users/JaeHakSon/software/yahs/juicer_tools.2.20.00.jar pre ./yahs_hic_manual/out_JBAT.txt ./yahs_hic_manual/out_JBAT.hic.part <(cat ./yahs_hic_manual/out_JBAT.log  | grep PRE_C_SIZE | awk '{print $2" "$3}')
mv ./yahs_hic_manual/out_JBAT.hic.part ./yahs_hic_manual/out_JBAT.hic

# Step 7: Post review
/project/meisel/users/JaeHakSon/software/yahs/juicer post -o ./yahs_hic_manual/out_JBAT ./yahs_hic_manual/out_JBAT.review.assembly ./yahs_hic_manual/out_JBAT.liftover.agp ../purge_dups/purged.fa
mkdir ./yahs_hic_manual/post_review
mv ./yahs_hic_manual/out_JBAT.FINAL* ./yahs_hic_manual/post_review

#!/bin/bash
#SBATCH --job-name=5M_p_matix
#SBATCH --nodes=1                 # Number of nodes you require
#SBATCH --ntasks=1                # Total # of tasks across all nodes
#SBATCH --cpus-per-task=36         # Cores per task (>1 if multithread tasks)
#SBATCH --mem=128000                # Real memory (RAM) required (MB)
#SBATCH --time=24:00:00           # Total run time limit (HH:MM:SS)

cd /project/meisel/users/JaeHakSon/HouseflyGenomeHiFiHiC/M5/HiC/primary/fastq/

source ~/.bashrc
conda activate WaveTAD	# Github repository: https://github.com/ryanpellow84/WaveTAD

awk '{print $1}' /project/meisel/users/JaeHakSon/HouseflyGenomeHiFiHiC/M5/yahs/yahs_hic_manual/post_review3/out_JBAT.FINAL.fa.fai | head -8 > 5M_chrom_only.txt
readarray CHR_ARR < 5M_chrom_only.txt
rm 5M_chrom_only.txt

#Index Reference:
bwa index /project/meisel/users/JaeHakSon/HouseflyGenomeHiFiHiC/M5/yahs/yahs_hic_manual/post_review3/out_JBAT.FINAL.fa

#Map Reads:
bwa mem -t 36 -E 50 -L 0 /project/meisel/users/JaeHakSon/HouseflyGenomeHiFiHiC/M5/yahs/yahs_hic_manual/post_review3/out_JBAT.FINAL.fa /project/meisel/HouseflyGenomeHiFiHiC/35-4_M5/35-4males-2_1469108/Hi-C/35-4males-2_1469109_S3HiC_R1.fastq.gz | samtools view --threads 36 -bS - -o houseflyM5_1469109_S3HiC_1.bam
bwa mem -t 36 -E 50 -L 0 /project/meisel/users/JaeHakSon/HouseflyGenomeHiFiHiC/M5/yahs/yahs_hic_manual/post_review3/out_JBAT.FINAL.fa /project/meisel/HouseflyGenomeHiFiHiC/35-4_M5/35-4males-2_1469108/Hi-C/35-4males-2_1469109_S3HiC_R2.fastq.gz | samtools view --threads 36 -bS - -o houseflyM5_1469109_S3HiC_2.bam

#Build Matrices:
cd /project/meisel/users/JaeHakSon/HouseflyGenomeHiFiHiC/M5/HiC/primary/fastq/
hicFindRestSite -f /project/meisel/users/JaeHakSon/HouseflyGenomeHiFiHiC/M5/yahs/yahs_hic_manual/post_review3/out_JBAT.FINAL.fa -p GATC -o 5M_dpnII_sites.bed
hicBuildMatrix -s houseflyM5_1469109_S3HiC_1.bam houseflyM5_1469109_S3HiC_2.bam -bs 50000 -o 5M_50kb.cool --skipDuplicationCheck --QCfolder ../qcdir/qcfolder_50kb_mat_add --threads 24 -rs 5M_dpnII_sites.bed -seq GATC --danglingSequence GATC

rm *.bam
rm *_dpnII_sites.bed

#Correct Matrices:
hicCorrectMatrix correct --matrix 5M_50kb.cool --chromosomes ${CHR_ARR[@]} -o 5M_50kb_corrected.cool

#Normalize Matrices:
hicNormalize -m 5M_50kb_corrected.cool --normalize norm_range -o 5M_50kb_norm_corrected.cool

rm *kb_corrected.cool
rm *kb.cool

conda deactivate

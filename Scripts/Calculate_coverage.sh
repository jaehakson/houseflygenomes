#!/bin/bash
#SBATCH --job-name=coverage
#SBATCH --nodes=1                 # Number of nodes you require
#SBATCH --ntasks=1                # Total # of tasks across all nodes
#SBATCH --cpus-per-task=32         # Cores per task (>1 if multithread tasks)
#SBATCH --mem=256000                # Real memory (RAM) required (MB)
#SBATCH --time=5:00:00           # Total run time limit (HH:MM:SS)
#SBATCH --output=/project/meisel/users/JaeHakSon/HouseflyGenomeHiFiHiC/M5/yahs/yahs_hic_manual/post_review3/coverage.%N.%j.out  # STDOUT output file
#SBATCH --error=/project/meisel/users/JaeHakSon/HouseflyGenomeHiFiHiC/M5/yahs/yahs_hic_manual/post_review3/coverage.%N.%j.err   # STDERR output file (optional)

module add SAMtools
module add BEDTools
module add minimap2

##### 5M
cd /project/meisel/users/JaeHakSon/HouseflyGenomeHiFiHiC/M5/yahs/yahs_hic_manual/post_review3/
minimap2 -t 32 -ax map-hifi out_JBAT.FINAL.fa /project/meisel/HouseflyGenomeHiFiHiC/M5/*.fastq.gz | sed '7699475d' - | samtools view -h -q 20 - | samtools view -b - > 5M_hifi.sort.bam
samtools index 5M_hifi.sort.bam
bedtools coverage -a <(bedtools makewindows -g out_JBAT.FINAL.fa.fai -w 10000 | sort -k1,1 -k2,2n) -b 5M_hifi.sort.bam -mean | awk '{OFS="\t"; print $1,($2+$3)/2,$4}' > M5_hifi_all_mean.tsv

rm 5M_hifi.sort.bam 5M_hifi.sort.bam.bai

grep -wE "scaffold_1|scaffold_2|scaffold_3|scaffold_4|scaffold_5" M5_hifi_all_mean.tsv | cut -f3 | sort -n | awk '{
    a[i++] = $1;
} END {
    if (i % 2 == 0) {
        # Even number of elements: median is the average of the two middle elements
        print (a[i/2 - 1] + a[i/2]) / 2;
    } else {
        # Odd number of elements: median is the middle element
        print a[int(i/2)];
    }
}' # 60.9332008
awk '{OFS="\t"; print $1,$2,$3/60.9332008}' M5_hifi_all_mean.tsv > M5_hifi_all_mean_autosome_norm.tsv

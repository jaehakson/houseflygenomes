#!/bin/bash
#SBATCH --job-name=variant_M1
#SBATCH --nodes=1                 # Number of nodes you require
#SBATCH --ntasks=32                # Total # of tasks across all nodes
#SBATCH --mem=256000                # Real memory (RAM) required (MB)
#SBATCH --time=96:00:00           # Total run time limit (HH:MM:SS)
#SBATCH --output=/project/meisel/users/JaeHakSon/HouseflyGenomeHiFiHiC/M1/VariantCall/deepvariant_M1.%N.%j.out  # STDOUT output file
#SBATCH --error=/project/meisel/users/JaeHakSon/HouseflyGenomeHiFiHiC/M1/VariantCall/deepvariant_M1.%N.%j.err   # STDERR output file (optional)

module load minimap2
module load SAMtools

cd /project/meisel/users/JaeHakSon/HouseflyGenomeHiFiHiC/M1/genomes/
minimap2 -t 28 -ax map-hifi ./MdomM1-pctg_v3.fa /project/meisel/HouseflyGenomeHiFiHiC/JaeHiFi_aabys/VPGRU_rawHiFi/HiFi_08262025/m84156_250117_204443_s1.hifi_reads.fastq.gz | \
samtools view -@ 16 -b - | samtools sort -@ 16 - -o 1M_hifi.sort.bam
samtools index 1M_hifi.sort.bam

cd /project/meisel/users/JaeHakSon/HouseflyGenomeHiFiHiC/M1/VariantCall/
export datadir=/project/meisel/
apptainer exec --bind $datadir,$TMPDIR /project/meisel/users/JaeHakSon/software/deepvariant_1.9.0.sif /opt/deepvariant/bin/run_deepvariant \
    --model_type=PACBIO --reads /project/meisel/users/JaeHakSon/HouseflyGenomeHiFiHiC/M1/genomes/1M_hifi.sort.bam \
    --ref=/project/meisel/users/JaeHakSon/HouseflyGenomeHiFiHiC/M1/genomes/MdomM1-pctg_v3.fa \
    --output_vcf=/project/meisel/users/JaeHakSon/HouseflyGenomeHiFiHiC/M1/VariantCall/M1_p_v3.vcf.gz \
    --intermediate_results_dir=/project/meisel/users/JaeHakSon/HouseflyGenomeHiFiHiC/M1/VariantCall/intermidiates

rm /project/meisel/users/JaeHakSon/HouseflyGenomeHiFiHiC/M1/genomes/1M_hifi.sort.bam
rm /project/meisel/users/JaeHakSon/HouseflyGenomeHiFiHiC/M1/genomes/1M_hifi.sort.bam.bai

## USE genotype quality 20
module load VCFtools
zcat M1_p_v3.vcf.gz | grep -E "#|0/1" | grep -E "PASS" | vcftools --vcf - --min-alleles 2 --max-alleles 2 --recode --stdout | \
awk '/^#/ || $6>=20' | vcftools --vcf - --SNPdensity 1000000 --out M1_het_density_GQ20_1Mb

#!/bin/bash
#SBATCH --job-name=M1_kaks_cal
#SBATCH --nodes=4                 # Number of nodes you require
#SBATCH --ntasks=24               # Total # of tasks across all nodes
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=1         # Cores per task (>1 if multithread tasks)
#SBATCH --mem-per-cpu=32000                # Real memory (RAM) required (MB)
#SBATCH --time=144:00:00          # Total run time limit (HH:MM:SS)
#SBATCH --output=/project/meisel/users/JaeHakSon/HouseflyGenomeHiFiHiC/M1/compare_haps/kaks_cal_parallel.%N.%j.out  # STDOUT output file
#SBATCH --error=/project/meisel/users/JaeHakSon/HouseflyGenomeHiFiHiC/M1/compare_haps/kaks_cal_parallel.%N.%j.err   # STDERR output file (optional)


#### Extract representative coding sequences, if multiple isoforms exist in a singel gene, using 'gffread' and 'NCBI annotation table'


#### Split multiple sequences of a fasta file up into individual sequences
cd /project/meisel/users/JaeHakSon/HouseflyGenomeHiFiHiC/M1/compare_haps
mkdir ./MACSE/polished_cds/hap1
mkdir ./MACSE/polished_cds/hap2
/project/meisel/users/JaeHakSon/software/exonerate-2.2.0-x86_64/bin/fastaexplode -f ./liftoff_hap1/M1_hifi_hic_hap1_polished_CDS_longest.fasta -d ./MACSE/polished_cds/hap1
/project/meisel/users/JaeHakSon/software/exonerate-2.2.0-x86_64/bin/fastaexplode -f ./liftoff_hap2/M1_hifi_hic_hap2_polished_CDS_longest.fasta -d ./MACSE/polished_cds/hap2


#### Replace header with hap
cd /project/meisel/users/JaeHakSon/HouseflyGenomeHiFiHiC/M1/compare_haps/MACSE/polished_cds/

for id in $(ls hap1/ | sort | uniq | grep ".fa" )
do
sed -i 's/>/>hap1_/g' hap1/$id
done

for id in $(ls hap2/ | sort | uniq | grep ".fa" )
do
sed -i 's/>/>hap2_/g' hap2/$id
done


#### Aggregate/Concatenate each CDS
mkdir aggregate
for cds in $(ls hap*/ | sort | uniq -d | grep ".fa")
do
cat hap*/$cds > ./aggregate/$cds
done


#### CDS alignment using MACSE
module load jdk 
# or module load Java

mkdir align_NT
mkdir align_AA
for fa in $(ls aggregate/ | grep ".fa" | rev | cut -c4- | rev)
do
java -Xmx112g -jar /project/meisel/users/JaeHakSon/software/macse_v2.07.jar -prog alignSequences -seq aggregate/${fa}.fa -out_NT align_NT/${fa}_NT.fasta -out_AA align_AA/${fa}_AA.fasta
done


#### Post filtering by using MACSE "reportMaskAA2NT" program to report this AA masking/filtering at the nucleotide level
mkdir out_filter_NT
for fa in $(ls aggregate/ | grep ".fa" | rev | cut -c4- | rev)
do
java -jar /project/meisel/users/JaeHakSon/software/macse_v2.07.jar -prog reportMaskAA2NT -align_AA align_AA/${fa}_AA.fasta -align align_NT/${fa}_NT.fasta -mask_AA ! -out_NT out_filter_NT/${fa}_NT.fasta
done


#### parallel run conversion to AXT & KaKs calculator
cd /project/meisel/users/JaeHakSon/HouseflyGenomeHiFiHiC/M1/compare_haps/MACSE/
mkdir ./polished_KaKs_single/
cd ./polished_KaKs_single/

parallel -j 4 perl /project/meisel/users/JaeHakSon/software/KaKs_Calculator3.0/parseFastaIntoAXT.pl {} {/} ::: /project/meisel/users/JaeHakSon/HouseflyGenomeHiFiHiC/M1/compare_haps/MACSE/polished_cds/out_filter_NT/*_NT.fasta

parallel -j 24 /project/meisel/users/JaeHakSon/software/KaKs_Calculator3.0/src/KaKs -i {} -o {}.kaks ::: *.fasta.axt

cat *.fasta.axt.kaks | head -1 > M1_all.kaks && tail -n+2 -q *.fasta.axt.kaks >> M1_all.kaks

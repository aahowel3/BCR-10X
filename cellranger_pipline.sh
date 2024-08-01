#!/bin/sh
   #SBATCH --nodes=1
   #SBATCH --ntasks=1
   #SBATCH --cpus-per-task=4
   #SBATCH --mem=7000
   #SBATCH --time=240:00:00

cd /gpfs/home/ahowell/ca_av_0124_zhu_20240626_0_10x5prime-vdj-only_mouse_AH/analysis2

module load cellranger/8.0.0

cellranger vdj --id=MouseB_Cell \
    --reference=/gpfs/home/ahowell/refdata-cellranger-vdj-GRCm38-alts-ensembl-7.0.0 \
    --fastqs=../seqdata \
    --sample=01id01JZ \
    --localcores=8 \
    --localmem=64 \

# 10X Based B-cell Immune Profiling in UFC RSV-F Trimer Immunized Mouse Splenic Cells

Respiratory syncytial virus (RSV) is a significant contributor to cases of acute lower respiratory infection (ALRI) and primarily affects young children and the elderly (Nam and Ison, [2019](https://pubmed.ncbi.nlm.nih.gov/31506273/); Branche and Falsey, [2015](https://pubmed.ncbi.nlm.nih.gov/25851217/); Shi et al., [2017](https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(22)00478-0/fulltext); Savic et al., [2023](https://pubmed.ncbi.nlm.nih.gov/36369772/)). Worldwide, an estimated 3.2 million hospitalizations and 118,2000 deaths were attributed to RSV in children under 5 years of age in 2015 (Shi et al., [2017](https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(22)00478-0/fulltext)), and in 2019 there were approximately 470,000 RSV-related hospitalizations and 33,000 in-hospital deaths globally due to RSV in adults older than 60 (Savic et al., [2023](https://pubmed.ncbi.nlm.nih.gov/36369772/)). Therefore, the development of an RSV vaccine is of critical importance to global health. Here, we present a 10X Genomics based protocol for isolating and single-cell sequencing monoclonal antibodies (mAbs) after immunization with an uncleaved prefusion-closed (UFC) trimer of the RSV F protein (Lee et al., [2024](https://pubmed.ncbi.nlm.nih.gov/38496645/)) in order the evaluate the vaccine’s ability to induce robust antibody responses with high neutralizing titers. 

## Sample Preperation and Sequencing 

Mouse splenocyte samples (2X 2mL cryopreserved cells approx 30 million cells) were thawed into PBS + 1% BSA + 2mM EDTA buffer and stained with Fixable Aqua dead cell stain (Invitrogen), 2.4G2 Fc receptor block, biotinylated trimer probe, and premium-grade allophycocyanin (APC)-labeled streptavidin (Thermo Fisher) for sorting on a MoFloAstrios EQ (Beckman Coulter) using double positive gating. Approx. 6000 antigen-specific mouse B cells were captured from sorting and loaded onto a 10X Chromium Controller for single-cell GEM bead emulsification. 

## Data Processing with Cellranger

`cellranger vdj` pipeline assembles full V(D)J contigs from fastq reads and determines viable single-cells from GEMs using the supporting number of UMIs and barcodes. <br />

Output files are contained in the folder `MouseB_cell`

## Analysis and Plotting

Using `samtools` MD tags were added to the file `concat_ref.bam` which contained the alignments of the full assembled light and heavy V(D)J contigs relative to the concatenated reference V(D)J contigs. The MD tags contain details not only on the number of matches/mistmatches/indels relative to the reference but also the exact substitutions. <br />

`concat_ref.bam` is one of the standard outputs of the `cellranger vdj` pipeline <br />

`samtools calmd --threads 15 -rb concat_ref.bam concat_ref.fasta > concat_ref_MDtag.bam` <br />

In order to determine the exact mutations and calculate the germline divergence of the V gene per complete cell barcode, the corresponding region of the concatenated reference genome (combined VDJ segements) that an assembled contig aligned to was extracted using the tool `sam2pairwise`. <br />

`samtools view concat_ref_MDtag.bam | sam2pairwise > concat_ref_pairwise_MDtag.out` <br />

Use script `cellranger_mouseRSV_spleen_dataanalysis_plots.R` to select candidate antibodies for gene synthesis and functional analysis based on their clonotype abundance, CDR3 length, and somatic hypermutation of the V gene segments.

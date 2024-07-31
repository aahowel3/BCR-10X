# Cellranger 10X Genomics data processing from RSV trimer immunized mouse BCR

### 10X Based B-cell Immune Profiling in UFC RSV-F Trimer Immunized Mouse Splenic Cells 
https://github.com/aziele/phirbo <br />

Respiratory syncytial virus (RSV) is a significant contributor to cases of acute lower respiratory infection (ALRI) and primarily affects young children and the elderly (Nam and Ison, [2019](https://pubmed.ncbi.nlm.nih.gov/31506273/)); Branche and Falsey, 2015; Shi et al., 2017; Savic et al., 2023). Worldwide, an estimated 3.2 million hospitalizations and 118,2000 deaths were attributed to RSV in children under 5 years of age in 2015 (Shi et al,, 2017), and in 2019 there were approximately 470,000 RSV-related hospitalizations and 33,000 in-hospital deaths globally due to RSV in adults older than 60 (Savic et al., 2023). Therefore, the development of an RSV vaccine is of critical importance to global health. Here, we present a 10X Genomics based protocol for isolating and single-cell sequencing monoclonal antibodies (mAbs) after immunization with an uncleaved prefusion-closed (UFC) trimer of the RSV F protein (Lee et al., 2024) in order the evaluate the vaccine’s ability to induce robust antibody responses with high neutralizing titers. 

Antigen-specific mouse B cells were isolated via 
Download taxdb.btd, taxdb.bti, taxdb.tar.gz from https://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz <br />
Prepare rank-biased overlap blast results via `phirbo_preprocessing.txt` <br />
Run phribo by providing two input directories (i.e., for phages [`phage_virusblast/`] and bacteria [`phage_hostblast/`]) containing ranked lists from blast output, and an output file name (`phage_phirbo/predictions.csv`) <br />

`python phirbo/phirbo.py phage_virusblast/ phage_hostsblast/ phage_phirbo/predictions.csv`
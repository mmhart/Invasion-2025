METADATA AND DATA DICTIONARY — INVASION STUDY, 2025 AMF COMMUNITY SAMPLES

Associated manuscript: Long term persistence of a commercial arbuscular mycorrhizal biofertilizer in a natural grassland
Authors: Joyalea Carson-Austin; Dane Goh; Miranda M. Hart
Affiliation: Biology, University of British Columbia Okanagan, Kelowna, BC, Canada
Corresponding author: miranda.hart@ubc.ca
Journal: Ecosphere
Data file: InvasMeta_Expanded_18S.tsv (tab-separated, UTF-8; one row per listed soil sample, 12 rows)
Original file: InvasMeta.tsv. Original five fields and all 12 data values were retained without alteration. Newly appended fields are study annotations drawn from the manuscript and Appendix S1, Table S11. Seven further columns specify the 18S rRNA metabarcoding performed on all 12 samples.

DATA SCOPE
This file lists twelve soil samples collected in 2025 from Quail (UBC Okanagan West Campus Lands) and KM5 (Knox Mountain Provincial Park), selected for AMF community analysis after determining presence/absence and abundance of the introduced inoculant by ddPCR. Six samples are listed per plot: three with the highest inoculant abundance and three with no detectable inoculant, as described in the manuscript. These are selected samples, not all soil samples collected in the longitudinal field experiment.

STUDY DESIGN / CONTEXT
Two sites in Kelowna, British Columbia, Canada, in the Ponderosa Pine biogeoclimatic zone, featuring open-canopy forest and bluebunch wheatgrass (Pseudoroegneria spicata).
Site-level coordinates (NOT precise sample coordinates): Knox Mountain 49.91331 degrees N, 119.48065 degrees W; UBC Okanagan West Campus Lands 49.93982 degrees N, 119.396764 degrees W.
Original field introduction: July 2012, live or twice-autoclaved commercial Rhizophagus irregularis DAOM197198 inoculum (MYKE PRO GREENHOUSE G) applied to transplanted Artemisia ludoviciana and Achillea millefolium seedlings. The 2025 dataset is 13 years after introduction. Earlier (2014–2017) samples were plant-centred distance cores; 2022 and 2025 used plot transects after some original plants/markers became unavailable. Thus, 2022/2025 data were excluded from the longitudinal statistical models but included in descriptive plots.
2025 sampling, as reported: a 10-m transect across each plot, with soil cores at 1-m intervals, described as ten samples per plot; samples analysed here were selected for 2025 AMF community comparison.
Laboratory: DNA extracted using Qiagen DNeasy PowerSoil kit per manufacturer's instructions. Inoculant quantified with Bio-Rad droplet digital PCR (ddPCR) and genotype-targeted cox3-rnl mitochondrial intergenic region primer/probe (Kokkoris et al. 2019); manuscript describes control reactions, thermal cycling and Bio-Rad QX100 analysis. The manuscript also contains an inconsistent cox1-rnl spelling; see verification notes below.
AMF communities: semi-nested PCR targeting SSU 18S rRNA (AML1/AML2 followed by WANDA/AML2), Illumina MiSeq 2 × 300 bp, QIIME 2 v2024.5, cutadapt and DADA2, taxonomic assignment using SILVA and MaarjAM. Raw reads accession: NCBI BioProject PRJNA1364738 (as stated in manuscript).
Analysis scripts (as stated in manuscript): https://github.com/mmhart/Invasion-2025. This GitHub link is not a substitute for a versioned permanent archive under ESA policy.

COLUMN DICTIONARY
SampleID — String; original sample identifier ending in INV. Identity mappings to Appendix S1 Table S11 were checked using plot and numeric position (e.g., Q/14 corresponds to Q14). Do not interpret the numerical portion of SampleID as a transect position.
position_transect — Numeric, original values 2–28; original field retained unchanged. The manuscript states 1-m intervals along a 10-m transect, but does not establish the coordinate system or units for values up to 28. UNITS AND DEFINITION NEED AUTHOR VERIFICATION; do not assume metres or sample indices.
Inoculant — Categorical y/n; in the provided original table 'y' corresponds to detected inoculant and 'n' to undetected. This is a detection classification for the sample, NOT the 2012 inoculation-treatment assignment.
site — Categorical Q/KM5; original plot codes: Q=Quail at UBC Okanagan; KM5=KM5 at Knox Mountain. This variable identifies a plot, not just a broad site.
CopiesPerGSoil — Numeric; ddPCR-derived estimate of target mitochondrial DNA copies per gram of soil; 0 indicates no detection in the table and should not be construed as a biological guarantee of absence. Original values were preserved.
CollectionYear — Integer; 2025 for all twelve rows (manuscript/Appendix Table S11).
StudyLocation — Text; location/site corresponding to the plot code (manuscript).
PlotName — Text; Quail or KM5 (Appendix Table S11).
InoculantDetected — Binary 1/0; 1 corresponds to 'y' and 0 to 'n' in original field (Appendix Table S11).
AMFCommunitySequencing — Yes/No; Yes for these twelve soil samples (Appendix Table S11 'Community analysis' column). This means AMF community sequencing by 18S rRNA metabarcoding, not genotype-specific ddPCR.
TargetOrganism — Text; Rhizophagus irregularis isolate DAOM197198 (manuscript).
Assay — Text; ddPCR for inoculant detection and quantification (manuscript and Appendix Table S11).
SampleMaterial — Text; soil (manuscript and Appendix Table S11).


Metabarcoding18S — Yes/No; Yes for all twelve samples. Indicates 18S rRNA gene amplicon metabarcoding was performed for resident AMF community characterization.
MetabarcodingTarget — Text; Glomeromycota-specific small-subunit (SSU) 18S rRNA gene region. The sequenced nested-amplicon region is WANDA–AML2.
MetabarcodingPCR1Primers — Text; AML1/AML2, first PCR round.
MetabarcodingPCR2Primers — Text; WANDA/AML2, second PCR round.
MetabarcodingPlatform — Text; Illumina MiSeq.
MetabarcodingReadFormat — Text; 2 x 300 bp paired-end sequencing (MiSeq Reagent Kit v3).
MetabarcodingBioProject — Text; PRJNA1364738, manuscript-stated NCBI BioProject for sequence reads. This is project-level, NOT a per-sample SRA run accession; the sample/run crosswalk remains to be supplied or verified.

18S ILLUMINA METABARCODING: LIBRARY GENERATION AND BIOINFORMATICS
Purpose: Characterize resident arbuscular mycorrhizal fungal (AMF) communities in the 2025 subset, comparing soil with abundant DAOM197198 detection against soil without detected DAOM197198. These community measurements are distinct from ddPCR estimates of inoculant abundance.
DNA: Extracted from soil with Qiagen DNeasy PowerSoil according to manufacturer instructions.
Target and primers: Semi-nested PCR targeting a Glomeromycota-associated portion of SSU (18S) rRNA. First round AML1 (5'-ATCAACTTTCGATGGTAGGATAGA-3') / AML2 (5'-GAACCCAAACACTTTGGTTTCC-3'); second round WANDA (5'-CAGCCGCGGTAATTCCAGCT-3') / AML2. Sequences are transcribed from the manuscript with spaces/hyphen formatting standardized only.
First-round PCR: 50 uL reactions including 2 uL DNA; 30 cycles of 94 C 1 min, 50 C 1 min, 72 C 1 min after 94 C 3 min initial denaturation; 72 C 10 min final extension. Positive and negative controls were included.
Second-round PCR: 50 uL reactions; 30 cycles of 95 C 30 s, 54 C 40 s, 72 C 1 min after 95 C 3 min initial denaturation; 72 C 10 min final extension. Positive and negative controls included. Products cleaned with QIAquick PCR Purification Kit and sent to Dalhousie University's IMR Sequencing Centre for barcoding and sequencing.
Sequencing: Illumina MiSeq, paired-end 2 x 300 bp (MiSeq Reagent Kit v3).
Processing: QIIME 2 v2024.5; adapters/primers removed with cutadapt; DADA2 denoising with forward/reverse truncation at 289/279 bp, maximum expected errors 2/3; chimera filtering during denoising; MAFFT alignment; taxonomy with QIIME 2 feature-classifier classify-sklearn, SILVA 99 v138 and MaarjAM reference databases customized to WANDA-AML2; classification confidence 0.2. ASVs with <=2 reads, non-Glomeromycotan sequences, and singletons removed; 256 ASVs retained, as reported in manuscript.
Diversity: Data rarefied to 3,344 reads/sample for diversity analyses. Alpha diversity: observed richness, Shannon, Pielou's evenness, abundance-based coverage estimator and Faith's phylogenetic diversity. Beta diversity: Bray-Curtis, Jaccard and weighted/unweighted UniFrac. These are methods and dataset-level processing specifications, not sequence/count columns in this sample annotation table.
Sequence archive: NCBI BioProject PRJNA1364738 (stated in manuscript). 

# Invasion 2025

## Study Overview

This study tracks the persistence of an arbuscular mycorrhizal fungal (AMF) biofertilizer over 13 years. It examines the fate of an introduced AMF, *Rhizophagus irregularis* DAOM197198, in natural systems and addresses the following research questions:

- Does the inoculant establish?
- Does it persist over multiple growing seasons?
- Does host plant identity affect establishment or persistence?
- Do uninoculated plants become colonized? If so, when?
- Does inoculation affect resident AMF communities after 13 years?

We planted four plants at each of five sites within two blocks near Kelowna, British Columbia, Canada, in 2013. Plants were inoculated either in the greenhouse or in the field.

Between 2014 and 2017, soil cores were collected at distances of 10, 25, and 50 cm from individual plants. DNA was extracted and pooled by plant, and inoculant abundance was measured using DAOM197198-specific primers and digital droplet PCR (ddPCR).

In 2022 and 2025, remaining plots were resampled to assess inoculant persistence. In 2025, changes in resident AMF communities were also investigated using Illumina MiSeq metabarcoding targeting the SSU 18S rRNA gene (WANDA/AML2; Dumbrell et al. 2011; Lee et al. 2008).

## Sample Information

- **Land type:** Grassland
- **Geographic region:** Okanagan, southern British Columbia, Canada
- **Soil types (pooled):** Loam and sandy loam
- **Sampling locations:**
  - UBC Okanagan: 49.9394° N, 119.3948° W
  - Knox Mountain: 49.911° N, 119.482° W
- **Collection years:** 2014, 2015, 2016, 2017, 2022, and 2025

## DNA Extraction and PCR

- **DNA extraction:** QIAGEN DNeasy PowerSoil Kit (0.2 g soil samples)
- **PCR cleanup:** QIAquick PCR Purification Kit (Qiagen Inc.)

### Digital Droplet PCR (ddPCR)

#### Primer Region

The ddPCR assay targeted the single-copy *cox3–rnl* intergenic region of mitochondrial DNA in *Rhizophagus irregularis* DAOM197198.

The isolate-specific primer/probe assay was developed to detect DAOM197198 while excluding other *R. irregularis* genotypes (Kokkoris et al. 2019).

#### PCR Reaction and Cycling Conditions

Digital droplet PCR was performed using a Bio-Rad QX100 system.

**Reaction composition (24 µL total):**

| Component | Volume |
|---|---|
| ddPCR Supermix for Probes (Bio-Rad) | 12 µL |
| 20× primer/probe mix (500 nM primers, 250 nM probe) | 1.2 µL |
| DNA template | 2.4 µL |
| DNase-free water | 8.4 µL |

For droplet generation, 22 µL of reaction mixture and 72 µL of droplet generation oil were used.

**Thermal cycling conditions:**

| Step | Temperature | Duration | Cycles |
|---|---|---|---|
| Initial denaturation | 95°C | 10 min | 1 |
| Denaturation | 94°C | 30 s | 40 |
| Annealing/extension | 59°C | 1 min | 40 |
| Final incubation | 98°C | 10 min | 1 |

A ramp rate of 2°C/s was used between steps.

**Controls and analysis:**

Each plate included three no-template controls and three positive controls containing DAOM197198 DNA derived from in vitro cultures.

Droplets were read using a Bio-Rad QX100 droplet reader, and results were processed using QuantaLife software version 1.7.4. Fluorescence thresholds were determined using the [Define the Rain](https://definetherain.org.uk/) approach.

### 18S rRNA Metabarcoding

#### Target Region and Primers

- **Target region:** Small subunit (SSU) 18S rRNA gene
- **Amplification method:** Semi-nested PCR
  - **Round 1:** AML1 (forward) and AML2 (reverse) (Lee et al. 2008)
  - **Round 2:** WANDA (forward) and AML2 (reverse) (Dumbrell et al. 2011)

#### PCR Conditions

- **Round 1:** Followed Lee et al. (2008).
- **Round 2:** Followed protocols outlined by the Integrated Microbiome Resource (IMR), Dalhousie University.

#### Illumina Sequencing

- **Platform:** Illumina MiSeq
- **Sequencing chemistry:** 2 × 300 bp paired-end reads
- **Facility:** Integrated Microbiome Resource (IMR), Dalhousie University
- **Protocol:** https://imr.bio/protocols.html

## Repository Files

This repository contains the study data, sample metadata, associated data dictionaries, and analysis files.

### ddPCR Data and Metadata

**`Invasion_ddPCR_AllYears.csv`**

Original master ddPCR data covering the study sampling years. This file contains the underlying ddPCR measurements. Consult its column headings and year-specific records before analysis.

**`Invasion_ddPCR_AllYears_metadata.txt`**

Combined ddPCR dataset and metadata for 2014, 2015, 2016, 2017, 2022, and 2025. Includes soil sample measurements and separate assay-control records.

The latest version contains 559 records: 503 soil sample records and 56 assay controls.

The 2025 records include samples from Quail, KM2, KM4, and KM5, together with instrument-level measurements and control wells.

**`Invasion_ddPCR_AllYears_README.txt`**

Data dictionary and provenance documentation for the combined ddPCR dataset, including:

- Original source sheets and row identifiers
- ddPCR assay variables and units
- Soil sample and assay-control classifications
- Sampling design and transect information
- Outstanding unit and protocol verification notes

### 18S Metabarcoding Data and Metadata

**`InvasMeta_18S.txt`**

Metadata for the 12 selected 2025 soil samples used in resident AMF 18S rRNA metabarcoding analyses, together with their associated inoculant ddPCR measurements.

These 12 samples are a subset of the 2025 ddPCR sampling and do not represent additional independent soil samples.

**`InvasMeta_README_18S.txt`**

Data dictionary and methodological documentation for the 2025 metabarcoding dataset, including sample identifiers, sample selection procedures, inoculant quantification, PCR primers, Illumina MiSeq sequencing, and bioinformatic processing.

### Analysis Files

**`Community_analysis_MARKDOWN`**

Analysis of resident AMF community structure and diversity.

**`UBCO_analysis_v1_MARKDOWN`**

Statistical analyses of inoculant establishment, persistence, and abundance at the UBC Okanagan field site.

**`Knox_analysis_v2_MARKDOWN`**

Statistical analyses of inoculant establishment, persistence, and abundance at the Knox Mountain field site.

### Additional Sequencing Files

**`Invasion 2025 Sample Key for sequences`**

Sample key linking Illumina sequencing identifiers to the community metadata.

**`Feature-table-Silva-L7`**

AMF amplicon sequence variant (ASV) feature table with SILVA taxonomic assignments.

## References

Kokkoris, V., Li, Y., Hamel, C., Hanson, K., & Hart, M. (2019). Site specificity in establishment of a commercial arbuscular mycorrhizal fungal inoculant. *Science of the Total Environment*, 660, 1135–1143.

Lee, J., Lee, S., & Young, J. P. W. (2008). Improved PCR primers for the detection and identification of arbuscular mycorrhizal fungi. *FEMS Microbiology Ecology*, 65(2), 339–349.

Dumbrell, A. J., et al. (2011). Distinct seasonal assemblages of arbuscular mycorrhizal fungi revealed by massively parallel pyrosequencing. *New Phytologist*, 190(3), 794–804.

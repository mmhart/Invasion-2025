Read me

# Invasion 2025
Study tracking the persistence of a AMF (arbuscular mycorrhizal fungi) biofertlizer over 13 years. This study looks at the fate of an introduced AMF, Rhizophagus irregularis DAOM 197198 in natural systems and asks the questions:
	•	Does it establish?
	•	Does it persist over multiple growing seasons?
	•	Does it matter what the host plant was?
	•	Did uninoculated plants get colonized? When?
	•	Did inoculation affect resident AMF after 13 years?
We planted 4 plants each at 5 sites within 2 blocks near Kelowna, BC Canada in 2013.  Plants were either inoculated in the greenhouse, or in the field.  Over four years, between 2014 and 2017 we took soil cores at 10, 25 and 50cm around the plants, extracted DNA before pooling by plant and measured inoculant abundance using primers specific to DAOM197198 and ddPCR.  In 2022 and 2025, we resampled remaining plots for inoculant abundance and looked at changes to resident AMF communities using Illumina Miseq (Wanda/AML2)(Dumbrell et al. 2011, Young and Lee 2009)(2025 only).

Sample Information
	•	Land Type: Grassland
	•	Geographic Region: Okanagan, Southern British Columbia, Canada
	•	Soil Types (pooled): Loam and sandy-loam
	•	Sampling Locations:
	•	UBC Okangan: : 49.9394° N, 119.3948° W
	•	Knox Mountain: 49°54?40?N 119°28?55?W? / ?49.911°N 119.482°W
Collection Dates: 2014, 2015, 2016, 2017, 2022, 2025
DNA Extraction and PCR
	•	DNA Extraction: QIAGEN DNeasy PowerSoil Kit (0.2 g samples)
	•	PCR Cleanup Kit: QIAquick PCR Purification Kit (Qiagen Inc.)

Primer region for digitial droplet PCR (ddPCR)
PCR conditions for ddPCR

Primer Region for metabarcoding
	•	Targeted Region: SSU rRNA
	•	Semi-nested PCR:
	•	Round 1: AML1 (forward) & AML2 (reverse) - Lee et al., 2008
	•	Round 2: WANDA (forward) & AML2 (reverse) - Dumbrell et al., 2011
PCR Conditions for metabarcoding
	•	Round 1 followed Lee et al. (2008)
	•	Round 2 followed protocols outlined by the Integrated Microbiomae Resource (IMR), Dalhousie University
Sequencing
	•	Platform: Illumina MiSeq (2?300 bp)
	•	Facility: Integrated Microbiome Resource (IMR), Dalhousie University
	•	Protocol: https://imr.bio/protocols.html

Files 
The repository contains the study data, sample metadata, accompanying data dictionaries, and analysis files. The five files added in October 2026 are listed first. Names below are as shown in the GitHub project folder.
Invasion_197198_Master_Data.txt — Original master ddPCR data export covering the study sampling years. This is the source data file; consult its column headings and year-specific records before analysis.
Invasion_ddPCR_AllYears_metadata.txt — Combined ddPCR dataset/metadata for 2014, 2015, 2016, 2017, 2022 and 2025. Contains soil-sample measurements plus separate assay-control records. The latest version contains 559 rows: 503 soil records and 56 controls. The 2025 records cover Quail, KM2, KM4 and KM5; instrument-level measurements and the control wells are included.
Invasion_ddPCR_AllYears_README.txt — Data dictionary and provenance notes for the combined ddPCR file, including source sheet/row tracking, assay variables, sample/control classification, sampling design, and outstanding unit or protocol checks.
InvasMeta_18S.txt — Metadata for the twelve selected 2025 soil samples used in resident AMF 18S rRNA metabarcoding analyses, alongside their inoculant ddPCR measurements. These twelve samples are a subset of the 2025 ddPCR sampling, not twelve additional independent soil samples.
InvasMeta_README_18S.txt — Data dictionary and methods for the 2025 metabarcoding samples: identifiers, sample selection, inoculant quantification, SSU 18S primers (AML1/AML2 then WANDA/AML2), Illumina MiSeq sequencing (2 ? 300 bp), and bioinformatic processing.
Community_analysis_MARKDOWN — Analysis for resident AMF communities.
UBCO_analysis_v1_MARKDOWN — Analysis for the UBC Okanagan field site.
Knox_analysis_v2_MARKDOWN — Analysis for the Knox Mountain field site.
Invasion 2025 Sample Key for sequences — Key linking sequencing sample IDs to the community metadata.
Feature-table-Silva-L7 — AMF feature table with SILVA taxonomic assignments.
References
Lee, J., Lee, S., & Young, J. P. W. (2008). Improved PCR primers for the detection and identification of arbuscular mycorrhizal fungi. FEMS Microbiology Ecology, 65(2), 339–349.
Dumbrell, A. J., et al. (2011). Distinct seasonal assemblages of arbuscular mycorrhizal fungi revealed by massively parallel pyrosequencing. New Phytologist, 190(3), 794–804.


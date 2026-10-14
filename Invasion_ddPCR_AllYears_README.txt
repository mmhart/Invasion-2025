INVASION STUDY — COMBINED ddPCR DATA AND METADATA

Title: Long term persistence of a commercial arbuscular mycorrhizal biofertilizer in a natural grassland
Authors: Joyalea Carson-Austin; Dane Goh; Miranda M. Hart
Study: Kelowna, British Columbia, Canada; University of British Columbia Okanagan West Campus Lands (Aberdeen, Quail, Tutt) and Knox Mountain (KM2, KM4, KM5).
Target: Introduced Rhizophagus irregularis isolate DAOM197198 (197198).

FILES AND PROVENANCE
The companion Invasion_ddPCR_AllYears_metadata_UPDATED_2025.tsv combines every populated record from the five source workbook sheets T1 2014, T2 (2015), T3 (2016), T4 (2017), T5 (2022), plus the separately uploaded dd2025.csv for 2025, augmented with the subsequently supplied 2025 raw ddPCR instrument table. Each row has SourceSheet and SourceRow identifiers that permit tracing it to the original source. Numeric ddPCR assay values have not been calculated, normalized, or substituted. The workbook also includes no-call/assay control entries, which are retained. The 2025 raw ddPCR instrument table adds 3 no-template controls and 2 positive controls. The 2025 raw instrument table also includes 30 additional soil records from KM2 and KM4, incorporated as 2025 field-sample rows with their assay values preserved. Soil-normalized copies/g and presence categories for these 30 samples were not supplied and are blank. The combined table is NOT an analysis-ready replicate-aggregated dataset.

CRITICAL: This is the separate ddPCR dataset and DOES NOT assert that soil samples had no 18S sequencing; 2025 sequencing is documented in a separate metabarcoding dataset.

STUDY DESIGN (from manuscript and researcher clarification)
2014–2017: Samples collected in July using soil cores at 10, 25, and 50 cm from target plants; manuscript describes three replicate cores per distance, pooled for ddPCR, and repeated sampling across years. The 2014 source sheet also includes 1, 2, 10, 25, and 50 as distance field values; the raw codes are preserved and the 1 and 2 codes should be verified, not assumed to be cm. Dates of sampling are year only in this dataset; exact days were not provided.
2022: Following loss of some original plant markers, the manuscript describes transect-based sampling. The workbook's Transect # is retained as a POSITION INDEX, not treated as metres; physical spacing requires confirmation.
2025: The supplied dd2025.csv covers Quail and KM5. The newly supplied full instrument table ALSO includes 15 positions each for KM2 and KM4, included as additional soil samples in the master dataset. Researcher confirms Position = metres from the origin along a 30-m transect (2, 4, 6,...,30 m). THIS CORRECTS THE MANUSCRIPT, which currently describes a 10-m transect and 10 samples per plot; revise Methods to match the actual sampling scheme. The 2025 csv contains ddPCR metrics also found for some samples in the separate 18S metadata; these are not necessarily independent samples.

LAB METHODS (from manuscript)
Soil DNA extracted using Qiagen DNeasy PowerSoil kit. Isolate-specific Bio-Rad droplet digital PCR assay targeted an intergenic mitochondrial region with a specific primer/probe. Manuscript uses conflicting designations cox3-rnl and cox1-rnl; confirm before publishing the target designation. Reactions used Bio-Rad ddPCR Supermix for probes and QX100 droplet reader, QuantaLife 1.7.4; manuscript provides cycling and controls. Results include droplet counts, concentration and fluorescence amplitudes. The listed 'Concentration' values are believed to be instrument reported, but precise export units should be confirmed from original assay files before labeling as absolute copies/µL.

FIELD DEFINITIONS
RecordID: unique data-record identifier created during metadata compilation (not a biological sample ID).
CollectionYear: calendar year of sample collection.
SourceSheet: original workbook tab or 2025 CSV filename.
SourceRow: physical original spreadsheet/CSV row, with header row = 1.
RecordType: soil sample versus assay control; controls are not field replicates.
StudyLocation: broader field site inferred from plot naming convention.
Plot: subsite/plot label retained from source; earlier years are forward-filled from merged/blank grouping cells.
PlantNumber: source plant identifier carried within its plot for 2014–2017.
DistanceFromPlant_cm: 2014–2017 source Distance (interpreted in cm where standard 10/25/50 values are recorded; special 2014 values 1/2 require verification).
TransectPositionIndex: original 2022 Transect # (physical distance not assumed).
DistanceFromTransectOrigin_m: 2025 metres along 30-m transect, confirmed by researcher.
OriginalLocation, OriginalArea, OriginalPlant, OriginalDistanceOrPosition, OriginalSampleID: original source fields; blank source cells remain blank here. OriginalSampleID may contain suffix (ND), carried unchanged.
Concentration_original_reported: numerical values transcribed from the source workbook's Concentration field; UNIT UNVERIFIED, despite column label—do not treat as validated copies/µL without confirmation.
CopiesPer20uLWell: original measured/reported copy number per 20-µL well.
Positives, Negatives, AcceptedDroplets: droplet counts from source workbook.
Threshold: plate fluorescence threshold where provided, otherwise blank.
MeanAmplitudeofPositives, MeanAmplitudeofNegatives, MeanAmplitudeTotal: original fluorescence amplitude values.
CopiesPer1gSoil_reported: 2025 CSV reported copies per gram soil; conversion cannot be independently validated from provided data.
Presence_reported: original 2025 binary presence label, unmodified; no detection threshold has been imposed on other years.
Assay: droplet digital PCR.
Target: introduced R. irregularis DAOM197198.
SampleMaterial: soil for field samples, assay control otherwise.
Metabarcoding18S: indicates that 18S sequence measurements are not included in this separate ddPCR dataset; it does not make a claim about whether individual samples were sequenced elsewhere.
InstrumentWell: 2025 assay plate well, e.g., H02; not provided for 2014–2022.
InstrumentSampleNumber: 2025 raw sample number or assay-control label (NTC / 197198 +).
PoissonFractionalAbundanceMin: original 2025 raw instrument field; supplied blank for all 2025 records.
ReferenceAssayNumber, TargetAssayNumber: raw 2025 export assay identifiers (1 in supplied records); not proof of distinct targets.
PoissonConfMax68, PoissonConfMin68: 2025 raw instrument export fields named for Poisson confidence limits at 68%; preserve values as supplied. The minimum is zero across supplied records, even when abundance is positive, so meaning/validity should be confirmed before interpretation.
Notes: row-specific caveats.

2025 RAW INSTRUMENT MATCHING
2025 Quail and KM5 samples (n=30) were matched by Sample Name interpreted as plot and metre position; no conflicting CopiesPer20uLWell values were encountered. Instrument concentration and quality metrics were added as reported, without modifying soil-normalized copies/g or presence calls from dd2025.csv.
The additional KM2 (n=15) and KM4 (n=15) field samples appearing in the newly supplied 2025 instrument table were added as separate field-sample rows. Values of CopiesPer1gSoil_reported and Presence_reported remain blank because these fields were not in the supplied instrument table. If needed, confirm any site-specific exclusions for manuscript analyses before reusing all records. All five assay-control rows are included in the combined file. No template controls have zero positive droplets; positive controls have >0 positive droplets.

ROW COUNTS BY YEAR AND RECORD TYPE
2014: Assay control: 12
2014: Soil sample: 126
2015: Assay control: 12
2015: Soil sample: 90
2016: Assay control: 10
2016: Soil sample: 90
2017: Assay control: 11
2017: Soil sample: 87
2022: Assay control: 6
2022: Soil sample: 50
2025: Soil sample: 60
2025: Assay control: 5
Total: 559 records

![](fmbnlogo.png)

**FoodMicrobionet 5: Structure of edges and nodes tables**.

Prof. Eugenio Parente. Scuola di Scienze Agrarie, Forestali, Alimentari ed Ambientali, Università degli Studi della Basilicata.

This document describes the structure of data tables for FoodMicrobionet (FMBN) version 5 (v5) 

The structure of the database was changed and is now more coherent with that of a relational database, to facilitate extraction of data. It now includes:

1. **studies table**: data on the studies included in FMBN. Includes info on the target (ITS region or 16S DNA or cDNA or both, region) sequencing platform, bioinformatic pipeline, sequence accession number, number of samples, bibliographic information, primers, study type; linked to the other tables by studyID[^1]
1. **primers table**: include information on the primers used in the studies
1. **samples table** (sample nodes): data on the samples included in FMBN, see below; linked to the study table by studyId and to edge table by sampleId
1. **taxa table** (OTU nodes, formerly lineages table): data on taxa included in FMBN, see below; linked to the edge table by taxonId 
1. **edges table** (edges, the source is a sample, the target a taxon, formerly edge table): data on sample - taxon edges, with abundances; linked to the study table by studyId, the sample table by sampleId, the taxon table by taxonId; the table is available in two versions:
   1. **edges\_b**: data for bacteria
   1. **edges\_f**: data for fungi
1. **abstracts table:** for convenience only, in the future it might be used for free text searches
1. **FoodEx2 table**: for convenience only; this table includes all codes from Exposure hierarchy from revision 2
1. **Version history**: for convenience only; a table including information on the past versions of FMBN, including publication history and DOI links

The structure of the database is illustrated in Figure 1.

The emphasis is now on data extraction and processing using a range of R Shiny apps, and most of the new features were designed to facilitate processing with other software packages (Microsoft Excel, Gephi, Cytoscape, the CoNet App of Cytoscape, the bipartite package and the phyloseq package of R). To learn more on FMBN: <http://web.unibas.it/parente/?page_id=761> . 

All scripts are available at <https://github.com/ep142/FoodMicrobionet>.


![](fmbnstructure.jpeg)

**Figure 1**. Structure and relationships in FoodMicrobionet. FoodMicrobionet main tables are in dark blue. FMBN secondary tables are in turquoise. External resources are in pale blue

**The studies table**.

This table lists all relevant information for the study. The table has the following fields (each study is a record)

1. **study**: a numerical (integer) id for the study
1. **studyId**: an ID for the study (in the format STx, where x is an integer with the study number); this is assigned as studies are added to FMBN (see samples table)	
1. **FMBN\_version**: the version of FMBN in which the study appeared for the first time (unless stated otherwise, the study will appear in all further versions)	
1. **target**: the gene target(s) used for amplification	
1. **region**: the region of the 16S RNA or ITS target used for amplification
1. **platform**: the sequencing platform	
1. **read\_length\_bp**: the average read length	
1. **seq\_center**: the sequencing centre	
1. **bioinf\_software**: the bioinformatics software/pipeline used for the analysis	
1. **OTU\_picking**: the OTU picking method	
1. **assign\_tax\_method**: the taxonomy assignment method	
1. **tax\_database**: the taxonomic database(s) used for annotation	
1. **seq\_accn**: the accession number for the SRA study, if available, empty otherwise (will become NA in R)
1. **seq\_accn\_link**: the hyperlink to the run selector for the SRA study, if available, otherwise will point to a page in the FMBN website
1. **bioproject**: SRA bioproject id if available	
1. **samples**: the number of samples included in FMBN	
1. **food\_group**: 	the food group (for tabulation purposes)
1. **short\_descr**: 	a short description of the study including FoodId (see samples table)
1. **DOI**: the DOI of the study, if published, "unpublished" otherwise	
1. **DOI\_link**: the DOI hyperlink	
1. **ref\_short**: an abbreviated reference to the published paper, if available	
1. **year**: publication year; if missing, the year the experiment was carried out	
1. **ref\_complete**: a complete reference to the paper describing the study, "unpublished" otherwise
1. **corr\_author**: corresponding author surname
1. **corr\_author\_mail**: corresponding author e-mail address
1. **geoloc**: a list of the countries in which the samples were obtained; see samples table
1. **study\_type**: the study type; cross-sectional, if there are not multiple time points for each treatment/sample type, longitudinal, if there are at least 5 time points for each treatment/sample type, mixed otherwise; this is followed by two numbers in the format *mxn* where *m* is the number of treatments/sample types and *n* the number of time points[^2]; the number of replicates, if any, can be estimated from these figures and the total number of samples in the study
1. **primer\_f**: forward primer used for amplification (see primers table)
1. **primer\_r**: reverse primer used for amplification (see primers table)	
1. **overlapping**: logical flag indicating if during processing paired end sequences overlapped correctly (if not a 10N spacer was added between forward and reverse, see documentation for DADA2)	
1. **paired\_end**: logical flag indicating if paired end sequences were available[^3]

**The primers table**.

The **primers** table includes some information on the primers (when available). This can be useful if a user wants to repeat the processing of raw sequences. The following fields are included

1. **primer**: the name of the primer (must be unique)
1. **sequence**: the sequence of the primer without adapters and spacers	
1. **region**: the 16S region	
1. **start**: start nucleotide	
1. **length**: length in bp (calculated from sequence)	
1. **note:**  notes

**The samples table**.

The **samples** table (which is present in two versions, for bacteria and fungi, due to incompatibility in the information of some fields, like number of sequences, issues, etc.)[^4] includes metadata for samples/nodes. The current version of FMBN includes the following fields:

1. **studyId**: an ID for the study (in the format STx, where x is an integer with the study number), link to the same field in the studies table
1. **sampleId**: an integer with the number of the sample (samples are given a progressive number when added to FMBN), link to the same field in the edges table; the database for fungi and bacteria do not share sampleIds; however samples can be matched by biosample
1. **label\_1**: an ID for the sample (in the format SAy, where y is an integer with the number of the sample), 
1. **label\_2**: character label for the sample, as given when the sample was added to FMBN. For compatibility with older versions of FMBN, the original sample label (i.e. the label as it appears in the original study or in SRA RunInfo tables under SampleName, whichever was available when the sample was added to FMBN) may have a prefix in the format xxx\_
1. **label\_3**: this column contains the node display label (used as node display label in interactive visualisations produced with Gephi); built by combining either label\_2 and foodId (if biosample is not available) or biosample and foodId
1. **llabel**: incorporates information on foodId and food status (see below) in a 10-character format (5 characters for the foo did, 3 for food status, see below for abbreviations, 1 for target, either d for 16S RNA gene or r for 16S RNA).[^5]
1. **s\_type**: sample type: "Sample" for nodes corresponding to food samples, and "Environment" for samples from food environments.
1. **n\_reads**: number of reads after clean-up, empty if not available
1. **n\_reads2**: n\_reads if available, otherwise a dummy value (3227, the median of the values of n\_reads for samples for which this information is available in FMBN v 2) to be used for samples (food or environment) if n\_reads is not available; this value is used to rebuild OTU tables with absolute frequencies starting from the edge and samples tables
1. **n\_issues**: number of issues during sequence processing (i.e. high number of chimeras/bimeras, high loss during quality filtering, low number of features, low number of sequences after processing); set to -1 for studies 1:33 to avoid issues with column parsing when importing with readr functions
1. **Ontology fields, samples/environment**: the hierarchical terms and IDs are derived from the food classification and description system FoodEx 2 (revision 2) adopted by the European Food Safety Authority (EFSA, 2015). The terms used in FoodMicrobionet are derived from the exposure hierarchy. For food environments they are set to the same values of the food processed in that environment.[^6] 
   1. **foodId**: the code from the FoodEx 2 classification. A few special codes were used for situations that did not fit in this classification: MOCK and BLAK were used for mock communities and PCR/reagents blanks, respectively, while ENVXX was used for food environments which could not be associated to a specific food category[^7]
   1. **description**: a description of the sample, including at least terms from L6
   1. **L1**: the top level in the hierarchy
   1. **L4**: level 4 in the hierarchy
   1. **L6**: level 6 in the hierarchy.
1. **Food status fields**: these columns describe the status of the food sample (empty environmental samples). They are inferred from what is known (from the original paper, from information provided by the contributor) on a sample[^8]
   1. **nature**: **Raw** (raw material or ingredient), **Intermediate**, **Finished** (including the products during storage). Shorthands: r, i, f; 0 for environmental samples or when this information is not available.
   1. **process**: **None** (raw materials and foods other than those covered by the other options); **Mild** (foods with no physical lethal treatment; includes foods with added preservatives, mild heat treatments and fermentation); **Type2** (foods submitted at least to a lethal treatment resulting potentially in a 6D decrease in *L. monocytogenes*, but whose treatment is lower than that of type 1 foods; in these foods the heat treatment is between 70°C for 2 min and 90°C for 10 min); **Type1** (foods submitted at least to a lethal treatment resulting potentially in a 6D decrease in *C. botulinum* type E spores, approximately 90°C for 10 min or equivalent)[^9]. Shorthands: n, m, 2, 1; 0 for environmental samples or when this information is not available.
   1. **spoilage**:[^10] **NA** (no information available), **Unspoiled** (all fresh or processed, non-spoiled foods, u); **Fermented** (all unspoiled fermented foods, f), **Spoiled** (all non-fermented foods which are spoiled or are at the end of shelf-life, b), **Fermented+spoiled** (fermented foods which are also spoiled or at the end of shelf-life, b). 0 is used for environmental samples or when this information is not available.
1. **target1**: the gene used as a target
1. **target2**: the region used as a target
1. **links to NCBI SRA database**:
   1. biosample: SRA biosample, if available, otherwise empty
   1. SRA\_sample: SRA sample, if available, otherwise empty
   1. SRA\_run: SRA run, if available, otherwise empty
1. **geolocation information**: geolocation information derived from SRA metadata tables, when possible or from the published papers
   1. **geo\_loc\_country**: the country where the sample was collected	
   1. **geo\_loc\_continent**: the continent where the sample was collected	
   1. **lat\_lon**: the coordinates of the place where the sample was collected
1. **matchId**: the matching Id in the bacterial or fungal table: this is necessary because matching by biosample is not always possible

Future iterations of this table may contain further metadata on environmental variables (pH, a<sub>W</sub>, etc.) as measured or estimated values.

**The taxa table**.

The **taxa** table includes information on the lineages of taxa included in FMBN and outlinks to external resources ([NCBI taxonomy](file:///Users/eugenio/Library/CloudStorage/GoogleDrive-eugenio.parente@unibas.it/My Drive/FMBN/FMBNv4/\(NCBI taxonomy), [List of Prokaryotic Names with Standing in Nomenclature](https://lpsn.dsmz.de/)<sub>,</sub>  [Omnicrobe](https://omnicrobe.migale.inrae.fr/)). It is linked to the edges table by the taxonId field.

1. **taxonId**: an integer with the number of the taxon (taxa are given a progressive number when added to FMBN as a result of the addition of a sample-taxon relationship)
1. **label**: a label for the taxon; must be unique; corresponds to the binomial for the species (if available) or to a higher taxon (genus, family, etc. + \_unidentified) when the OTU was not identified at the species level
1. **lineage fields**: the order is not used in this version; when the information is not available the field is left empty
   1. **domain**
   1. **phylum**
   1. **class**
   1. **order**
   1. **family**
   1. **genus**
   1. **species**
1. **id\_L6** the taxonomic assignment including the lineage in the format[^11]: Root:k\_\_Bacteria\_\_Firmicutes:c\_\_Bacilli:f\_\_Lactobacillaceae:g\_\_Lactobacillus:s\_\_Lactobacillusdelbrueckii
1. **taxonomy**: a lineage built for compatibility with the CoNet app of Cytoscape (this will most likely removed in future versions, as there are much better ways of estimating microbial association networks)
1. **idelevel**: the taxonomic level of identification (domain, …, species)
1. **NCBI\_outlink**: outlink to NCBI taxonomy
1. **Secondary\_outlink**: outlink to LPSN for bacteria
1. **Omni\_habitat**: outlink to the Florilege database (taxon lives in habitat tab)
1. **Omni\_pheno**: outlink to the Florilege database (taxon exhibits phenotype tab)
1. **Omni\_use**: outlink to the Florilege database (taxon studied for use tab)

**The edges table**.

The **edges table** includes a list of all relationships between taxa and samples (food or food environment). The fields included in the current version of FoodMicrobionet are:

1. **sampleId** (Food/food environment sample): an integer with the number of the sample (samples are given a progressive number when added to FMBN), link to the same field in the samples table
1. **taxonId**: an integer with the number of the taxon, link to the same field of the taxa table
1. **weight**: the OTU abundance in the sample: for (FMBN) historical reasons this is in %[^12] but is converted in frequency (0-1) or absolute abundance during filtering and selection in the ShinyFMBN app.

**

**The Abstracts table**.

This table is included for reference purposes, since a few fields are already included in the Samples table. In the future I plan to implement a free text search in ShinyFMBN. The table is only included in FMBNplus

1. **StudyId**: the study id, link to the Study and Sample table
1. **DOI\_link**: for convenience only, taken from the Study table
1. **Abstract**: a long text field with the abstract for the publication describing the study.

**The FoodEx2 exposure table**.

This table is included for reference purposes, since a few fields are already included in the Samples table. However, users may conceivably want to use level of aggregation (L2, L3, L5) which are not included in FoodMicrobionet or pool FoodId codes under a code of higher hierarchical level.

1. **Code**: the unique identifier for a food or food category, matches FoodId in the samples table (unique)
1. **L1**: the level 1 in the classification, matches L1 in the samples table
1. **L2**: the level 1 in the classification
1. **L3**: the level 1 in the classification
1. **L4**: the level 1 in the classification, matches L4 in the Samples table
1. **L5**: the level 1 in the classification
1. **L6**: the level 1 in the classification, unique, matches L6 in the Samples table


**The version history table**.  

A very simple table with information on publication history for each version of FMBN

1. version: the version number  

1. release notes: information on features of the version

1. citation: information on the first publication for each version, if available.

**References**.

- EFSA (2015). The food classification and description system FoodEx 2 (revision 2). Supporting publications EFSA Supporting publication 2015:EN-804. http://www.efsa.europa.eu/en/supporting/pub/215e.htm
- Parente, E., Cocolin, L., De Filippis, F., Zotta, T., Ferrocino, I., O'Sullivan, O., Neviani, E., De Angelis, M., Cotter, P. D., Ercolini D. (2016). FoodMicrobionet: A database for the visualisation and exploration of food bacterial communities based on network analysis. *International Journal of Food Microbiology*, *219*, 28–37. <http://doi.org/10.1016/j.ijfoodmicro.2015.12.001>
- Parente, E., De Filippis, F., Ercolini, D., Ricciardi, A., Zotta, T., 2019. Advancing integration of data on food microbiome studies: FoodMicrobionet 3.1, a major upgrade of the FoodMicrobionet database. *International Journal of Food Microbiology*, 305, 108249. [https://doi.org/10.1016/j.ijfoodmicro.2019.108249]()
- Parente, E., Zotta, T., Ricciardi, A., 2022. FoodMicrobionet v4: A large, integrated, open and transparent database for food bacterial communities. Int J Food Microbiol 372, 109696. <https://doi.org/10.1016/j.ijfoodmicro.2022.109696>
- Parente, E., Zotta, T., Ricciardi, A., 2022. FoodMicrobionet v4: A large, integrated, open and transparent database for food bacterial communities. Int J Food Microbiol 372, 109696. https://doi.org/10.1016/j.ijfoodmicro.2022.109696
3

FMBN5Tablespecs v1.0, 19/04/2024, EP

[^1]: A study has a unique accession number; in some cases, sequences for samples belonging to the same study had been deposited with different bioprojects. In this case two entries were created in the study table but the occurrence of a second study was noted and every effort was made to manually match the samples
[^2]: If n is in the format x-y it indicates the range of time points
[^3]: This is not necessarily true for all sequences obtained by Illumina, because in some cases only forward sequences were deposited in SRA
[^4]: Due to errors in the uploading of sequences and their metadata it might indeed happen that the same samples have different biosample accession; correcting this issue is error prone; therefore we have decided to use the same description for samples which have both sequences for bacteria and for fungi only when all or part of the biosample accessions match
[^5]: this is one of the fields used for filtering and labelling (in graphs)
[^6]: there is an issue with food samples obtained during processing of, for example, a hard, ripened cheese. A cheese curd is quite obviously a completely different environment compared to a mature hard cheese. This issue is dealt with by the "Nature" field; when the correct classification is not found in FoodEx or information are partially missing the closest relative is used.
[^7]: I am also considering the addition of a field with FoodOn classification, form compatibility purposes
[^8]: in the future we may add info on intrinsic factors (pH, a<sub>W</sub>) or extrinsic factors (temperature used during storage/ripening, atmosphere); however, this information is (alas!) often missing and may need to be inferred by knowledge on the process; as an alternative one may add facets from FoodEx2, see table 22 and 25 of the document. However, this may be difficult because sample description in published studies is often incomplete.
[^9]: this is applicable to ready to eat foods. Both type1 and type2 may include fermented foods which have been submitted to a heat treatment ACMSF, 2007. Report on safe cooking of burgers. Retrieved from http://www.food.gov. uk/multimedia/pdfs/acmsfburgers0807.pdf.
[^10]: includes foods at the end of shelf life, even if they do not show clear spoilage.
[^11]: It is important to understand that the original taxonomic lineage obtained when directly analysing sequences from SRA/ENA or from contributors is manually edited for coherence. Contributors have submitted in the past OTU abundance tables with lineages in different formats and obtained from different version of taxonomic reference databases; as a result, the same taxon often had a different lineage, especially for *incertae sedis* taxa. When a new lineage is added several steps are carried out: lineages pointing to the same (at the genus, family, class, phylum level) taxon are merged; the lineage is checked against SILVA classification, bad taxa (taxa missing one level of classification) are fixed and a few genera are edited to allow linking to external database (*Corynebacterium\_1* becomes *Corynebacterium*). When any of the taxonomic terms (k\_\_; p\_\_, c\_\_, etc.) the corresponding column is left blank.
[^12]: this was needed to set correctly node size in Gephi; proportions or number of sequences can be obtained simply by dividing by 100 or by using information in the Samples table.
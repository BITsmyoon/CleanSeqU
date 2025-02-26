# CleanSeqU

**CleanSeqU** is an R package designed to effectively remove contaminants from low-biomass 16S rRNA data. It is specifically developed to decontaminate 16S rRNA sequencing data, particularly for low-biomass samples such as catheterized urine.

## Key Features

* **Innovative Decontamination Algorithm:** Implements the CleanSeqU algorithm to enhance the accuracy of 16S rRNA sequencing data analysis in low-biomass samples.
* **Characteristic Decontamination Approach:** Leverages a single blank extraction control, an in-house blacklist, and curated non-biologic taxa information for precise contaminant removal.
* **Clear Input and Output:** Uses an Amplicon Sequence Variant (ASV) count table as input and provides a decontaminated ASV count table as output.
* **User-Friendly Interface:** Offers a user-friendly interface for straightforward decontamination analysis and yields high-quality data for downstream microbiome research.

CleanSeqU is an essential tool for researchers seeking to address contamination challenges in low-biomass samples and obtain more accurate microbial community analysis results.

## Package Installation

The CleanSeqU package can be installed from GitHub.

### Install from GitHub

```R
# Install devtools package first if not installed
# install.packages("devtools")

devtools::install_github("BITsmyoon/CleanSeqU")
```
Basic Usage
The core functionality of the CleanSeqU package is provided through the **cleanseq_u_decontam()** function.

### Required Input Files

To use the `cleanseq_u_decontam()` function, the following input files are required:

* `input_asv_count`: Path to the ASV count table file
* `input_meta_data`: Path to the sample metadata file
* `input_taxa_table`: Path to the ASV taxonomy table file
* `input_taxa_blacklist`: Path to the taxonomy blacklist file
* `input_in_house_asv_blacklist`: Path to the in-house ASV blacklist file

Please refer to the function documentation (`?cleanseq_u_decontam`) or the package manual for the format of each input file.

### Example Code

The following is a simple example code to run the cleanseq_u_decontam() function.  Actual file paths should be modified according to the user environment.

#### Set file paths (Replace with your actual file paths)
asv_count_path <- "path/to/asv_count_table.tsv"
meta_data_path <- "path/to/meta_data.tsv"
taxa_table_path <- "path/to/taxa_table.tsv"
taxa_blacklist_path <- "path/to/taxa_blacklist.tsv"
in_house_asv_blacklist_path <- "path/to/in_house_asv_blacklist.tsv"

# Run CleanSeqU function
result <- cleanseq_u_decontam(
  input_asv_count = asv_count_path,
  input_meta_data = meta_data_path,
  input_taxa_table = taxa_table_path,
  input_taxa_blacklist = taxa_blacklist_path,
  input_in_house_asv_blacklist = in_house_asv_blacklist_path
)

# View Results
print(names(result))
head(result$decontamed_asv_count)
head(result$exclude_asv)
head(result$include_asv)
Detailed Explanation of Input Parameters
The cleanseq_u_decontam() function provides various input parameters for contaminant removal analysis. Detailed descriptions for each parameter are as follows:

input_asv_count: Path to the ASV count table file. The ASV count table is the basic data for microbial community analysis, containing counts (frequency) information for each ASV (Amplicon Sequence Variant) in each sample.

Column names: Should consist of asv_id column and sample ID columns.
First column: Must be the asv_id column.
Format: Tab-separated values text file (.tsv) format.
input_meta_data: Path to the sample metadata file. The sample metadata file contains sample information. In the CleanSeqU package, it is used to check negative control (NTC) sample information.

Column names: Should consist of sample_id column and ntc_check column.
sample_id column: Must match the sample names in the input_asv_count table.
ntc_check column: Negative control samples should be marked as "ntc", and other samples should be marked as "-".
Format: Tab-separated values text file (.tsv) format.
input_taxa_table: Path to the ASV taxonomy table file. The ASV taxonomy table contains taxonomic information (Taxonomy) for each ASV.

Required column names: asv_id and Taxon columns are required.
Taxon column: Taxonomic lineage information should be written in the format k_Kingdom;p_Phylum;c_Class;o_Order;f_Family;g_Genus;s_Species. Example: k_Bacteria;p_Firmicutes;c_Bacilli;o_Lactobacillales;f_Lactobacillaceae;g_Lactobacillus;s_iners
Format: Tab-separated values text file (.tsv) format.
input_taxa_blacklist: Path to the Taxonomy blacklist file. Contains blacklist information of taxa frequently detected in control samples.

Column names: Should consist of level and name columns.
level column: Taxonomic level should be written in lowercase. Example: k (kingdom), p (phylum), g (genus), etc.
name column: Taxon names to be blacklisted should be written to match the names used in the Taxon column of input_taxa_table.
Format: Tab-separated values text file (.tsv) format.
input_in_house_asv_blacklist: Path to the in-house ASV blacklist file. Contains a list of ASV IDs frequently detected in the user's experimental setting.

Column names: inhouse_asv_blacklist column is required.
inhouse_asv_blacklist column: ASV ID information to be blacklisted should be written.
Format: Tab-separated values text file (.tsv) format.
Output
The cleanseq_u_decontam() function returns a list containing three data frames:

decontaminated_asv_table: Decontaminated ASV count data frame. This is a table with counts of ASVs that are highly likely to be contaminants set to 0 by applying the contaminant removal algorithm. This table can be used for downstream analysis.
removed_contaminant_asv_ids: Data frame of removed contaminant ASV IDs. Contains ASV IDs identified as contaminants and the case (case_2, case_3, etc.) in which the ASV was removed. Useful for checking which ASVs were removed and why.
non_removed_contaminant_asv_ids: Data frame of non-removed ASV IDs (potential true signal). Contains ASV IDs that were likely contaminants but were not removed. These ASVs may be potential true biological signals, and further review may be needed in downstream analysis.
The structure and column information of each data frame can be checked through the help("cleanseq_u_decontam") command.

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

## In R

### Run CleanSeqU function

```R
library(CleanSeqU)

result <- cleanseq_u_decontam(
  input_asv_count = asv_count_path,
  input_meta_data = meta_data_path,
  input_taxa_table = taxa_table_path,
  input_taxa_blacklist = taxa_blacklist_path,
  input_in_house_asv_blacklist = in_house_asv_blacklist_path
)
```

### Required Input Files

To use the `cleanseq_u_decontam()` function, the following input files are required:

* `input_asv_count`: Path to the ASV count table file
* `input_meta_data`: Path to the sample metadata file
* `input_taxa_table`: Path to the ASV taxonomy table file
* `input_taxa_blacklist`: Path to the taxonomy blacklist file
* `input_in_house_asv_blacklist`: Path to the in-house ASV blacklist file

Please refer to the function documentation (`?cleanseq_u_decontam`) or the package manual for the format of each input file.
**Please make sure to check the structure of each input file using toy data or function documentation. (`?cleanseq_u_decontam`)**

### Example

```R
asv_count_path <-
  system.file("extdata", "toy_asv_count.txt", package = "CleanSeqU")
meta_data_path <-
  system.file("extdata", "toy_meta_data.txt", package = "CleanSeqU")
taxa_table_path <-
  system.file("extdata", "toy_taxa_table.txt", package = "CleanSeqU")
taxa_blacklist_path <-
  system.file("extdata", "toy_asv_count.txt", package = "CleanSeqU")
in_house_asv_blacklist_path <-
  system.file("extdata", "toy_taxa_blacklist.txt", package = "CleanSeqU")

result <- cleanseq_u_decontam(
  input_asv_count = asv_count_path,
  input_meta_data = meta_data_path,
  input_taxa_table = taxa_table_path,
  input_taxa_blacklist = taxa_blacklist_path,
  input_in_house_asv_blacklist = in_house_asv_blacklist_path
)
```

### Results

The `cleanseq_u_decontam()` function returns a list containing three data frames. Example outputs are as follows:

1. Names of the returned list:

```R
print(names(result))
[1] "decontamed_asv_count" "exclude_asv" "include_asv"
```

2. Head of decontamed_asv_count data frame:

```R
head(result$decontamed_asv_count)
```

| asv_id                           | sample1 | sample2 | sample3 | sample4 | sample5 | sample6 | sample7 | sample8 | sample10 | sample11 | sample12 | sample13 | sample14 | sample15 | sample16 |
|------------------------------------|---------|---------|---------|---------|---------|---------|---------|---------|----------|----------|----------|----------|----------|----------|----------|
| b9fcd7d71b74853248517b892d03a745   | 61      | 0       | 0       | 0       | 0       | 0       | 0       | 0       | 0        | 0        | 0        | 0        | 429      | 595      | 148      |
| 6ea8228cb56f8a62f932f0e613bca40e   | 61      | 0       | 0       | 0       | 0       | 0       | 0       | 0       | 0        | 0        | 0        | 0        | 0        | 0        | 0        |
| 8c8839583ce4564dbf5c6e917cafb82e   | 28      | 0       | 0       | 0       | 0       | 0       | 0       | 0       | 0        | 0        | 0        | 0        | 0        | 0        | 0        |
| a8e53bd66879cd51ecb4e7daf552accb   | 56      | 0       | 0       | 0       | 316     | 0       | 0       | 0       | 0        | 0        | 0        | 0        | 0        | 0        | 0        |
| aee0dacb66bb3b44e8986f3859fd19d6   | 52      | 0       | 0       | 0       | 0       | 17      | 0       | 0       | 0        | 0        | 0        | 0        | 0        | 0        | 0        |
| 4850a238b7fbbd5ffa204c5a5a396f2c   | 95      | 0       | 0       | 0       | 0       | 0       | 0       | 0       | 0        | 0        | 0        | 0        | 0        | 687      | 0        |

3. 

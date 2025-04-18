% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cleanseq_u_decontam.R
\name{cleanseq_u_decontam}
\alias{cleanseq_u_decontam}
\title{Contaminant Removal for Low-Biomass 16S rRNA Data}
\usage{
cleanseq_u_decontam(
  input_asv_count,
  input_meta_data,
  input_taxa_table,
  input_taxa_blacklist,
  input_in_house_asv_blacklist
)
}
\arguments{
\item{input_asv_count}{Path to the ASV count table file. The table
should have columns named 'asv_id' and sample IDs, with the first
column being 'asv_id'.}

\item{input_meta_data}{Path to the metadata file for samples. The table
should have columns named 'sample_id' and 'ntc_check'. 'sample_id'
values should match the sample names in 'input_asv_count', and
'ntc_check' should be marked as 'ntc' for negative control samples
and '-' for other samples.}

\item{input_taxa_table}{Path to the taxonomy table file for ASVs. The
table must have columns named 'asv_id' and 'Taxon'. 'Taxon' should be
in the format
'k_Bacteria;p_Firmicutes;c_Bacilli;o_Lactobacillales;f_Lactobacillaceae;g_Lactobacillus;s_iners'.}

\item{input_taxa_blacklist}{Path to the blacklist file containing taxa
frequently detected in control samples. The table should have columns
named 'level' and 'name'. 'level' should be taxonomic level in
lowercase (e.g., 'k', 'p', 'g') and 'name' should be the taxon name,
matching names in 'input_taxa_table'.}

\item{input_in_house_asv_blacklist}{Path to the blacklist file containing
ASV IDs frequently detected in user's experimental setting. The table
should have a column named 'inhouse_asv_blacklist' with ASV IDs.}
}
\value{
A list containing three data frames:
\itemize{
\item{\code{decontaminated_asv_table}:} Decontaminated ASV count data frame.
\item{\code{removed_contaminant_asv_ids}:} Data frame of removed contaminant ASV IDs.
\item{\code{non_removed_contaminant_asv_ids}:} Data frame of non-removed ASV IDs (potential true signal).
}
}
\description{
CleanSeqU is an R package designed to effectively remove contaminants
from 16S rRNA sequencing data, particularly for low-biomass samples
such as catheterized urine. This package offers:
}
\details{
\enumerate{
\item Innovative Decontamination Algorithm: Implements the CleanSeq-U
algorithm, specifically developed to enhance the accuracy of 16S rRNA
sequencing data analysis in low-biomass samples.
\item Characteristic Decontamination Approach: Leverages a single blank
extraction control, an in-house blacklist, and curated non-biologic
taxa information for precise contaminant removal.
\item Clear Input and Output: Utilizes an Amplicon Sequence Variant (ASV)
count table as input and provides a decontaminated ASV count table
as the output. This user-friendly interface allows for straightforward
decontamination analysis, yielding high-quality data ready for
downstream microbiome research.
}

CleanSeqU is an essential tool for researchers seeking to address
contamination challenges in low-biomass samples and obtain more accurate
microbial community analysis results.

Decontaminate ASV count table using CleanSeqU algorithm

This function implements the CleanSeqU decontamination algorithm to remove
contaminants from Amplicon Sequence Variant (ASV) count tables, particularly
for low-biomass samples in 16S rRNA sequencing data. It utilizes a combination
of blank extraction controls, an in-house ASV blacklist, and a curated
taxonomic blacklist to identify and remove contaminants.
}
\examples{
cleanseq_u_decontam(
 input_asv_count = "asv_count_path.tsv",
 input_meta_data = "meta_data_path.tsv",
 input_taxa_table = "taxa_table_path.tsv",
 input_taxa_blacklist = "taxa_blacklist_path.tsv",
 input_in_house_asv_blacklist = "in_house_asv_blacklist_path.tsv"
)
}

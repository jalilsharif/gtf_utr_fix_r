library(rtracklayer)
library(tidyverse)
library(data.table)
gtf_original <- rtracklayer::import("gencode.v45.annotation.gtf") %>% as.data.frame()

setDT(gtf_original)
# Filter for only 'CDS' and 'UTR'
gtf <- gtf_original[type %in% c("CDS", "UTR")]

# Function to calculate the min and max CDS locations per gene_id
cds_start_end <- function(dt) {
  cds <- dt[type == "CDS", .(CDS_Start = min(start), CDS_End = max(end)), by = transcript_id]
  return(cds)
}

# Annotate the UTRs based on their position relative to the CDS
annotate_utr <- function(dt, cds) {
  dt <- merge(dt, cds, by = "transcript_id", all.x = TRUE)
  
  dt[type == "UTR", `:=` (type = ifelse((strand == "+" & start < CDS_Start) | (strand == "-" & end > CDS_End),
                                        "five_prime_utr",
                                        ifelse((strand == "+" & end > CDS_End) | (strand == "-" & start < CDS_Start),
                                               "three_prime_utr",
                                               NA_character_)))]
  return(dt)
}

# Calculate CDS start and end
cds_info <- cds_start_end(gtf)

# Annotate UTRs
annotated_gtf <- annotate_utr(gtf, cds_info)
annotated_gtf <- annotated_gtf %>% dplyr::select(-CDS_Start,-CDS_End)
gtf_without_cds_utr <- gtf_original[!type %in% c("CDS", "UTR")]

final_gtf <- rbind(gtf_without_cds_utr,annotated_gtf)
rtracklayer::export("gencode.v45.annotation_utr_corrected.gtf")

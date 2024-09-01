library(readxl)
library(data.table)
library(tidyr)
library(ggplot2)
library(viridis)

cmp_assw <- read_xlsx("/Users/xx20081/Library/CloudStorage/OneDrive-SharedLibraries-VUMC/Shaffer, Lauren L - Phenotype Workstream/phenotyping_code_consensus/nih_cardiomyopathy_phenotyping_combined.xlsx", sheet="CMP_AssociatedWith", skip=13) |> as.data.table()
cmp_assw <- cmp_assw[Source=="ICD10" & !is.na(Concensus)]
# cmp_assw <- rbind(cmp_assw, 
#                   data.table(Code        = "I255", 
#                              Description = "Ischemic cardiomyopathy"), fill=TRUE)


ukbb <- fread("/Users/xx20081/Desktop/hf_mvmr/output/tables/ukbb_individual_codes.gz")

icd10 <- fread("/Users/xx20081/Library/CloudStorage/OneDrive-UniversityofBristol/phenotyping/ICD10_Edition5_20160401/Content/ICD10_Edition5_CodesAndTitlesAndMetadata_GB_20160401.txt")

ukbb[, dcm := grepl("I420", code)]
ukbb[, dcm_pt := any(dcm, na.rm=T), by=eid]
ukbb <- ukbb[dcm_pt==TRUE & code %in% cmp_assw$Code]

summary <- ukbb[, .N, by=.(eid,code)]
summary <- complete(summary, eid, code, fill=list(N=0)) |> as.data.table()
summary[, pct := N/sum(N), by=eid]
summary[icd10, desc := i.DESCRIPTION, on=c("code"="ALT_CODE")]
get_lvls <- summary[code!="I420"][order(-pct)]
summary[, `:=`(eid = factor(eid, levels=unique(get_lvls[,eid])), 
               desc= factor(desc, levels=c("Dilated cardiomyopathy",unique(get_lvls[,desc]))))]
summary <- summary[order(eid)]
summary[, facet := .GRP, by=eid][, facet := cut(facet, breaks=seq(0,max(facet),ceiling(max(facet)/25)))]


ggplot(summary, aes(y = desc, x = eid, fill = pct)) +
  geom_tile() +
  scale_fill_viridis(option = "plasma") +
  theme(axis.text.x = element_blank(),
        strip.text = element_blank(), 
        axis.title.y = element_blank(), 
        legend.position = "top") +
  labs(x = paste0("UKBB participant ID (n=", length(unique(summary$eid)), ")"), 
                  fill = "% of cardiomyopathy codes per individual")


  

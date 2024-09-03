library(readxl)
library(data.table)
library(tidyr)
library(ggplot2)
library(viridis)

cmp <- read_xlsx("/Users/xx20081/Library/CloudStorage/OneDrive-SharedLibraries-VUMC/Shaffer, Lauren L - Phenotype Workstream/phenotyping_code_consensus/nih_cardiomyopathy_phenotyping_combined.xlsx", sheet="CMP_AssociatedWith", skip=13) |> as.data.table()
hf  <- read_xlsx("/Users/xx20081/Library/CloudStorage/OneDrive-SharedLibraries-VUMC/Shaffer, Lauren L - Phenotype Workstream/phenotyping_code_consensus/nih_cardiomyopathy_phenotyping_combined.xlsx", sheet="HF_IsA", skip=13) |> as.data.table()
codes <- rbind(cmp, hf[Code=="I501"])[Source=="ICD10" & !is.na(Concensus)]

ukbb <- fread("/Users/xx20081/Desktop/hf_mvmr/output/tables/ukbb_individual_codes.gz")

icd10 <- fread("/Users/xx20081/Library/CloudStorage/OneDrive-UniversityofBristol/phenotyping/ICD10_Edition5_20160401/Content/ICD10_Edition5_CodesAndTitlesAndMetadata_GB_20160401.txt")

ukbb[, lvf := grepl("I420", code)]
ukbb[, lvf_pt := any(lvf, na.rm=T), by=eid]
ukbb <- ukbb[lvf_pt==TRUE & (code %in% codes$Code | grepl("lvef",code))]

summary <- ukbb[, .N, by=.(eid,code)]
summary <- complete(summary, eid, code, fill=list(N=0)) |> as.data.table()
summary[, pct := N/sum(N), by=eid]
summary[icd10, desc := i.DESCRIPTION, on=c("code"="ALT_CODE")]
get_lvls <- summary[code!="I501"][order(-pct)]
summary[, `:=`(eid = factor(eid, levels=unique(get_lvls[,eid])),
               desc= factor(desc, levels=c("Left ventricular failure",unique(get_lvls[,desc]))))]
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
                  fill = "% of codes per LVF+ve (I501) individual")


# LVEF
lvef_cohort <- ukbb[, has_lvef := any(grepl("lvef",code), na.rm = TRUE), by=eid][has_lvef==TRUE]
lvef_cohort[, any_lt50 := any(grepl("lvef",code) & value<50), by=eid]
summary_lvef <- lvef_cohort[, .(any_lt50 = any(grepl("lvef",code) & value<50),
                                lvef     = min(value[grepl("lvef",code)]),
                                num_I501 = sum(code=="I501")), by=eid]

ggplot(summary_lvef, aes(x=lvef, y=num_I501, color=any_lt50)) +
  geom_point() +
  labs(x="LVEF", y="LVF I501 code occurences", color="LVEF<50%") +
  theme(legend.position = "top")



library(readxl)
library(data.table)
library(tidyr)
library(ggplot2)
library(viridis)

dcm <- read_xlsx("/Users/xx20081/Library/CloudStorage/OneDrive-SharedLibraries-VUMC/Shaffer, Lauren L - Phenotype Workstream/phenotyping_code_consensus/nih_cardiomyopathy_phenotyping_combined.xlsx", sheet="DCM_IsA", skip=13) |> as.data.table()
mi  <- read_xlsx("/Users/xx20081/Library/CloudStorage/OneDrive-SharedLibraries-VUMC/Shaffer, Lauren L - Phenotype Workstream/phenotyping_code_consensus/nih_cardiomyopathy_phenotyping_combined.xlsx", sheet="Myocardial_infarction", skip=13) |> as.data.table()
codes <- rbind(dcm, mi)[Source=="ICD10" & !is.na(Concensus)]

ukbb <- fread("/Users/xx20081/Desktop/hf_mvmr/output/tables/ukbb_individual_codes.gz")

icd10 <- fread("/Users/xx20081/Library/CloudStorage/OneDrive-UniversityofBristol/phenotyping/ICD10_Edition5_20160401/Content/ICD10_Edition5_CodesAndTitlesAndMetadata_GB_20160401.txt")

ukbb[, dcm := grepl("I420", code)]
ukbb[, `:=`(dcm_pt   = any(dcm),
            dcm_date = min(date[code=="I420"], na.rm=TRUE)), by=eid]
ukbb <- ukbb[dcm_pt==TRUE & code %in% codes$Code & ((date <= dcm_date & code != "I420") | code == "I420")]
summary <- ukbb[, .(dcm_count = sum(code=="I420"),
                    mi_count  = sum(code!="I420"),
                    dcm_pct   = sum(code=="I420") / .N), by=eid]
summary <- summary[order(-mi_count)]
summary[, `:=`(eid = factor(eid, levels=unique(eid)))]

ggplot(summary, aes(x = eid, y = dcm_pct, group = 1)) +
  geom_point() +
  geom_line(aes(y = mi_count / max(mi_count)), color = "blue") +
  geom_hline(yintercept = 1 / (max(summary$mi_count) / max(summary$dcm_pct)), color = "red", linetype = "dashed") +
  geom_vline(xintercept = "")
  scale_y_continuous(name = "DCM Percentage",
                     sec.axis = sec_axis(~ . * max(summary$mi_count) / max(summary$dcm_pct),
                                         breaks = seq(0,max(summary$mi_count),10),
                                         labels = seq(0,max(summary$mi_count),10), name = "MI Count")) +
  labs(x = "EID") +
  theme_minimal()


ggplot(complete(summary[, .N, by=.(dcm_count, mi_count)], dcm_count, mi_count, fill=list(N=0)),
       aes(x=as.factor(dcm_count), y=as.factor(mi_count), fill=log(N))) +
  geom_tile() +
  geom_text(aes(label=N), size=3, color="white") +
  scale_fill_viridis(option="plasma", end="0.95") +
  labs(x = "DCM code count", y = "MI code count", fill = "Num patients") +
  theme(legend.text = element_blank(),
        legend.position = "top")

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
                  fill = "% of codes per DCM+ve individual")







dcm <- read_xlsx("/Users/xx20081/Library/CloudStorage/OneDrive-SharedLibraries-VUMC/Shaffer, Lauren L - Phenotype Workstream/phenotyping_code_consensus/nih_cardiomyopathy_phenotyping_combined.xlsx", sheet="DCM_IsA", skip=13) |> as.data.table()
hf  <- read_xlsx("/Users/xx20081/Library/CloudStorage/OneDrive-SharedLibraries-VUMC/Shaffer, Lauren L - Phenotype Workstream/phenotyping_code_consensus/nih_cardiomyopathy_phenotyping_combined.xlsx", sheet="HF_IsA", skip=13) |> as.data.table()

codes <- rbind(dcm, hf)[Source=="ICD10" & !is.na(Concensus)]

ukbb <- fread("/Users/xx20081/Desktop/hf_mvmr/output/tables/ukbb_individual_codes.gz")

icd10 <- fread("/Users/xx20081/Library/CloudStorage/OneDrive-UniversityofBristol/phenotyping/ICD10_Edition5_20160401/Content/ICD10_Edition5_CodesAndTitlesAndMetadata_GB_20160401.txt")

ukbb[, dcm := grepl("I420", code)]
ukbb[, `:=`(dcm_pt = any(dcm, na.rm=T),
            dcm_date = min(date[code=="I420"], na.rm=TRUE)),by=eid]
ukbb <- ukbb[dcm_pt==TRUE & code %in% codes$Code] # & ((date <= dcm_date & code != "I420") | code == "I420")]

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
       fill = "% of codes per DCM+ve individual")


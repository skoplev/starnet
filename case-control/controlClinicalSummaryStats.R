rm(list=ls())

library(plyr)

setwd("~/GoogleDrive/projects/STARNET/cross-tissue")

# Load main phenotype table
pheno = fread(
	"~/GoogleDrive/projects/STARNET/phenotype/data/current/STARNET_main_phenotype_table.2017_12_03.tsv"
)

# Exclude typo outlier
pheno[["ALAT(U/l)"]][pheno[["ALAT(U/l)"]] > 1000] = NA

# Load additional SYNTAX, Duke, lesions, and number of vessels, based on Jason's angiograph assessment
pheno_control_angio = fread(
	"~/GoogleDrive/projects/STARNET/phenotype/data/Jason/20170101 Controls results.csv"
)

# Parse control angiograph data for subsequent merging tables
pheno_control_angio = rename(pheno_control_angio,
	"SYNTAX Score"="syntax_score",
	"Duke CAD Index"="DUKE",
	"# diseased vessels"="ndv",
	"# lesions"="lesions"
)

pheno_control_angio$starnet.ID = tolower(pheno_control_angio$Angio)


# General formating and calculation
# --------------------------------
pheno[["LDL/HDL"]] = pheno[["fP-LDL-Chol(mmol/l)"]] / pheno[["fP-HDL-Chol(mmol/l)"]]

pheno$CC[pheno$CC == "x"] = NA
pheno$CC = as.integer(pheno$CC)

pheno$NYHA[pheno$NYHA == "x"] = NA
pheno$NYHA = as.integer(pheno$NYHA)


features = c("Smoker", "PresHyperlipid", "PresHypertony", "PresDiabetes", "β-blocker", "LipidLowerer", "OralAntiDiabetics", "Insulin",
	"PriorStroke", "TIA", "PriorClaudicatio", "InherHeartAttackBef60")

# Lower case and "" -> NA
for (feat in features) {
	pheno[[feat]] = tolower(pheno[[feat]])
	pheno[[feat]][pheno[[feat]] == ""] = NA  # missing value encoding
	pheno[[feat]][pheno[[feat]] == "x"] = NA  # missing value encoding
}


# Extract control subset
pheno_control = pheno[pheno$CAD.status == "control", ]


# Merge 
# TODO
pheno_control = merge(
	pheno_control,
	pheno_control_angio[, c("starnet.ID", "syntax_score", "DUKE", "ndv", "lesions")],
	by="starnet.ID",
	all=TRUE
)

# Overwrite from merged table
pheno_control$syntax_score = pheno_control$syntax_score.y
pheno_control$DUKE= pheno_control$DUKE.y
pheno_control$ndv= pheno_control$ndv.y
pheno_control$lesions= pheno_control$lesions.y


feature  = c(
	"Age", "BMI(kg/m2)", "Waist/Hip", "SBP", "DBP", "Pulse",
	"Hgb(g/l)", "WBC(109/l)", "PLT(109/l)", "Creat(μmol/l)", "CRP(mg/l)", "TSH(mU/l)", "HbA1c(%)", "bl.glucose", "ASAT(U/l)", "ALAT(U/l)", "GGT(U/l)",
	"P-Chol(mmol/l)", "fP-LDL-Chol(mmol/l)", "fP-HDL-Chol(mmol/l)", "LDL/HDL", "fP-TG(mmol/l)",
	"CC", "NYHA",
	"syntax_score", "DUKE", "ndv", "lesions"

)

summary_stats = sapply(feature, function(feat) {
	message(feat)

	norm_test = tryCatch(
		shapiro.test(pheno_control[[feat]]),
		error=function(cond) {

			return(list(p.value=NA))
		}
	)

	print(norm_test)
	c(
		mean=mean(pheno_control[[feat]], na.rm=TRUE),
		sd=sd(pheno_control[[feat]], na.rm=TRUE),
		median=median(pheno_control[[feat]], na.rm=TRUE),
		n=sum(!is.na(pheno_control[[feat]])),
		normal=norm_test$p.value > 0.01
	)
})
t(summary_stats)


feature_catagorical = c(
	"Sex", 
)


prop.table(table(pheno_control$Sex)) * 100
prop.table(table(pheno_control$Smoker)) * 100
prop.table(table(pheno_control$PresHyperlipid)) * 100
prop.table(table(pheno_control$PresHypertony)) * 100
prop.table(table(pheno_control$PresDiabetes)) * 100
prop.table(table(pheno_control[["β-blocker"]])) * 100
prop.table(table(pheno_control$LipidLowerer)) * 100
prop.table(table(pheno_control$OralAntiDiabetics)) * 100
prop.table(table(pheno_control$Insulin)) * 100

prop.table(table(pheno_control$WhatHeartDisease)) * 100  # count up multiple
prop.table(table(pheno_control$PriorStroke)) * 100
prop.table(table(pheno_control$TIA)) * 100
prop.table(table(pheno_control$PriorClaudicatio)) * 100
prop.table(table(pheno_control$InherHeartAttackBef60)) * 100


# as.data.frame(pheno_control)[, grep("Aspirin", pheno_control)]


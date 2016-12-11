# Load STARNET phenotype data
pheno = fread(file.path(
	"/Volumes/SANDY/phenotype_data",
	"STARNET_main_phenotype_table.cases.Feb_29_2016.tbl"
))

pheno$Smoking.Years = as.numeric(pheno$Smoking.Years)


y = pheno$syntax_score
# x = pheno

lm(syntax_score ~ BMI + Smoking.Years, pheno)

fit = lm(syntax_score ~ BMI + Smoking.Years + LDL + HDL + alcohol + Age, pheno)

fit = lm(syntax_score ~ BMI + Smoking.Years + LDL + HDL + alcohol + Age + PresDiabetes + PresHypertension + SBP, pheno)

summary(fit)

# lm(predicted_syntax_score ~ BMI + Smoking.Years + LDL + HDL, pheno)


plot(pheno$BMI, pheno$syntax_score)
plot(pheno$LDL, pheno$syntax_score)

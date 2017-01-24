library(qvalue)


plotRegrVolcano = function(coef, pthresh=0.05, lab_offset=0.05) {
	plot(coef[,3], -log10(coef[,4]))
	# plot(coef[,3], -log10(coef$lfdr))

	abline(h=-log10(pthresh), col=gray.colors(5)[4], lty=3)

	sig_idx = which(coef[,4] < pthresh)

	text(coef[sig_idx ,3], -log10(coef[sig_idx ,4]) + lab_offset,
		cex=0.5,
		labels=rownames(coef)[sig_idx])
}


# Load STARNET phenotype data
pheno = fread(file.path(
	"/Volumes/SANDY/phenotype_data",
	"STARNET_main_phenotype_table.cases.Feb_29_2016.tbl"
))

pheno$Smoking.Years = as.numeric(pheno$Smoking.Years)
pheno$Smoking.Years[is.na(pheno$Smoking.Years)] = 0.0

pheno = as.data.frame(pheno)

# Load Brainshake phenotype data
brainshake = fread(file.path(
	"/Volumes/SANDY/phenotype_data",
	"tz.mat"
))

brainshake_matched = brainshake[match(pheno$starnet.ID, brainshake$id), ]




y = pheno$syntax_score
# x = pheno

# lm(syntax_score ~ BMI + Smoking.Years, pheno)

fit = lm(syntax_score ~ BMI + Smoking.Years + LDL + HDL + alcohol + Age, pheno)
fit = lm(syntax_score ~ BMI + Smoking.Years + LDL + HDL + Age, pheno)
# fit = lm(syntax_score ~ BMI + LDL + HDL + Age, pheno)
# fit = lm(syntax_score ~ BMI + Smoking.Years + LDL + HDL + alcohol + Age + PresDiabetes + PresHypertension + SBP, pheno)
# fit = lm(DUKE ~ BMI + Smoking.Years + LDL + HDL + alcohol + Age, pheno)
# fit = lm(lesions ~ BMI + Smoking.Years + LDL + HDL + alcohol + Age, pheno)

summary(fit)

coef = summary(fit)$coefficients

plotRegrVolcano(coef, 0.1, 0.02)

# lm(predicted_syntax_score ~ BMI + Smoking.Years + LDL + HDL, pheno)


plot(pheno$BMI, pheno$syntax_score)
plot(pheno$LDL, pheno$syntax_score)


fit = lm(pheno$syntax_score, brainshake_matched)


# All brainshake regression on SYNTAX score
# --------------------------------------------------------------

x = brainshake[, 2:(ncol(brainshake) - 3)]
x = scale(x)
x[is.na(x)] = 0
x = as.data.frame(x)

fit = lm(syntax_score ~ ., brainshake[, 2:(ncol(brainshake) - 3)])

# fit = lm(syntax_score ~ ., x)

# summary(fit)
coef = summary(fit)$coefficients
coef = as.data.frame(coef)

coef$lfdr = qvalue(coef[,4])$lfdr

plotRegrVolcano(coef)


plot(brainshake$ASAT, brainshake$syntax_score)
plot(brainshake$S_VLDL_CE__perc, brainshake$syntax_score)



# # All phenotype regression analysis
# # ------------------------------------------------------------
# include_features = apply(pheno, 2, function(feature) length(unique(feature))) > 3
# include_features[names(include_features) %in% c("starnet.ID", "DUKE", "lesions", "ndv")] = FALSE

# # include_features = include_features[!names(include_features) %in% c("starnet.ID")]


# # fit = lm(syntax_score ~ ., pheno[, 2:ncol(pheno)])

# fit = lm(syntax_score ~ ., pheno[, include_features])
# coef = summary(fit)$coefficients

# plotRegrVolcano(coef, 0.2)
# # plot(coef[, 3], -lof10)
# # rownames(coef)

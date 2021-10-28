# Load clinical data

rm(list=ls())

library(data.table)
library(readxl)
library(RColorBrewer)
library(GA)  # genetic algorithms

library(compiler)
enableJIT(3)

# Plots estimated 1D densities of 3 set of numbers.
# z, y is considered the limit, and a grey area is rendered in between.
plotCmpDensityCurves = function(x, y, z, xlim, cols, ...) {
	dens_x = density(x,
		na.rm=TRUE,
		from=xlim[1], to=xlim[2])
	dens_y = density(y,
		na.rm=TRUE,
		from=xlim[1], to=xlim[2])
	dens_z = density(z,
		na.rm=TRUE,
		from=xlim[1], to=xlim[2])

	# Set up viewport
	plot(
		dens_x,
		type="n",
		xlim=xlim,
		main="",
		...
	)

	# shade area between x and y
	polygon(
		c(dens_x$x, rev(dens_y$x)),
		c(dens_x$y, rev(dens_y$y)),
		col=rgb(0.9, 0.9, 0.9),
		border=NA
	)

	# Plot density curves
	lines(dens_x, col=cols[1], lwd=1.5)
	lines(dens_y, col=cols[2], lwd=1.5)
	lines(dens_z, col=cols[3], lwd=1.5)
}

setwd("/Users/sk/Google Drive/projects/STARNET-controls")


# Load data
# -------------------------------------------------------
# Read encrypted covariate and phenotype data.
# Does not currently contain clinical data on all controls
starnet_info = fread(file.path(
	"/Volumes/SANDY/phenotype_data",
	"covariates.cases_and_controls.April_12_2016.txt"
))

starnet_info$is_ctrl = sapply(starnet_info$subject, function(str) substr(str, 1, 1) == "C")

# Phenotype data for STARNET cases
pheno = fread(file.path(
	"/Volumes/SANDY/phenotype_data",
	"STARNET_main_phenotype_table.cases.Feb_29_2016.tbl"
))

# Load RNA control availability data
controls_avail = read_excel("control_info/CONTROL RNA C1 TO C250 OVERALL.xlsx")
controls_avail = controls_avail[1:250,]

# Trim whitespace
controls_avail[["ID NO"]] = trimws(controls_avail[["ID NO"]])

# Load additional clinical information on.
# Recieved from Raili Ermel
control_info = read_excel("control_info/c1 to c250.xlsx")


# Get the variables of interest

# STARNET population
starnet = data.frame(
	id=pheno$starnet.ID,
	bmi=pheno$BMI,
	age=pheno$Age,
	gender=pheno$Sex
)
# Exclude STARNET IDs 1-100. Assumes complete and ordered
starnet = starnet[101:nrow(starnet), ]

# Controls
ctrl = data.frame(
	id=control_info[["starnet-ID"]],
	age=control_info$Age,
	bmi=control_info[["BMI(kg/m2)"]],
	gender=tolower(control_info$Sex)
)

# Run genetic algorithm to optimize selection of STARNET samples
# Tries to balance BMI, age, and gender distributions.
k_opt = nrow(ctrl)  # optimal number of selected patients from STARNET
size_pen = 0.0005  # penalty for selection sizes different from k_opt

# Fitness function to optimize
# sel is a binary vector containing values (0, 1)
f = function(sel) {
	sel_bool = sel == 1  # selected starnet patients, boolean vector
	k = sum(sel_bool)  # size of selection

	# Kolmogorov-Smirnoff estimate for distribution similarity
	ks_age = ks.test(ctrl$age, starnet$age[sel_bool])
	ks_bmi = ks.test(ctrl$bmi, starnet$bmi[sel_bool])

	# Penalty terms
	D_age = ks_age$statistic
	D_bmi = ks_bmi$statistic
	D_gender = abs(
		# Frequency difference for gender 1
		prop.table(table(ctrl$gender))[1] -
		prop.table(table(starnet$gender[sel_bool]))[1]
	)
	D_sel = size_pen * abs(k - k_opt)  # selection size

	fitness = -D_age - D_bmi - D_gender - D_sel
	return(fitness)
}

# Run optimization
sol = ga(type="binary",
	fitness=f,
	nBits=nrow(starnet),
	maxiter=2000,
	optim=TRUE  # local optimization (Hybrid GA)
)

save(sol, file="solutions/without100.RData")


# Analyze output
# -------------------------------------------
# Multiple solutions, test length
# apply(sol@solution, 1, sum)

# Matched selection
starnet_sel = which(sol@solution[1,] == 1)  # first solution returned

length(starnet_set)

# List starnet IDs of selection
starnet$id[starnet_sel]

all_id_ordered = c(as.character(starnet$id[starnet_sel]), as.character(ctrl$id))

# Randomize selected samples and control patient IDs
all_id = sample(all_id_ordered)

write.table(all_id, file="rand_ids/rand_sample_ctrls.txt",
	row.names=FALSE, col.names=FALSE, quote=FALSE )

write.table(all_id_ordered, file="rand_ids/sample_ctrls.txt",
	row.names=FALSE, col.names=FALSE, quote=FALSE )


# Plot of opimization objective function
pdf("distributions/genetic_optimization.pdf")
plot(sol)
dev.off()

red = brewer.pal(9, "Set1")[1]
blue = brewer.pal(9, "Set1")[2]

# Plot of age and BMI distributions
pdf("distributions/age_bmi.pdf", height=4, width=4)
par(mfrow=c(2, 1), mar=c(4, 4, 1, 1))

plotCmpDensityCurves(
	starnet$age,
	ctrl$age,
	starnet$age[starnet_sel],
	xlim=c(10, 100),
	xlab="Age (years)",
	cols=c(red, "black", blue)
)

plotCmpDensityCurves(
	starnet$bmi,
	ctrl$bmi,
	starnet$bmi[starnet_sel],
	xlim=c(10, 50),
	cols=c(red, "black", blue),
	xlab="BMI",
	ylim=c(0, 0.09)
)

legend("topright", legend=c("STARNET", "Match", "Ctrl"), col=c(red, blue, "black"),
	pch=15, cex=0.8)
dev.off()

# Barplots of gender distributions
pdf("distributions/gender.pdf", height=2.0, width=3.5)
par(mfrow=c(1, 3))
barplot(table(starnet$gender),
	main="STARNET",
	las=2,
	col=red)
barplot(table(starnet$gender[starnet_sel]),
	main="Match",
	las=2,
	ylab="Individuals",
	col=blue)
barplot(table(ctrl$gender),
	main="Ctrl",
	las=2,
)
dev.off()

library(data.table)

rm(list=ls())


setwd("~/GoogleDrive/projects/STARNET/cross-tissue")

kd = fread("co-expression/annotate/grn_vamsi_eqtl/kda/modules.results.txt")

kd$MODULE = factor(kd$MODULE, levels=1:224)

df = data.frame(table(kd$MODULE[kd$FDR < 0.05]))
colnames(df) = c("ID", "n_key_drivers")

write.csv(df,
	row.names=FALSE,
	file="co-expression/annotate/grn_vamsi_eqtl/kd_per_module.csv")


rm(list=ls())
#Directed
kd = fread("co-expression/annotate/grn_vamsi_eqtl/kda/modules.directed.results.txt")

kd$MODULE = factor(kd$MODULE, levels=1:224)

df = data.frame(table(kd$MODULE[kd$FDR < 0.05]))
colnames(df) = c("ID", "n_key_drivers")

write.csv(df,
	row.names=FALSE,
	file="co-expression/annotate/grn_vamsi_eqtl/kd_per_module_directed.csv")


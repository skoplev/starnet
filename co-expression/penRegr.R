
# Regression models
library(glmnet)

mat[is.na(mat)] = 0

i = 2
fit = glmnet(
	t(mat[1:nrow(mat) != i,]),
	t(mat[i,, drop=FALSE])
)

lambda = 0.1
sum(coef(fit, s=lambda) > 0)

coef(fit, s=lambda)

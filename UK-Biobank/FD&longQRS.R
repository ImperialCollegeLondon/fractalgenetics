# UKBB data
# Logistic regression of mean global FD as a predictor of ECG QRS interval > 120ms

longQRS<- which(dataI$ECG_QRS_DURATION>120)
dataI$longQRS_binary<- rep(0, nrow(dataI))
dataI$longQRS_binary[longQRS]<- 1
fit<- glm(longQRS_binary ~ scale(MEANGLOBALFD), data=dataI, family=binomial(link='logit'))

exp(fit$coefficients)
exp(confint(fit))

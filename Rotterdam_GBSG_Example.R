####################################
# Rotterdam Tumor Bank
####################################
library(data.table)     # Store and process data
library(ggplot2)        # Plot

library(simsurv)        # Simulate data

library(discSurv)       # reshape data
library(riskRegression) # measure performance
library(parallel)       # speed up the bootstrap

library(survival)       # Models
library(glmnet)
library(gbm)
library(glinternet)
library(randomForestSRC)
library(ranger)
library(ANN2)

# Change the path to the "Helpers" file on your machine!
source("Helpers.R")

################################################################################
# Helpers
################################################################################
SAVE_RESULTS_HERE = "~/Desktop/rotterdam/"
path <- function(model) paste0(SAVE_RESULTS_HERE, "results_", model, ".csv") 
saveme <- function(dt, model) {
  dt[, modelname := model]
  write.csv(dt, path(model), row.names = F)
}

################################################################################
# Data processing
################################################################################
# File downloaded here:
# https://www.stata-press.com/data/fpsaus.html
# linked from the paper:
### Royston, Patrick, and Douglas G. Altman. 
### "External validation of a Cox prognostic model: principles and methods." 
### BMC medical research methodology 13 (2013): 1-15.
rd <- foreign::read.dta('~/Desktop/rott2.dta')
rd <- data.table(rd)

# This is from the survival package - the German Breast Cancer Study Group:
gb <- data.table(gbsg)

# Subset to node-positive:
rd <- rd[nodes > 0, ]
gb <- gb[nodes > 0,]

# German data from days to months, and make the column names match the Rotterdam data:
gb[, event.time := rfstime/30.437] 
setnames(gb, "pgr", "pr")
setnames(gb, "status", "event")
setnames(gb, "pid", "id")
# median survival: 21 months
# gb[status == 1, median(event.time)]


# Compute survival columns. Here are the column namess:
# rf - relapse free interval measured in months.
# rfi - relapse indicator.
# m - metastasis free.
# mfi - metastasis status.
# os - overall survival
# osi - overall survival indicator
rd[osi == "deceased", event.time := os]
rd[mfi == "yes", event.time := mf]
rd[rfi == 1, event.time := rf]
rd[, event := 1 * ( (event.time < 84) & ( (rfi == 1) | (osi == "deceased") | (mfi == "yes") ) )]
rd[event == 0, event.time := pmin(84, rf, os, mf)]

rd[, hormon := ifelse(hormon == 'yes', 1, 0)]
rd[, meno   := ifelse(meno == 'post', 1, 0)]

# Discretize time:
rd[, event.time.rounded := as.integer(100*round(event.time, 1))]
gb[, event.time.rounded := as.integer(100*round(event.time, 1))]

############################################################
# Do survival stacking
############################################################
# use discSurv to get the long form of our data.
cols <- c("age", "grade", "nodes", "pr", "er", "meno", "hormon")
rdLong <- dataLong(dataShort = data.frame(rd), 
                   timeColumn = "event.time.rounded", 
                   eventColumn = "event")
# discSurv creates a risk set for every time observed in the data. 
# But, we only want risk sets for every EVENT time observed. 
# So, we keep risk sets corresponding only to those event times.
keep.times <- sort(unique(rd[event == 1, event.time.rounded]))
rdLong <- rdLong[rdLong$timeInt %in% keep.times, ]
rdLong <- rdLong[, c("y", "age", "grade", "nodes", "pr", "er", "meno", "hormon", "timeInt")]
colnames(rdLong) <- gsub("timeInt", "time", colnames(rdLong))

# Now, get the long form of test data.
# We can't use discSurv, because we want to have risk sets for each 
# event time observed in the training data.
gbLong <- cbind(gb[, c(cols, "id"), with=F], time = keep.times[1])
for(tt in keep.times[2:length(keep.times)]){
  gbLong <- rbind(gbLong, 
                  cbind(gb[, c(cols, "id"), with=F], time = tt)
  )
}

X      <- rdLong[, c(cols, "time")]
X.test <- gbLong[, c(cols, "time"), with=F]

# Finally, let's standardize our data:
cms <- colMeans(X)
sds <- apply(X, 2, sd)
X.scaled <- sweep(X, 2, cms, FUN = "-")
X.scaled <- sweep(X.scaled, 2, sds, FUN = "/")
y <- rdLong$y

X.scaled.test <- sweep(X.test, 2, cms, FUN = "-")
X.scaled.test <- sweep(X.scaled.test, 2, sds, FUN = "/")

# Some models will require matrices, some want data frames:
df   <- data.frame(X.scaled)
df$y <- y

X.test.df <- data.frame(X.test)

# standardize for survival methods
rd.scaled <- as.matrix(rd[, cols, with=F])
gb.scaled <- as.matrix(gb[, cols, with=F])
cms <- colMeans(rd.scaled)
sds <- apply(rd.scaled, 2, sd)
rd.scaled <- sweep(rd.scaled, 2, cms, FUN = "-")
rd.scaled <- sweep(rd.scaled, 2, sds, FUN = "/")
gb.scaled <- sweep(gb.scaled, 2, cms, FUN = "-")
gb.scaled <- sweep(gb.scaled, 2, sds, FUN = "/")

gb.y      <- Surv(gb$event.time.rounded, gb$event)

# We will need this to estimate survival probabilities for each subject:
all.ids <- unique(gbLong$id)

f################################################################################
# Cox
################################################################################
cph <- cv.glmnet(
  as.matrix(rd[, .(age, grade, nodes, pr, er, meno, hormon)]),
  Surv(rd$event.time.rounded, rd$event),
  family = "cox"
)

baseline <- predict(cph, as.matrix(rd[, cols, with=F]), type = "response", s="lambda.1se")
baseline <- sapply(keep.times, function(t) 1/sum(baseline[rd$event.time.rounded >= t]))

# Predict with test data:
preds <- predict(cph, as.matrix(gb[, cols, with=F]), type='response', s="lambda.1se")
stacked.coxl1.pred.matrix <- t(sapply(preds, function(pred) exp(-cumsum(baseline) * pred)))

ci.coxl1    <- get.bootstrap.ci(gb.y, stacked.coxl1.pred.matrix, keep.times)
point.coxl1 <- get.one.ibs(gb.y, stacked.coxl1.pred.matrix, keep.times)

saveme(data.table(bs = ci.coxl1, point = point.coxl1), "Cox_Regularized")

################################################################################
# COX GBM
################################################################################
model  <- gbm(Surv(event.time.rounded, event)~., 
              data = cbind(rd.scaled, rd[, .(event, event.time.rounded)]),
              n.trees = 20,  
              shrinkage = .01,
              distribution="coxph",
              interaction.depth = 2);

# Compute the baseline hazard for the model (using train data)
baseline <- exp(predict(model, data.frame(rd.scaled)))
baseline <- sapply(keep.times, function(t) 1/sum(baseline[rd$event.time.rounded >= t]))

# Predict with test data:
preds <- predict(model, data.frame(gb.scaled))
preds <- c(sapply(preds, function(pred) exp(-cumsum(baseline) * exp(pred))))
pred.matrix <- matrix(preds, ncol=length(keep.times), nrow=length(all.ids), byrow = T)

ci.gbm    <- get.bootstrap.ci(gb.y, pred.matrix, keep.times)
point.gbm <- get.one.ibs(gb.y, pred.matrix, keep.times)

saveme(data.table(bs = ci.gbm, point = point.gbm), "GBM")

################################################################################
# random survival forest
################################################################################
rfmodel <- rfsrc(Surv(event.time.rounded, event)~., 
                 data = rd[, c("event", "event.time.rounded", cols), with=F], 
                 ntree = 500,
                 ntime = length(keep.times))
preds <- predict(rfmodel, gb[, cols, with=F])$survival

ci.rf    <- get.bootstrap.ci(gb.y, preds, keep.times)
point.rf <- get.one.ibs(gb.y, preds, keep.times)

saveme(data.table(bs = ci.rf, point = point.rf), "RF")

################################################################################
# stacked GBM
################################################################################
stacked.gbm <- gbm.fit(X.scaled,
                       y,
                       n.trees=10000,
                       verbose = T,
                       distribution="bernoulli")

preds   <- predict(stacked.gbm, X.scaled.test)
preds   <- 1/(1 + exp(-preds))

stacked.gbm.pred.matrix <- group.preds(preds, gbLong$id, all.ids, keep.times)

ci.gbm.st    <- get.bootstrap.ci(gb.y, stacked.gbm.pred.matrix, keep.times)
point.gbm.st <- get.one.ibs(gb.y, stacked.gbm.pred.matrix, keep.times)

saveme(data.table(bs = ci.gbm.st, point = point.gbm.st), "GBM_Stacked")

################################################################################
# stacked rf
################################################################################
stacked.rf <- ranger(y ~ ., data = df,
                     num.trees=100,
                     max.depth = 3,
                     mtry=5)

preds   <- predict(stacked.rf, X.test.df)$predictions

stacked.rf.pred.matrix <- group.preds(preds, gbLong$id, all.ids, keep.times)
ci.rf.st    <- get.bootstrap.ci(gb.y, stacked.rf.pred.matrix, keep.times)
point.rf.st <- get.one.ibs(gb.y, stacked.rf.pred.matrix, keep.times)

saveme(data.table(bs = ci.rf.st, point = point.rf.st), "RF_Stacked")

################################################################################
# stacked glinternet
################################################################################
# VERY SLOW!
# ~30 minutes on my MacBook Pro
start = proc.time()
stacked.gl <- glinternet(X.scaled, y,
                         numLevels = rep(1, 8),
                         interactionCandidates = 1:8,
                         family = "binomial"
)
proc.time() - start
preds   <- predict(stacked.gl, X.scaled.test, type="response")
preds   <- preds[, 50] # least regularization

stacked.glint.pred.matrix <- group.preds(preds, gbLong$id, all.ids, keep.times)
ci.glint    <- get.bootstrap.ci(gb.y, stacked.glint.pred.matrix, keep.times)
point.glint <- get.one.ibs(gb.y, stacked.glint.pred.matrix, keep.times)

saveme(data.table(bs = ci.glint, point = point.glint), "glinternet_Stacked")

################################################################################
# stacked nn
################################################################################
nn <- neuralnetwork(X.scaled,
                    y,
                    hidden=c(100),
                    activ.functions = "tanh",
                    loss.type="log",
                    regression=F,
                    batch.size=2000,
                    optim.type = 'adam',
                    verbose=T,
                    n.epochs=50,
                    standardize=F)

preds   <- predict(nn, newdata = X.scaled.test)
preds   <- preds$probabilities[, 2]

stacked.preds <- group.preds(preds, gbLong$id, all.ids, keep.times)
ci.st.nn     <- get.bootstrap.ci(gb.y, stacked.preds, keep.times)
point.st.nn  <- get.one.ibs(gb.y, stacked.preds, keep.times)

saveme(data.table(bs = ci.st.nn, point = point.st.nn), "neuralnet_Stacked")


################################################################################
# NNet-survival
# done in Python
################################################################################
write.csv(rd, "~/Desktop/forpython_rd.csv", row.names=F)
write.csv(gb, "~/Desktop/forpython_gb.csv", row.names=F)

# # Packages
# import sys
# # https://github.com/MGensheimer/nnet-survival cloned to my machine
# sys.path.insert(0, './nnet-survival/')
# 
# import numpy as np
# import matplotlib.pyplot as plt
# import matplotlib
# from tensorflow.keras.models import Sequential, Model
# from tensorflow.keras import optimizers, layers, regularizers
# from tensorflow.keras.callbacks import EarlyStopping, ModelCheckpoint
# from tensorflow.keras.layers import Input, Dense, Activation
# from tensorflow.keras.models import load_model
# from lifelines import KaplanMeierFitter
# from lifelines import CoxPHFitter
# from lifelines.utils import concordance_index
# import nnet_survival
# 
# from sklearn.preprocessing import StandardScaler
# 
# import pandas as pd
# 
# # Data
# rd = pd.read_csv("/Users/erincraig/Desktop/forpython_rd.csv")
# gb = pd.read_csv("/Users/erincraig/Desktop/forpython_gb.csv")
# 
# breaks = rd[rd.event == 1]["event.time.rounded"].unique()
# breaks.sort()
# 
# n_intervals = len(breaks) - 1
# 
# y = nnet_survival.make_surv_array(rd["event.time.rounded"],
#                                   rd["event"],
#                                   breaks)
# 
# x = rd[["age", "grade", "nodes", "pr", "er", "meno", "hormon"]].values
# 
# ss = StandardScaler()
# xx = ss.fit_transform(X = x)
# 
# xtest = ss.transform(X = gb[["age", "grade", "nodes", "pr", "er", "meno", "hormon"]].values)
# 
# # Model
# model = Sequential()
# model.add(Dense(100, input_shape=(7,)))
# model.add(Dense(n_intervals))
# model.add(Activation('sigmoid'))
# model.compile(loss=nnet_survival.surv_likelihood(n_intervals), 
#               optimizer=optimizers.RMSprop(learning_rate = .001))
# early_stopping = EarlyStopping(monitor='loss', patience=2)
# 
# history=model.fit(xx, y, batch_size=32, epochs=1000, verbose = 2, callbacks=[early_stopping])
# y_pred=model.predict_proba(x, verbose=0)
# 
# 
# # Predict
# predicted = model(xtest)
# 
# predicted = np.concatenate(
#   (
#     np.ones(predicted.shape[0]).reshape(-1, 1),
#     predicted),
#   1
# )
# 
# predicted = np.cumprod(predicted, 1)
# 
# pd.DataFrame(predicted).to_csv("/Users/erincraig/Desktop/nnet_survival_predicted_for_R", index=False)

preds <- fread("/Users/erincraig/Desktop/nnet_survival_predicted_for_R")
preds <- preds[2:nrow(preds), ]

ci.nnet    <- get.bootstrap.ci(gb, as.matrix(preds), keep.times)
point.nnet <- get.one.ibs(gb, as.matrix(preds), keep.times)
point.nnet

nnet.ci <- data.table(bs = ci.nnet, point = point.nnet)
saveme(nnet.ci, "nnet_survival")

################################################################################
# plot
################################################################################
all <- fread(paste0(SAVE_RESULTS_HERE, "results_Cox_Regularized.csv"))
all <- rbind(all, fread(paste0(SAVE_RESULTS_HERE, "results_GBM.csv")))
all <- rbind(all, fread(paste0(SAVE_RESULTS_HERE, "results_nnet_survival.csv")))
all <- rbind(all, fread(paste0(SAVE_RESULTS_HERE, "results_RF.csv")))
all <- rbind(all, fread(paste0(SAVE_RESULTS_HERE, "results_GBM_Stacked.csv")))
all <- rbind(all, fread(paste0(SAVE_RESULTS_HERE, "results_glinternet_Stacked.csv")))
all <- rbind(all, fread(paste0(SAVE_RESULTS_HERE, "results_RF_Stacked.csv")))
all <- rbind(all, fread(paste0(SAVE_RESULTS_HERE, "results_neuralnet_Stacked.csv")))

all[, data := "original"]
all[grepl("Stacked", modelname), data := "stacked"]
all[, modelname.processed := modelname]
all[, modelname.processed := gsub("_Stacked", "", modelname)]
all[modelname == "Cox_Regularized", modelname.processed := "Cox (LASSO)"]
all[modelname == "Cox_Regularized", modelname.processed := "Cox (LASSO)"]

alpha = .05
all.summary <- all[, .(lo = quantile(bs, alpha/2), 
                       hi = quantile(bs, 1 - (alpha/2)),
                       mn = mean(bs)), 
                   by=.(point, data, modelname, modelname.processed)]
all.summary$modelname.factor <- factor(all.summary$modelname, 
                                       levels = c(
                                         "GBM",
                                         "Cox_Regularized",
                                         "nnet_survival",
                                         "RF",
                                         "neuralnet_Stacked",
                                         "RF_Stacked",
                                         "glinternet_Stacked",
                                         "GBM_Stacked"
                                       )
)

ggplot(all.summary, aes(x=modelname.factor, y=mn, ymin=lo, ymax=hi, color=data)) + 
  geom_errorbar(width=.2) +
  geom_point(size=1.5, shape=19) +
  scale_x_discrete(labels=c("Cox_Regularized" = "Cox PH\n(LASSO)",
                            "RF" = "Random\nsurvival\nforest",
                            "GBM" = "GBM\nCox",
                            "nnet_survival" = "NNet\nSurvival",
                            "neuralnet_Stacked" = "Neural\nnet",
                            "GBM_Stacked" = "GBM", "RF_Stacked" = "RF", 
                            "glinternet_Stacked" = "glinternet")) +
  theme_minimal() + 
  labs(x="", y="Integrated Brier Score", color="Data",
       title = "Method comparison: Rotterdam & GBSG data",
       subtitle="(with 95% stratified bootstrap confidence intervals)",
       caption="(lower scores are preferred)") +
  scale_color_manual(breaks=c("original", "stacked"), values=c("#820000", "#F4795B")) #+
ggsave(paste0(SAVE_RESULTS_HERE, "Figure1.pdf"), width=6, height=3)


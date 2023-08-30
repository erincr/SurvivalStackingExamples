####################################
# Simulations
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
SAVE_RESULTS_HERE = "~/Desktop/simulated/"
path <- function(model) paste0(SAVE_RESULTS_HERE, "results_", model, ".csv") 
saveme <- function(dt, model) {
  dt[, modelname := model]
  write.csv(dt, path(model), row.names = F)
}

################################################################################
# Simulate data
################################################################################
# mostly from the simsurv examples.
ntrain = 2000
ntest  = 1000
ntot   = ntrain + ntest
p      = 5
covs   = matrix(rnorm((ntrain + ntest) * p), nrow = ntot, ncol = p)

X <- data.frame(id = 1:ntot, covs)
betas = c(-0.08,  -0.06, 0.02,  0, 0)
names(betas) = names(X)[2:ncol(X)]
y <- simsurv(lambdas = 1, gammas = 1.5, betas = betas,
             x = X, tde = c(X1 = 5),
             maxt = 1)
mean(y$status)

y$eventtime <- round(100*y$eventtime) # discSurv requires integer value times

X.train <- X[X$id <= ntrain, ]
y.train <- y[y$id <= ntrain, ]

X.test <- X[X$id > ntrain, ]
y.test <- y[y$id > ntrain, ]

df.train <- merge(X.train, y.train)
df.train <- df.train[, colnames(df.train)[2:ncol(df.train)]]
df.test  <- merge(X.test,  y.test)
df.test  <- df.test[, colnames(df.test)[2:ncol(df.test)]]

############################################################
# Do survival stacking! 
############################################################
# use discSurv to get the long form of our data.
cols = colnames(df.train)[1:p]
train.Long <- dataLong(dataShort = df.train, 
                       timeColumn = "eventtime", 
                       eventColumn = "status")
# discSurv creates a risk set for every time observed in the data. 
# But, we only want risk sets for every EVENT time observed. 
# So, we keep risk sets corresponding only to those event times.
keep.times <- sort(unique(df.train[df.train$status == 1, "eventtime"]))
train.Long <- train.Long[train.Long$timeInt %in% keep.times, ]
train.Long <- train.Long[, c("y", cols, "timeInt")]
colnames(train.Long) <- gsub("timeInt", "time", colnames(train.Long))


# Now, get the long form of test data.
# We can't use discSurv, because we want to have risk sets for each 
# event time observed in the training data.
test.Long <- cbind(X.test[, c(cols, "id")], time = keep.times[1])
for(tt in keep.times[2:length(keep.times)]){
  test.Long <- rbind(test.Long, 
                     cbind(X.test[, c(cols, "id")], time = tt)
  )
}

# We will need this to estimate survival probabilities for each subject:
all.ids <- (ntrain+1):ntot
id.ix   <- test.Long$id

############################################################
# Standardize
############################################################
# Standardize stacked matrices
cms <- colMeans(train.Long[, c(cols, "time")])
sds <- apply(train.Long[, c(cols, "time")], 2, sd)
long.scaled <- sweep(train.Long[, c(cols, "time")], 2, cms, FUN = "-")
long.scaled <- sweep(long.scaled, 2, sds, FUN = "/")
long.y <- train.Long$y

long.scaled.test <- sweep(test.Long[, c(cols, "time")], 2, cms, FUN = "-")
long.scaled.test <- sweep(long.scaled.test, 2, sds, FUN = "/")

# standardize for survival methods
short.scaled <- as.matrix(df.train[, cols])
short.scaled.test <- as.matrix(df.test[, cols])
cms <- colMeans(short.scaled)
sds <- apply(short.scaled, 2, sd)
short.scaled <- sweep(short.scaled, 2, cms, FUN = "-")
short.scaled <- sweep(short.scaled, 2, sds, FUN = "/")
short.y <- Surv(df.train$eventtime, df.train$status)

short.scaled.test <- sweep(short.scaled.test, 2, cms, FUN = "-")
short.scaled.test <- sweep(short.scaled.test, 2, sds, FUN = "/")
short.y.test <- Surv(df.test$eventtime, df.test$status)

################################################################################
# Cox
################################################################################
cph <- cv.glmnet(short.scaled, short.y, family = "cox")

baseline <- predict(cph, short.scaled, type = "response", s="lambda.1se")
baseline <- sapply(keep.times, function(t) 1/sum(baseline[short.y[, 1] >= t]))

# Predict with test data:
preds <- predict(cph, short.scaled.test, type='response', s="lambda.1se")
stacked.coxl1.pred.matrix <- t(sapply(preds, function(pred) exp(-cumsum(baseline) * pred)))

# Bootstrap to get the CI around the IBS:
ci.coxl1    <- get.bootstrap.ci(short.y.test, stacked.coxl1.pred.matrix, keep.times)
# Point estimate of the IBS:
point.coxl1 <- get.one.ibs(short.y.test, stacked.coxl1.pred.matrix, keep.times)
# Save these results in a csv:
saveme(data.table(bs = ci.coxl1, point = point.coxl1), "Cox_Regularized")

################################################################################
# COX GBM
################################################################################
model  <- gbm(Surv(time, event)~., 
              data = data.frame(cbind(short.scaled, time = short.y[, 1], event = short.y[, 2])),
              n.trees = 5000,  
              shrinkage = .01,
              distribution="coxph",
              interaction.depth = 2);

# Compute the baseline hazard for the model (using train data)
baseline <- exp(predict(model, data.frame(short.scaled)))
baseline <- sapply(keep.times, function(t) 1/sum(baseline[short.y[, 1] >= t]))

# Predict with test data:
preds <- predict(model, data.frame(short.scaled.test))
preds <- c(sapply(preds, function(pred) exp(-cumsum(baseline) * exp(pred))))
pred.matrix <- matrix(preds, ncol=length(keep.times), nrow=length(all.ids), byrow = T)

# Measure performance
ci.gbm    <- get.bootstrap.ci(short.y.test, pred.matrix, keep.times)
point.gbm <- get.one.ibs(short.y.test, pred.matrix, keep.times)
saveme(data.table(bs = ci.gbm, point = point.gbm), "GBM")

################################################################################
# random survival forest
################################################################################
rfmodel <- rfsrc(Surv(time, event) ~ ., 
                 data = data.frame(cbind(short.scaled, time = short.y[, 1], event = short.y[, 2])), 
                 ntree = 500,
                 ntime = length(keep.times))
preds <- predict(rfmodel, data.frame(short.scaled.test))$survival

ci.rf    <- get.bootstrap.ci(short.y.test, preds, keep.times)
point.rf <- get.one.ibs(short.y.test, preds, keep.times)
saveme(data.table(bs = ci.rf, point = point.rf), "RF")

################################################################################
# stacked GBM
################################################################################
stacked.gbm <- gbm.fit(long.scaled,
                       long.y,
                       n.trees=5000,
                       interaction.depth = 2,
                       distribution="bernoulli")

preds   <- predict(stacked.gbm, long.scaled.test)
preds   <- 1/(1 + exp(-preds))

stacked.gbm.pred.matrix <- group.preds(preds, id.ix, all.ids, keep.times)
ci.gbm.st    <- get.bootstrap.ci(short.y.test, stacked.gbm.pred.matrix, keep.times)
point.gbm.st <- get.one.ibs(short.y.test, stacked.gbm.pred.matrix, keep.times)

saveme(data.table(bs = ci.gbm.st, point = point.gbm.st), "GBM_Stacked")

################################################################################
# stacked rf
################################################################################
stacked.rf <- ranger(y ~ ., 
                     data = cbind(long.scaled, y = long.y),
                     num.trees=100,
                     max.depth = 6
)

preds   <- predict(stacked.rf, long.scaled.test)$predictions

stacked.rf.pred.matrix <- group.preds(preds, id.ix, all.ids, keep.times)
ci.rf.st    <- get.bootstrap.ci(short.y.test, stacked.rf.pred.matrix, keep.times)
point.rf.st <- get.one.ibs(short.y.test, stacked.rf.pred.matrix, keep.times)

saveme(data.table(bs = ci.rf.st, point = point.rf.st), "RF_Stacked")

################################################################################
# stacked glinternet
################################################################################
stacked.gl <- glinternet(long.scaled, long.y,
                         numLevels = rep(1, 6),
                         interactionCandidates = 6,
                         family = "binomial"
                         )
preds   <- predict(stacked.gl, long.scaled.test, type="response")
preds   <- preds[, 50] # arbitrary choice!

stacked.glint.pred.matrix <- group.preds(preds, id.ix, all.ids, keep.times)
ci.glint    <- get.bootstrap.ci(short.y.test, stacked.glint.pred.matrix, keep.times)
point.glint <- get.one.ibs(short.y.test, stacked.glint.pred.matrix, keep.times)

saveme(data.table(bs = ci.glint, point = point.glint), "glinternet_Stacked")

################################################################################
# stacked nn
################################################################################
nn <- neuralnetwork(long.scaled, long.y,
                     hidden=c(10),
                     activ.functions = "sigmoid",
                     loss.type="log",
                     regression=F,
                     batch.size=1000,
                     optim.type = 'adam',
                     verbose=T,
                     drop.last=F,
                     n.epochs=500,
                     standardize=F)

preds   <- predict(nn, newdata = long.scaled.test)$probabilities[, 2]

stacked.preds <- group.preds(preds, id.ix, all.ids, keep.times)
ci.st.nn <- get.bootstrap.ci(short.y.test, stacked.preds, keep.times)
point.st.nn  <- get.one.ibs(short.y.test, stacked.preds, keep.times)

saveme(data.table(bs = ci.st.nn, point = point.st.nn), "neuralnet_Stacked")


################################################################################
# NNet-survival
# done in Python
################################################################################
# Write data frames for use in Python:
write.csv(df.train, paste0(SAVE_RESULTS_HERE, "train.csv"), row.names=F)
write.csv(df.test,  paste0(SAVE_RESULTS_HERE, "test.csv"),  row.names=F)

# # Packages
# import sys
# # Path to  https://github.com/MGensheimer/nnet-survival cloned to my machine
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
# trn = pd.read_csv("/Users/erincraig/Desktop/simulated/train.csv")
# tst = pd.read_csv("/Users/erincraig/Desktop/simulated/test.csv")
# 
# breaks = trn[trn.status == 1]["eventtime"].unique()
# breaks.sort()
# 
# n_intervals = len(breaks) - 1
# 
# y = nnet_survival.make_surv_array(trn["eventtime"],
#                                   trn["status"],
#                                   breaks)
# 
# x = trn[["X1", "X2", "X3", "X4", "X5"]].values
# ss = StandardScaler()
# xx = ss.fit_transform(X = x)
# 
# xtest = ss.transform(X = tst[["X1", "X2", "X3", "X4", "X5"]].values)
# 
# # Model
# model = Sequential()
# model.add(Dense(5, input_shape=(5,)))
# model.add(Dense(n_intervals, input_shape=(3,)))
# model.add(Activation('sigmoid'))
# model.compile(loss=nnet_survival.surv_likelihood(n_intervals), 
#               optimizer=optimizers.RMSprop(learning_rate = .001))
# early_stopping = EarlyStopping(monitor='loss', patience=2)
# 
# history=model.fit(x, y, batch_size=64, epochs=1000, verbose = 2, callbacks=[early_stopping])
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
# pd.DataFrame(predicted).to_csv("/Users/erincraig/Desktop/simulated/nnet_survival_predicted_for_R", index=False)

preds <- fread(paste0(SAVE_RESULTS_HERE, "nnet_survival_predicted_for_R"))
preds <- preds[2:nrow(preds), ]

ci.nnet    <- get.bootstrap.ci(short.y.test, as.matrix(preds), keep.times)
point.nnet <- get.one.ibs(short.y.test, as.matrix(preds), keep.times)

saveme(data.table(bs = ci.nnet, point = point.nnet), "nnet_survival")

################################################################################
# plot our results
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
                                         "Cox_Regularized",
                                         "GBM",
                                         "nnet_survival",
                                         "RF",
                                         "RF_Stacked",
                                         "neuralnet_Stacked",
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
       title = "Method comparison: simulated data",
       subtitle="(with 95% stratified bootstrap confidence intervals)",
       caption="(lower scores are preferred)") +
  scale_color_manual(breaks=c("original", "stacked"), values=c("#820000", "#F4795B")) 
ggsave(paste0(SAVE_RESULTS_HERE, "Figure2.pdf"), width=6, height=3)




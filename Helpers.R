library(riskRegression)
library(survival)
library(data.table)

# Function to convert hazard estimates to a matrix of survival probabilities.
# - preds is a vector of hazard predictions.
# - id.ix is a vector of length len(preds). id.ix[i] tells us which subject is associated with preds[i].
#    We assume that preds is in time order!
# - all.ids tells us the unique IDs.
# - keep.times are the event times at which we are computing the survival probabilities.
# This returns a matrix n subjects by m time points. 
# The i^th row gives us survival probabilities for subject all.ids[i], in order of keep.times.
group.preds <- function(preds, id.ix, all.ids, keep.times){
  grouped.preds <- matrix(nrow = length(all.ids), ncol = length(keep.times))
  i = 1
  for(id in all.ids){
    grouped.preds[i, ] <- preds[id.ix == id]
    i = i+1
  } 
  grouped.preds <- t(apply(1 - grouped.preds, 1, cumprod)) # survival at each time
  return(grouped.preds)
}


# Compute the integrated Brier score:
get.one.ibs <- function(surv_object, pred.matrix, keep.times){
  if(ncol(surv_object) == 3){
    df <- data.frame(time = surv_object[, 2], event = surv_object[, 3])
  } else {
    df <- data.frame(time = surv_object[, 1], event = surv_object[, 2])
  }
  
  score <-  Score(
    object = list(prob.matrix = 1-pred.matrix),
    formula = Surv(time, event) ~ 1,
    data = df,
    metrics = "Brier",
    summary = "ibs",
    times = keep.times,
    use.event.times = F
  )
  
  score$Brier$score[model == "prob.matrix" & 
                      times == max(score$Brier$score$times), IBS]
}


# Do the stratified bootstrap to allow the computation of a confidence interval around the IBS.
get.bootstrap.ci <- function(surv_object, pred.matrix, keep.times, n.bs = 1000){
  event <- surv_object[, ncol(surv_object)]
  
  n0 <- sum(event == 0)
  s0 <- which(event == 0)
  n1 <- sum(event == 1)
  s1 <- which(event == 1)
  
  event.times <- surv_object[, ncol(surv_object)-1]
  
  bs <- mclapply(1:n.bs, 
                 function(x) {
                   spl <- c(sample(s0, n0, replace=T), sample(s1, n1, replace=T)) 
                   keep.cols = unique(which(keep.times <= max(event.times)))
                   get.one.ibs(surv_object[spl, ], 
                               pred.matrix[spl, keep.cols], 
                               keep.times[keep.cols])
                 }
  )
  
  unlist(bs)
}
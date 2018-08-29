## data file: first column is normalized probe intensity, second column is probe sequence
## linux OS version

data.file = "./AB_a1_block3_combinatorial.txt"
seed.pwm.file = "./AB_a1_block3_8mers_pwm.txt"

source("code.R")
source("PBM_util.R")

L = 6
all.seqs.01.l = make.kmers.01(L)

## list of lmers that does not include reverse complement of itself
i=1
seen = new.env(hash=T)
for(i in 1:nrow(all.seqs.01.l)) {
  word.revcomp.seen = exists(decode01(rev(all.seqs.01.l[i,])), envir=seen)
  if ( !word.revcomp.seen ) assign(decode01(all.seqs.01.l[i,]),i,envir=seen)
}
x = ls(envir=seen)
good.lmer.idx = sapply(x, function(idx) get(idx,envir=seen))

rm(i)
rm(x)
rm(seen)
rm(word.revcomp.seen)

pbm.data = read.table(data.file, as.is=T)

## process sequence. this only needs to be done once for each array design, takes ~8 minutes on an imac (2006 model)
## probe sequences from Badis et. al. have 36 long variable region, I take the sequences of the first 40 bases into account
## which includes some of the costant region.
var.len = 40
num.lmers = var.len - L + 1

var.seqs = sapply(pbm.data[,2], function(seq) substr(seq, 1, var.len))
var.seqs.01 = t(sapply(var.seqs, encode01))

## for each lmer, figure out which probe contains it
## pbm.lmer.probes.idx is a list, each element of which is the vector of probe indices that contain a particular lmer, in either orientation
num.lmers = var.len - L + 1
pbm.lmer.probes.idx = list()
pbm.lmer.probes.idx[4^L+1]=NA
start.pos = seq(1, by=4, length.out=num.lmers)
for(probe.idx in 1:nrow(var.seqs.01)) {
  for (start in start.pos) {
    word = var.seqs.01[probe.idx, start:(start+(4*L-1))]
    word.idx = find.idx.01(word)
    word.revcomp.idx = find.idx.01(rev(word))
    pbm.lmer.probes.idx[[word.idx]] = c(pbm.lmer.probes.idx[[word.idx]], probe.idx)
    pbm.lmer.probes.idx[[word.revcomp.idx]] = c(pbm.lmer.probes.idx[[word.revcomp.idx]], probe.idx)
  }
}
pbm.lmer.probes.idx = lapply(pbm.lmer.probes.idx[1:(4^L)],unique)
## for each probe, assign lmer index
pbm.probes.lmer.idx.f = t(apply(var.seqs.01, 1, function(x) sapply(start.pos, function(i) find.idx.01(x[i:(i+L*4-1)]))))
pbm.probes.lmer.idx.r = t(apply(var.seqs.01, 1, function(x) {y = rev(x);sapply(rev(start.pos), function(i) find.idx.01(y[i:(i+L*4-1)]))}))


## figure out the position weights, the same sequence seem to have more signal farther away from the glass
pbm.lmer.medians = unlist(lapply(pbm.lmer.probes.idx, function(x) median(pbm.data[x,1],na.rm=T)))
## the signal intensities of top lmers are mostly (or only) affected by positional effect
num.top.lmers = 50
## order according to median signal intensity, and find the top 50 high affinity motifs
top.lmers = apply(all.seqs.01.l[order(pbm.lmer.medians,decreasing=T)[1:num.top.lmers],],1,decode01)
## for var.len = 40, weights of positions 1-length(pbm.position.weights) were estimated,
## length(pbm.position.weights)-40 is represented by 30
pbm.position.weights = est.position.weights(0, 0, var.seqs, pbm.data[,1], kmers = top.lmers)
pbm.position.weights = c(pbm.position.weights, rep(pbm.position.weights[length(pbm.position.weights)],
                                                   40 - length(pbm.position.weights)))

## without BEEML-PBM fitting at all, how much variance of probe intensities do lmer median intensities explain?
print("Correlation: measured probe intensity - top lmer median intensity of that probe")
nlmers.pred = apply(cbind(pbm.probes.lmer.idx.f, pbm.probes.lmer.idx.r), 1, function(x) max(pbm.lmer.medians[x]))
round(cor(nlmers.pred, pbm.data[, 1])^2, 2)
## lmer median intensities explain ~80% of variance of probe intensities? Position weighted
print("Correlation: measured probe intensity - position weighted top lmer median intensity of that probe")
pbm.lmers.pred = apply(cbind(pbm.probes.lmer.idx.f, pbm.probes.lmer.idx.r), 1,
                       function(x) max(pbm.lmer.medians[x] %*% c(pbm.position.weights[1:ncol(pbm.probes.lmer.idx.f)],
                                                                 pbm.position.weights[1:ncol(pbm.probes.lmer.idx.r)])))
pbm.lmers.pred.rsqr = round(cor(pbm.lmers.pred, pbm.data[,1])^2,2)
pbm.lmers.pred.rsqr

## take L contigous columns with highest information from uniprobe pwm to use as starting position for optimization
seed.mtx = get.high.information.positions(as.matrix(read.table(seed.pwm.file,skip=1)[,-1]), L, F)
## change the initial values
## seed.mtx[seed.mtx == sample(seed.mtx, 16)] = 0.1
seed.mtx
## BEEML modeling for mono-nucleotide PWM
pbm.sol.mono = get.beeml.solution.sum(values = pbm.data[,1]/sd(pbm.data[,1],na.rm=T), idx.f = pbm.probes.lmer.idx.f,
                                 idx.r = pbm.probes.lmer.idx.r, position.weights = pbm.position.weights[1:ncol(pbm.probes.lmer.idx.f)],
                                 seed.mtx = seed.mtx, seqs.01 = all.seqs.01.l, lmer.len = L, palidromic = F, nprint = 1, lambda = 0.1)

print(pbm.sol.mono)
# print.logo(pbm.sol.mono$mtx)
# print.logo(rev(pbm.sol.mono$mtx))

## how well does the beeml mono-nucleotide regression model perform?
print("how well does the beeml mono-nucleotide regression model perform?")
## fit at the probe level
print("fit at the probe level")
pbm.beeml.pred = predict.occupancy(cbind(pbm.sol.mono$mtx, pbm.sol.mono$mu), all.seqs.01.l, pbm.probes.lmer.idx.f, pbm.probes.lmer.idx.r, pbm.position.weights[1:ncol(pbm.probes.lmer.idx.f)])
pbm.beeml.pred.rsqr = round(cor(pbm.beeml.pred, pbm.data[,1])^2,2)
print(pbm.beeml.pred.rsqr)

## fit at the lmer median intensity level, excluding reverse complements
## BEEML predicted lmer median intensities
print("BEEML predicted lmer median intensities")
pbm.beeml.pred.lmers = unlist(lapply(pbm.lmer.probes.idx, function(x) median(pbm.beeml.pred[x],na.rm=T)))
## correlate with measured lmer median intensities
pbm.beeml.pred.lmers.rsqr = round(cor(pbm.beeml.pred.lmers[good.lmer.idx], pbm.lmer.medians[good.lmer.idx], use="complete")^2,2)
print(pbm.beeml.pred.lmers.rsqr)

## how well does the PWM from seed-and-wobble perform?
print("how well does the PWM from seed-and-wobble perform?")
pbm.sol.sNw = get.sNw.solution.sum(values = pbm.data[,1]/sd(pbm.data[,1],na.rm=T), idx.f = pbm.probes.lmer.idx.f,
                                      idx.r = pbm.probes.lmer.idx.r, position.weights = pbm.position.weights[1:ncol(pbm.probes.lmer.idx.f)],
                                      seed.mtx = seed.mtx, seqs.01 = all.seqs.01.l, palidromic = F, nprint = 1, lambda = 0.1)

print(pbm.sol.sNw)

## fit at probe level
print("fit at the probe level")
pbm.beeml.pred.sNw = predict.occupancy(cbind(seed.mtx, pbm.sol.sNw$mu), all.seqs.01.l, pbm.probes.lmer.idx.f, pbm.probes.lmer.idx.r, pbm.position.weights[1:ncol(pbm.probes.lmer.idx.f)])
pbm.beeml.pred.rsqr.sNw = round(cor(pbm.beeml.pred.sNw, pbm.data[,1])^2,2)
pbm.beeml.pred.rsqr.sNw
## fit at the lmer median intensity level, excluding reverse complements
print("fit at the median intensity level")
pbm.beeml.pred.lmers.sNw = unlist(lapply(pbm.lmer.probes.idx, function(x) median(pbm.beeml.pred.sNw[x],na.rm=T)))
pbm.beeml.pred.lmers.rsqr.sNw = round(cor(pbm.beeml.pred.lmers.sNw[good.lmer.idx], pbm.lmer.medians[good.lmer.idx], use="complete")^2,2)
pbm.beeml.pred.lmers.rsqr.sNw


###########################################################################################################
## The replicate array measurements

data.file2 = "./AB_a2_block4_combinatorial.txt"

pbm.data2 = read.table(data.file2, as.is=T)

var.seqs2 = sapply(pbm.data2[,2], function(seq) substr(seq, 1, var.len))
var.seqs2.01 = t(sapply(var.seqs2, encode01))

pbm.lmer.probes2.idx = list()
pbm.lmer.probes2.idx[4^L+1]=NA
start.pos = seq(1, by=4, length.out=num.lmers)
for(probe.idx in 1:nrow(var.seqs2.01)) {
  for (start in start.pos) {
    word = var.seqs2.01[probe.idx, start:(start+(4*L-1))]
    word.idx = find.idx.01(word)
    word.revcomp.idx = find.idx.01(rev(word))
    pbm.lmer.probes2.idx[[word.idx]] = c(pbm.lmer.probes2.idx[[word.idx]], probe.idx)
    pbm.lmer.probes2.idx[[word.revcomp.idx]] = c(pbm.lmer.probes2.idx[[word.revcomp.idx]], probe.idx)
  }
}
pbm.lmer.probes2.idx = lapply(pbm.lmer.probes2.idx[1:(4^L)],unique)
## for each probe, assign lmer index
pbm.probes2.lmer.idx.f = t(apply(var.seqs2.01, 1, function(x) sapply(start.pos, function(i) find.idx.01(x[i:(i+L*4-1)]))))
pbm.probes2.lmer.idx.r = t(apply(var.seqs2.01, 1, function(x) {y = rev(x);sapply(rev(start.pos), function(i) find.idx.01(y[i:(i+L*4-1)]))}))
## figure out the position weights, the same sequence seem to have more signal farther away from the glass
pbm.lmer.medians2 = unlist(lapply(pbm.lmer.probes2.idx, function(x) median(pbm.data2[x,1],na.rm=T)))

## order according to median signal intensity, and find the top 50 high affinity motifs
top.lmers2 = apply(all.seqs.01.l[order(pbm.lmer.medians2,decreasing=T)[1:num.top.lmers],],1,decode01)
## for var.len = 40, weights of positions 1-length(pbm.position.weights) were estimated,
## length(pbm.position.weights)-40 is represented by 30
pbm.position.weights2 = est.position.weights(0, 0, var.seqs2, pbm.data2[,1], kmers = top.lmers)
pbm.position.weights2 = c(pbm.position.weights2, rep(pbm.position.weights2[length(pbm.position.weights2)],
                                                     40 - length(pbm.position.weights2)))

## Experimental reproducibility: correlation lmer median intensities between two arrays
print("Experimental reproducibility")
expr.reproducibility.rsqr = round(cor(pbm.lmer.medians[good.lmer.idx], pbm.lmer.medians2[good.lmer.idx], use="complete")^2, 2)
print(expr.reproducibility.rsqr)

## BEEML mono-nucleotide model performance
## how well does the beeml mono-nucleotide regression model perform on the replicate array?
## fit at the probe level using array 1 model
pbm.beeml.pred2 = predict.occupancy(cbind(pbm.sol.mono$mtx, pbm.sol.mono$mu), all.seqs.01.l, pbm.probes2.lmer.idx.f, pbm.probes2.lmer.idx.r, pbm.position.weights2[1:ncol(pbm.probes2.lmer.idx.f)])
## fit at the lmer median intensity level, excluding reverse complements
## BEEML predicted array 2 lmer median intensities
pbm.beeml.pred.lmers2 = unlist(lapply(pbm.lmer.probes2.idx, function(x) median(pbm.beeml.pred2[x],na.rm=T)))
## correlate with measured array2 lmer median intensities
print("Correlation: predicted array 2 lmer median intensities-measured array2 lmer median intensities")
pbm.beeml.pred.lmers.rsqr2 = round(cor(pbm.beeml.pred.lmers2[good.lmer.idx], pbm.lmer.medians2[good.lmer.idx], use="complete")^2,2)
print(pbm.beeml.pred.lmers.rsqr2)


###########################################################################################################
## add di-nucleotide interaction terms to model array 1

## generate frequency matrix with the highest information content
mtx.IC <- get.lmer.pwm(as.matrix(read.table(file.path(working_dir,seed.pwm.file),skip=1)[,-1]), L, F)

## generate a list of all possible matrices with one di-nucleotide interaction term
## generate the sequence matrix
all.seqs.di = list()
## generate the initial values
seed.mtx.di = list()
## store the names
name.di = character()
for (i in 1:(L-1)) {
  for (j in (i+1):L) {
    di.coding <- t(apply(all.seqs.01.l[, c((4*i-3):(4*i), (4*j-3):(4*j))], 1, function(x) as.integer(x[1:4] %o% x[5:8])))
    out <- cbind(all.seqs.01.l, di.coding)
    name <- paste(as.factor(i), as.factor(j), sep=".")
    all.seqs.di[[name]] = out
    
    ## calculate initial dinucleotide interaction frequencies
    v <- as.vector(mtx.IC[,i] %o% mtx.IC[,j])
    ## calculate log-odds for each position&base, and normalized each position
    v.di <- sapply(-log(v), function(x) x-min(-log(v)))
    seed.mtx.di[[name]] = c(seed.mtx, v.di)
    
    name.di = c(name.di, name)
  }  
}

fitting.cor <- matrix(ncol=3, nrow=length(name.di), dimnames=list(name.di, c("Fit.probe", "Fit.median", "Fit.array2")))
beeml.di.model <- function(di.idx) {
  
  cycle = paste("Cycle", di.idx, sep=": ")
  print(cycle)
  
  ## nonlinear regression with di-nucleotide interaction
  ## values: normalized intensities by dividing standard deviation of the whole dataset
  pbm.sol = get.beeml.solution.sum(values = pbm.data[,1]/sd(pbm.data[,1],na.rm=T), idx.f = pbm.probes.lmer.idx.f,
                                   idx.r = pbm.probes.lmer.idx.r, position.weights = pbm.position.weights[1:ncol(pbm.probes.lmer.idx.f)],
                                   seed.mtx = seed.mtx.di[[di.idx]], seqs.01 = all.seqs.di[[di.idx]], lmer.len = L, palidromic = F,
                                   nprint = 1, lambda = 0.1)
  
  print(pbm.sol)
  
  ## how well does the beeml regression model perform?
  ## fit at the probe level
  print("Fit at the probe level:")
  pbm.beeml.pred.di = predict.occupancy(cbind(pbm.sol$mtx, pbm.sol$mu), all.seqs.di[[di.idx]], pbm.probes.lmer.idx.f, pbm.probes.lmer.idx.r, pbm.position.weights[1:ncol(pbm.probes.lmer.idx.f)])
  pbm.beeml.pred.rsqr.di = round(cor(pbm.beeml.pred.di, pbm.data[,1])^2,2)
  print(pbm.beeml.pred.rsqr.di)
    
  ## fit at the lmer median intensity level, excluding reverse complements
  print("Fit at the median intensity level:")
  ## BEEML predicted lmer median intensities
  pbm.beeml.pred.lmers.di = unlist(lapply(pbm.lmer.probes.idx, function(x) median(pbm.beeml.pred.di[x],na.rm=T)))
  ## correlate with measured lmer median intensities
  pbm.beeml.pred.lmers.di.rsqr = round(cor(pbm.beeml.pred.lmers.di[good.lmer.idx], pbm.lmer.medians[good.lmer.idx], use="complete")^2,2)
  print(pbm.beeml.pred.lmers.di.rsqr)
  
  ## BEEML di-nucleotide model performance
  ## how well does the beeml di-nucleotide regression model perform on the replicate array?
  ## fit at the probe level
  print("How well model 1 predicts array 2 median intensity:")
  pbm.beeml.pred3 = predict.occupancy(cbind(pbm.sol$mtx, pbm.sol$mu), all.seqs.di[[di.idx]], pbm.probes2.lmer.idx.f, pbm.probes2.lmer.idx.r, pbm.position.weights2[1:ncol(pbm.probes2.lmer.idx.f)])
  # ## fit at the lmer median intensity level, excluding reverse complements
  # ## BEEML predicted lmer median intensities
  pbm.beeml.pred.lmers3 = unlist(lapply(pbm.lmer.probes2.idx, function(x) median(pbm.beeml.pred3[x],na.rm=T)))
  # ## correlate with measured array2 lmer median intensities
  pbm.beeml.pred.lmers.rsqr3 = round(cor(pbm.beeml.pred.lmers3[good.lmer.idx], pbm.lmer.medians2[good.lmer.idx])^2,2)
  print(pbm.beeml.pred.lmers.rsqr3)
  
  c(pbm.beeml.pred.rsqr.di, pbm.beeml.pred.lmers.di.rsqr, pbm.beeml.pred.lmers.rsqr3)
}

for (di.idx in name.di) {
  coeff.fitting = beeml.di.model(di.idx)
  fitting.cor[di.idx,] = coeff.fitting
}
#beeml.di.model.pred = lapply(name.di, function(x) beeml.di.model(x))
outf = sub("_combinatorial.*", "", data.file)
write.table(fitting.cor, file.path(working_dir,paste0(outf,"_dinucleotide_fitting.txt")), sep="\t")

## print.logo(pbm.sol$mtx)


## make a plot of the results
# par(mfrow=c(2,2),pch=16,col=rgb(0,0,1,0.8),bty="n")
# hist(pbm.data[,1],breaks="Scott",main="Histogram of Probe Intensities")
# plot(pbm.lmers.pred, pbm.data[,1], main=paste("Probe Intensities from lmer medians\n R^2 =", pbm.lmers.pred.rsqr), xlab="lmer predicted probe intensities", ylab="Probe Intensities")
# plot(pbm.beeml.pred.lmers, pbm.lmer.medians, main=paste("BEEML Fit on lmer Median Intensities\n R^2 =", pbm.beeml.pred.lmers.rsqr), xlab="BEEML predicted lmer medians", ylab="PBM lmer Median Intensities")
# plot(pbm.beeml.pred, pbm.data[,1], main=paste("Probe Intensities from BEEML model\n R^2 =", pbm.beeml.pred.rsqr), xlab="BEEML predicted probe intensities", ylab="PBM Probe Intensities")

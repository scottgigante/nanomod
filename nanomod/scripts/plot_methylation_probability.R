library(ggplot2)
library(jsonlite)
library(reshape2)
library(gplots)
library(ggtern)
library(grid)
library("e1071")
library(gridExtra)
library(pROC)
library(randomForest)
library(MASS)
library(misc3d)
library(plotly)
library(webshot)

read_file <- function(in_filename, truth, subsample=1) {
  idx = list()
  if (FALSE) {
    mat = read.table(in_filename, header=TRUE, sep="\t", col.names=c("chromosome", "start", "end", "proportion", "modified", "unmodified", "total"))
    idx$start = 2
    idx$end = 3
    idx$prop= 4
    idx$mods = 5
    idx$unmods = 6
    idx$cov = 7
  } else {
    mat = read.table(in_filename, header=TRUE, sep="\t", col.names=c("chromosome", "pos", "proportion", "modified", "unmodified", "total"))
    idx$start = 2
    idx$end = 2
    idx$prop= 3
    idx$mods = 4
    idx$unmods = 5
    idx$cov = 6
  }
  df = data.frame(position=(mat[,idx$start]+mat[,idx$end])/2, methylation=mat[,idx$prop])
  df$coverage = mat[,idx$cov]
  df$modcounts = mat[,idx$mods]
  df$unmodcounts = mat[,idx$unmods]
  df$misscounts = df$coverage - df$modcounts - df$unmodcounts
  subsample_counts <- function(v, s) { if (s == 1) {
    v 
  } else {
    sapply(v, FUN = function(x) { sum( runif(x) < s ) } )
  } }
  df$modcounts = subsample_counts(df$modcounts, subsample)
  df$unmodcounts = subsample_counts(df$unmodcounts, subsample)
  df$misscounts = subsample_counts(df$misscounts, subsample)
  df$coverage = df$modcounts + df$unmodcounts + df$misscounts
  df$mods = df$modcounts/df$coverage
  df$unmods = df$unmodcounts/df$coverage
  df$misses = df$misscounts/df$coverage
  df = df[which(!is.na(df$misses)),]
  # Truth: 1 if methylated, -1 if not
  df$truth = truth
  df$ratio = (df$methylation/(1-df$methylation))
  df$ll = log(df$ratio)
  df$ll_corrected = df$ll * df$truth
  df$correctness = 1 - abs(df$methylation-(df$truth+1)/2)
  df$correct = 0
  df$correct[which(df$correctness>0.5)] = 1
  df$correct[which(df$correctness<0.5)] = -1
  df
}

correct_for_misses <- function(df, weight, threshold) {
  df$methylation = (df$modcounts + weight * df$misscounts + 1/2) / (df$modcounts + df$unmodcounts + 1)
  df$correct = ifelse((df$methylation > threshold & df$truth == 1) | (df$methylation < threshold & df$truth == -1), 1, -1)
  if (sum(df$methylation == threshold) > 0) df[which(df$methylation == threshold),]$correct = 0
  df
}

check_weighted_correctness <- function(df, weight, threshold) {
  df = correct_for_misses(df, weight, threshold)
  sum(df$correct == 1) / (sum(df$correct == 1) + sum(df$correct == -1))
}

plot_weighted_correctness <- function(df, weight, threshold) {
  df$corrected_methylation = (df$modcounts + weight * df$misscounts + 1/2) / (df$modcounts + df$unmodcounts + 1)
  phred=-10*log10(1-check_weighted_correctness(df, weight, threshold))
  p = ggplot(data=df, aes(x=corrected_methylation)) + 
    geom_density(aes(fill=factor(truth)), alpha=0.4) +
    geom_vline(xintercept = threshold) +
    labs(title=paste0("Weight = ", weight, ", Threshold = ", threshold),
         x="Methylation",
         y="Density",
         caption=paste0("Phred = ",phred),
         fill="Truth") +
    scale_x_continuous(limits=c(0,1)) + 
    scale_fill_manual(labels=c("Unmethyl","Methyl"), values=c("blue","red"))
  p
}

get_pass_prop <- function(df, weight=NA, threshold=NA, binwidth=0.02) {
  x = seq(from=0, to=1, by=binwidth)
  if (!is.na(weight) && !is.na(threshold)) df = correct_for_misses(df, weight, threshold)
  df_pass = df[which(df$correct==1),]
  df_fail = df[which(df$correct==-1),]
  pass_prop = sapply(x, FUN = function(x) { 
    npass = sum(df_pass$misses <= x + binwidth/2 & df_pass$misses > x - binwidth/2) 
    nfail = sum(df_fail$misses <= x + binwidth/2 & df_fail$misses > x - binwidth/2) 
    npass / (npass + nfail)
  })
}

num_sites = 693340

#canonical_filename = "/stornext/Home/data/allstaff/g/gigante.s/python/nanomod/data/lucattini/fasta/canonical_rerun.methylation.txt"
canonical_filename = "canonical.subsampled.methylation.txt"
canonical_truth = -1
#title = "dcm-/dam- Mutant"
#modified_filename = "/stornext/Home/data/allstaff/g/gigante.s/python/nanomod/data/lucattini/fasta/modified_rerun.methylation.txt"
modified_filename = "modified_rerun.methylation.txt"
modified_truth = 1
title = "Treated with M.SssI"
# methyl = read_json(in_filename)
# in_filename = modified_filename
# out_filename = sub(".txt$", "", in_filename)
# truth = modified_truth
# df = read_file(in_filename, truth)

df = rbind(read_file(modified_filename, modified_truth), read_file(canonical_filename, canonical_truth))
#out_filename = "/stornext/Home/data/allstaff/g/gigante.s/python/nanomod/data/lucattini/fasta/combined_rerun"
out_filename = "combined_rerun"

threshold_df = data.frame(thresholds=seq(from=0,to=5,by=0.01))
threshold_df$num_calls = sapply(threshold_df$thresholds, FUN = function(x) { sum(df$ll > x | df$ll < -x, na.rm=TRUE) } )
threshold_df$error_rate = sapply(threshold_df$thresholds, FUN = function(x) { 1 - sum(df$ll_corrected > x, na.rm=TRUE)/threshold_df$num_calls[which(threshold_df$thresholds == x)]  } ) * 100
threshold_df$prop_called = threshold_df$num_calls / num_sites # number of sites
threshold_df$prop_not_called = 1 - threshold_df$prop_called
threshold_df$correct_rate = 100 - threshold_df$error_rate


fit <- lm(misses ~ methylation, data = df)

df_pass = df[which(abs(df$methylation-(df$truth+1)/2)<0.5),]
df_fail = df[which(abs(df$methylation-(df$truth+1)/2)>0.5),]
cdf = data.frame(z=seq(from=0, to=1, by=0.001))
cdf$pass = sapply(cdf$z, FUN = function(x) { 1.0/nrow(df_pass) * sum(df_pass$misses <= x) } )
cdf$fail = sapply(cdf$z, FUN = function(x) { 1.0/nrow(df_fail) * sum(df_fail$misses <= x) } )

smooth=50
cdf$pass_pdf = sapply(1:nrow(cdf), FUN = function(x) { (cdf$pass[min(nrow(cdf),x+smooth)]-cdf$pass[max(1,x-smooth)])/(cdf$z[min(nrow(cdf),x+smooth)]-cdf$z[max(1,x-smooth)]) })
cdf$fail_pdf = sapply(1:nrow(cdf), FUN = function(x) { (cdf$fail[min(nrow(cdf),x+smooth)]-cdf$fail[max(1,x-smooth)])/(cdf$z[min(nrow(cdf),x+smooth)]-cdf$z[max(1,x-smooth)]) })
cdf$fail_pdf = cdf$fail_pdf * nrow(df_fail) / (nrow(df_pass) + nrow(df_fail))

# check time lag dependence 
df_forward = df[which(sapply(1:(nrow(df)-1), FUN = function(x) { df$position[x+1] == df$position[x] + 1 })),]
df_reverse = df[which(sapply(2:(nrow(df)), FUN = function(x) { df$position[x-1] == df$position[x] - 1 })),]

pdf(paste(out_filename, "definetti_density", "pdf", sep="."))
ggtern(data = df, aes(mods, unmods, misses)) +
  stat_density_tern(aes(fill=..level..), bins=60, geom='polygon') + 
  scale_colour_discrete(guide = FALSE)
dev.off()

ggtern() +
  stat_density_tern(data = df_fail, aes(mods, unmods, misses, alpha=..level..), fill="red", geom = "polygon", bins=15) +
  stat_density_tern(data = df_pass, aes(mods, unmods, misses, alpha=..level..), fill="blue", geom = "polygon", bins=15) +
  scale_alpha_continuous(range=c(0.1, 0.8)) + 
  scale_colour_discrete(guide = FALSE)

pdf(paste(out_filename, "definetti_density_fault", "pdf", sep="."))
ggtern(data=df_pass[sample(nrow(df_pass), 10000),], aes(mods, unmods, misses)) +
  geom_point(alpha=0.1) +
  stat_density_tern(aes(fill=..level..), alpha=0.1,geom = "polygon", bins=25)
dev.off()

df_rounded = data.frame(mods = round(df$mods, digits=1), unmods = round(df$unmods, digits=1), misses=round(df$misses, digits=1))
levs = sapply(0:10, FUN = function(x) { sapply(0:10, FUN = function(y) { sapply(0:10, FUN = function(z) { paste(x/10, y/10, z/10, sep="_" ) })}) })
pos = factor(paste(df_rounded$mods, df_rounded$unmods, df_rounded$misses, sep="_"), levels=levs)
bubble_df = data.frame(mods=rep(seq(from=0, to=1, by=0.1), each=121), unmods=rep(seq(from=0, to=1, by=0.1), each=11), misses=seq(from=0, to=1, by=0.1))
bubble_df_indices = lapply(levs, FUN = function(x) { which(pos == x) })
bubble_df$count = sapply(bubble_df_indices, FUN = length)
bubble_df$log_count = log(bubble_df$count)
bubble_df$coverage = sapply(bubble_df_indices, FUN = function(x) { mean(df$coverage[x]) })
bubble_df = bubble_df[which(bubble_df$count != 0),]

pdf(paste(out_filename, "definetti_binned", "pdf", sep="."))
ggtern(data=bubble_df, aes(mods, unmods, misses)) + 
  geom_point(aes(size=count, colour=count), alpha=0.8) +
  labs(title=paste0("Density of count data normalised by coverage and binned to nearest 0.1 (n=",nrow(df),")"),
       x="n_M",
       y="n_C",
       z="n_X",
       size="Number of sites",
       colour="Number of sites") + 
  guides(colour = guide_legend())
dev.off()

# check choices of weighting of misses and call thresholds
weights_range = seq(from=0.035, to=0.055, length.out=11)
thresholds_range = seq(from=0.5, to=0.501, length.out=11)
thresholds_range = c(0.5, 0.5)
weight_threshold_matrix = sapply(weights_range, FUN = function(w) {
  sapply(thresholds_range, FUN = function(t) {
    check_weighted_correctness(df, w, t)
  })
})
rownames(weight_threshold_matrix) = thresholds_range
colnames(weight_threshold_matrix) = weights_range
heatmap.2(weight_threshold_matrix, Colv=FALSE, Rowv=FALSE, dendrogram="none", ylab="call threshold", xlab="misses weight")
p1 = plot_weighted_correctness(df, 0, 0.5)
# threshold only: 0.431
p2 = plot_weighted_correctness(df, 0, 0.431)
# weight only: 0.039
p3 = plot_weighted_correctness(df, 0.039, 0.5)
# ideal: weight 0.013, threshold 0.445
p4 = plot_weighted_correctness(df, 0.013, 0.445)
pdf(paste(out_filename, "weight_threshold_adjustment", "pdf", sep="."))
grid.arrange(p1, p2, p3, p4, ncol=2)
dev.off()

# check dependence of correctness on misses, stratified by coverage
binwidth = 0.02
plot_df = data.frame(misses=seq(from=0, to=1, by=binwidth))
plot_df$`1-20` = get_pass_prop(df[which(df$coverage - df$misscounts > 0 & df$coverage - df$misscounts <= 20),], 0, 0.5)
plot_df$`21-26` = get_pass_prop(df[which(df$coverage - df$misscounts > 20 & df$coverage - df$misscounts <= 26),], 0, 0.5)
plot_df$`27-32` = get_pass_prop(df[which(df$coverage - df$misscounts > 26 & df$coverage - df$misscounts <= 32),], 0, 0.5)
plot_df$`33-38` = get_pass_prop(df[which(df$coverage - df$misscounts > 32 & df$coverage - df$misscounts <= 38),], 0, 0.5)
plot_df$`39-47` = get_pass_prop(df[which(df$coverage - df$misscounts > 38 & df$coverage - df$misscounts <= 47),], 0, 0.5)
plot_df$`48+` = get_pass_prop(df[which(df$coverage - df$misscounts > 47),], 0, 0.5)
plot_df$all = get_pass_prop(df, 0, 0.5)
plot_df = melt(plot_df, id='misses')
plot_df = plot_df[which(!is.nan(plot_df$value)),]
colnames(plot_df) = c("misses", "coverage", "pass_prop")
pdf(paste(out_filename, "pass_prop_by_misses", "pdf", sep="."))
ggplot() +
  geom_line(data=plot_df, aes(x=misses, y=pass_prop, colour=coverage)) +  
  labs(y="Proportion of sites called correctly",
       x="Proportion of misses")
dev.off()

###########################################################
# Pass by proportion of X stratified by truth
###########################################################

plot_pass_by_prop_X <- function(df, corrected=FALSE, stratify=TRUE) {
  binwidth = 0.02
  plot_df = data.frame(misses=seq(from=0, to=1, by=binwidth))
  plot_df$correctness_all = get_pass_prop(df, 0, 0.5)
  if (corrected) plot_df$correctness_corrected = sapply(plot_df$misses, FUN = function(x) { s1 = x - binwidth/2; s2 = x+binwidth/2; sum(with(df, correct_corrected == 1 & misses >= s1 & misses < s2))/sum(with(df, misses >= s1 & misses < s2)) })
  if (stratify) plot_df$correctness_M = sapply(plot_df$misses, FUN = function(x) { 
    s1 = x - binwidth/2
    s2 = x+binwidth/2
    sum(with(df, 
             truth == 1 & 
               correct == 1 & 
               misses >= s1 & 
               misses < s2
    ))/sum(with(df, 
                truth == 1 & 
                  misses >= s1 & 
                  misses < s2)) })
  if (corrected && stratify) plot_df$correctness_corrected_M = sapply(plot_df$misses, FUN = function(x) { 
    s1 = x - binwidth/2
    s2 = x+binwidth/2
    sum(with(df, 
             truth == 1 & 
               correct_corrected == 1 & 
               misses >= s1 & 
               misses < s2
    ))/sum(with(df, 
                truth == 1 & 
                  misses >= s1 & 
                  misses < s2)) })
  if (stratify) plot_df$correctness_C = sapply(plot_df$misses, FUN = function(x) { s1 = x - binwidth/2; s2 = x+binwidth/2; sum(with(df, truth == -1 & correct == 1 & misses >= s1 & misses < s2))/sum(with(df, truth == -1 & misses >= s1 & misses < s2)) })
  if (corrected && stratify) plot_df$correctness_corrected_C = sapply(plot_df$misses, FUN = function(x) { s1 = x - binwidth/2; s2 = x+binwidth/2; sum(with(df, truth == -1 & correct_corrected == 1 & misses >= s1 & misses < s2))/sum(with(df, truth == -1 & misses >= s1 & misses < s2)) })
  plot_df = melt(plot_df, id="misses")
  p = ggplot() +
    geom_line(data=plot_df, aes(x=misses, y=value, colour=variable)) +  
    labs(y="Proportion of sites called correctly",
         x="Proportion of X's")
  p
}
pdf(paste(out_filename, "pass_prop_by_misses_by_truth", "pdf", sep="."))
print(plot_pass_by_prop_X(df))
dev.off()

###########################################################
# Dynamic calculation of threshold by misses
###########################################################
binwidth = 0.02
plot_df = data.frame(misses=seq(from=0, to=1, by=binwidth))
threshold=0.5
plot_df$original = sapply(plot_df$misses, FUN = function(loc) { 
  s1 = loc - binwidth/2; s2 = loc +binwidth/2
  sum(with(df, ((truth == -1 & methylation < threshold) | (truth == 1 & methylation > threshold)) & misses >= s1 & misses < s2))/sum(with(df, misses >= s1 & misses < s2)) 
})
plot_df$thresholds = unlist(sapply(plot_df$misses, FUN = function(x) { s1 = x - binwidth/2; s2 = x+binwidth/2; df_subsection = df[which(with(df,misses >= s1 & misses < s2)),]; coords(roc(df_subsection$truth, df_subsection$methylation), "best", ret="threshold")[1] }))
plot_df$dynamic = apply(plot_df, 1, FUN = function(x) { 
  loc = x[1]; threshold = x[3]; s1 = loc - binwidth/2; s2 = loc +binwidth/2
  sum(with(df, ((truth == -1 & methylation < threshold) | (truth == 1 & methylation > threshold)) & misses >= s1 & misses < s2))/sum(with(df, misses >= s1 & misses < s2)) 
})
plot_df$counts = sapply(plot_df$misses, FUN = function(x) { s1 = x - binwidth/2; s2 = x +binwidth/2; sum(with(df, misses >=s1 & misses < s2)) })
ggplot(data=plot_df, aes(x=misses)) + 
  geom_line(aes(y=original), colour="blue") + 
  geom_line(aes(y=dynamic), colour="red")

############################################################
# Basic pass by proportion of X
############################################################
binwidth = 0.02
plot_df = data.frame(misses=seq(from=0, to=1, by=binwidth))
plot_df$`W = 0, t = 0.5` = get_pass_prop(df, 0, 0.5)
pdf(paste(out_filename, "pass_prop_by_misses", "pdf", sep="."))
ggplot() +
  geom_line(data=plot_df, aes(x=misses, y=`W = 0, t = 0.5`)) +  
  labs(y="Proportion of sites called correctly",
       x="Proportion of X's")
dev.off()

# check dependence of correctness on misses with different thresholds and weights
plot_df$`W = 0.039, t = 0.5` = get_pass_prop(df, 0.039, 0.5)
plot_df$`W = 0, t = 0.431` = get_pass_prop(df, 0, 0.431)
plot_df$`W = 0.013, t = 0.445` = get_pass_prop(df, 0.013, 0.445)
plot_df = melt(plot_df, id="misses")
#plot_df$value = -10*log10(1-plot_df$value)
plot_df = plot_df[which(!is.nan(plot_df$value)),]
pdf(paste(out_filename, "pass_prop_by_misses_corrected", "pdf", sep="."))
ggplot() +
  geom_line(data=plot_df, aes(x=misses, y=value, colour=variable)) +  
  labs(y="Proportion of sites called correctly",
       x="Proportion of X's")
dev.off()

# get confidence scores by our various parameters of interest - does what does likelihood of correctness depend on?
correctness = array(0, dim=c(max(df$modcounts), max(df$unmodcounts), max(df$misscounts), 2))
apply(df, 1, FUN = function(x) {
  modcounts = x[4]
  unmodcounts = x[5]
  misscounts = x[6]
  correct = x[15]
  if (correct == 1) {
    correctness[modcounts, unmodcounts, misscounts, 1] <<- correctness[modcounts, unmodcounts, misscounts, 1] + 1
  } else if (correct == -1) {
    correctness[modcounts, unmodcounts, misscounts, 2] <<- correctness[modcounts, unmodcounts, misscounts, 2] + 1
  }
  NULL
})
correctness_df = melt(correctness)
colnames(correctness_df) = c("modcounts","unmodcounts","misscounts","correct","count")
correctness_df = correctness_df[which(correctness_df$correct==1),c("modcounts","unmodcounts","misscounts")]
correctness_df$phred = apply(correctness_df, 1, FUN = function(x) { 
  prop = correctness[x[1],x[2],x[3],1] / ( correctness[x[1],x[2],x[3],1] + correctness[x[1],x[2],x[3],2] )
  -10 * log10(1 - prop)
})
correctness_df[which(correctness_df$phred == Inf),"phred"] = max(correctness_df[which(!is.nan(correctness_df$phred) & correctness_df$phred != Inf),]$phred)
correctness_df = correctness_df[which(!is.nan(correctness_df$phred)),]
correctness_df$coverage = correctness_df$modcounts + correctness_df$unmodcounts + correctness_df$misscounts
correctness_df$mods = correctness_df$modcounts / correctness_df$coverage
correctness_df$unmods = correctness_df$unmodcounts / correctness_df$coverage
p1 = ggplot(data=correctness_df, aes(x=mods, y=mods, colour=phred)) + geom_point(alpha=0.1) + scale_colour_continuous(guide=FALSE)
p2 = ggplot(data=correctness_df, aes(x=mods, y=unmods, colour=phred)) + geom_point(alpha=0.1) + scale_colour_continuous(guide=FALSE)
p3 = ggplot(data=correctness_df, aes(x=mods, y=coverage, colour=phred)) + geom_point(alpha=0.1) + scale_colour_continuous(guide=FALSE)
p4 = ggplot(data=correctness_df, aes(x=unmods, y=mods, colour=phred)) + geom_point(alpha=0.1) + scale_colour_continuous(guide=FALSE)
p5 = ggplot(data=correctness_df, aes(x=unmods, y=unmods, colour=phred)) + geom_point(alpha=0.1) + scale_colour_continuous(guide=FALSE)
p6 = ggplot(data=correctness_df, aes(x=unmods, y=coverage, colour=phred)) + geom_point(alpha=0.1) + scale_colour_continuous(guide=FALSE)
p7 = ggplot(data=correctness_df, aes(x=coverage, y=mods, colour=phred)) + geom_point(alpha=0.1) + scale_colour_continuous(guide=FALSE)
p8 = ggplot(data=correctness_df, aes(x=coverage, y=unmods, colour=phred)) + geom_point(alpha=0.1) + scale_colour_continuous(guide=FALSE)
p9 = ggplot(data=correctness_df, aes(x=coverage, y=coverage, colour=phred)) + geom_point(alpha=0.1) + scale_colour_continuous(guide=FALSE)
grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,ncol=3)

# bin correctness over 5 counts
correctness_df = as.data.frame(do.call("rbind",lapply(seq(from=2.5, to=82.5, by=5), FUN = function(i) {
  do.call("rbind",lapply(seq(from=2.5, to=172.5, by=5), FUN = function(j) {
    t(as.data.frame(lapply(seq(from=2.5, to=157.5, by=5), FUN = function(k) {
      c(i,j,k,sum(correctness[(i-1.5):(i+2.5),(j-1.5):(j+1.5),(k-1.5):(k+2.5),1])/sum(correctness[(i-1.5):(i+2.5),(j-1.5):(j+2.5),(k-1.5):(k+2.5),]))
    })))
  }))
})))
colnames(correctness_df) = c("modcounts","unmodcounts","misscounts","correct")
correctness_df = correctness_df[which(!is.nan(correctness_df$correct)),]
plot_ly(correctness_df, x=~modcounts, y=~unmodcounts, z=~misscounts, type="scatter3d", mode="markers", color=~correct, marker=list(size=4, opacity=0.3))

###########################################
# build methylation by linear fit
###########################################
train_idx = sample(nrow(df), 100000)
train_df = df[train_idx,]
test_df = df[-train_idx,]
fit = lm(truth~modcounts + unmodcounts + misscounts, data=train_df)
train_df$methylation = (predict(fit, train_df) + 1)/2
roc_curve = roc(truth~methylation, data=train_df)
threshold = coords(roc_curve, "best", ret="threshold")
test_df$methylation = (predict(fit, test_df) + 1)/2
test_df$correct = ifelse((test_df$methylation > threshold & test_df$truth == 1) |
                           (test_df$methylation < threshold & test_df$truth == -1), 1, -1)
phred=-10*log10(sum(test_df$correct==-1)/nrow(test_df))
coefficients = signif(fit$coefficients, digits=2)
pdf(paste(out_filename, "linear_fit", "pdf", sep="."))
ggplot(data=test_df, aes(x=methylation)) + 
  geom_density(aes(fill=factor(truth)), alpha=0.4) +
  geom_vline(xintercept = threshold) +
  labs(title=paste0("P = ",coefficients[1]," + ",coefficients[2]," * n_M + ",coefficients[3]," * n_C + ", coefficients[4], " * n_X, Threshold = ", signif(threshold, digits=2)),
       x="Methylation",
       y="Density",
       caption=paste0("Phred = ",phred),
       fill="Truth") +
  #scale_x_continuous(limits=c(0,1)) + 
  scale_fill_manual(labels=c("Unmethyl","Methyl"), values=c("blue","red"))
dev.off()
linear_df = test_df

###########################################
# build methylation on polynomial fit
###########################################
#train_idx = sample(nrow(df), 200000)
train_df = df[train_idx,]
test_df = df[-train_idx,]
fit = nls(truth~a*modcounts*modcounts + b*unmodcounts*unmodcounts + c*misscounts*misscounts + d*modcounts*unmodcounts + e*modcounts*misscounts + f*unmodcounts*misscounts + g*modcounts + h * unmodcounts + i*misscounts + j, data=train_df, start=list(
  a=runif(1,min=-0.1, max=0.1), 
  b=runif(1,min=-0.1, max=0.1), 
  c=runif(1,min=-0.1, max=0.1), 
  d=runif(1,min=-0.1, max=0.1), 
  e=runif(1,min=-0.1, max=0.1), 
  f=runif(1,min=-0.1, max=0.1), 
  g=runif(1,min=-0.1, max=0.1), 
  h=runif(1,min=-0.1, max=0.1), 
  i=runif(1,min=-0.1, max=0.1), 
  j=runif(1,min=-0.1, max=0.1)))
train_df$methylation = (predict(fit, train_df) + 1)/2
roc_curve = roc(truth~methylation, data=train_df)
threshold = coords(roc_curve, "best", ret="threshold")
test_df$methylation = (predict(fit, test_df) + 1)/2
test_df$correct = ifelse((test_df$methylation > threshold & test_df$truth == 1) |
                           (test_df$methylation < threshold & test_df$truth == -1), 1, -1)
phred=-10*log10(sum(test_df$correct==-1)/nrow(test_df))
coefficients = signif(fit$m$getPars(),digits=2)
pdf(paste(out_filename, "nonlinear_fit", "pdf", sep="."))
ggplot(data=test_df, aes(x=methylation)) + 
  geom_density(aes(fill=factor(truth)), alpha=0.4) +
  geom_vline(xintercept = threshold) +
  labs(title=paste0("P = ",coefficients[10]," + ",coefficients[7]," * n_M + ",coefficients[8]," * n_C + ", coefficients[9], " * n_X +\n",
                    coefficients[4], " * n_M * n_C + ", coefficients[5], " * n_M * n_X  + ", coefficients[6], " * n_C * n_X +\n",
                    coefficients[1], " * n_M^2 + ", coefficients[2], " * n_C^2 + ", coefficients[3], " * n_X^2, Threshold = ", signif(threshold, digits=2)),
       x="Methylation",
       y="Density",
       caption=paste0("Phred = ",phred),
       fill="Truth") +
  #scale_x_continuous(limits=c(0,1)) + 
  scale_fill_manual(labels=c("Unmethyl","Methyl"), values=c("blue","red"))
dev.off()
polynomial_df = test_df

###########################################
# build methylation with SVM
###########################################
df$truth = factor(df$truth)
#train_idx = sample(nrow(df), 100000)
train_df = df[train_idx,]
test_df = df[-train_idx,]
#tuned_svm = tune(svm, truth~modcounts + unmodcounts + misscounts, data=train_df, ranges=list(cost=2^(-2:6), gamma=2^(-1:1)))
# cost=8, gamma=2
svm_fit = svm(truth~modcounts + unmodcounts + misscounts, data=train_df, probability=TRUE)
prediction = predict(svm_fit, test_df, probability=TRUE)
test_df$correct = ifelse(prediction == test_df$truth, 1, -1)
#test_df$phred=sapply(1:nrow(test_df), FUN = function(x) {
#  truth = test_df$truth[x]
#  probability = 1 - attr(prediction, "probabilities")[x, truth]
#  -10*log10(probability)
#})
test_df$methylation = attr(prediction, "probabilities")[,1]
svm_phred = -10*log10(sum(test_df$correct==-1)/nrow(test_df))
threshold=0.5
pdf(paste(out_filename, "svm_fit", "pdf", sep="."))
ggplot(data=test_df, aes(x=methylation)) + 
  geom_histogram(position = "dodge", aes(fill=factor(truth)), alpha=0.4) +
  geom_vline(xintercept = threshold) +
  labs(title="Support Vector Machine",
       x="Methylation",
       y="Density",
       caption=paste0("Phred = ",svm_phred),
       fill="Truth") +
  scale_fill_manual(labels=c("Unmethyl","Methyl"), values=c("blue","red")) + 
  scale_y_sqrt()
dev.off()
svm_df = test_df

###########################################
# build methylation with RF
###########################################
df$truth = factor(df$truth)
#train_idx = sample(nrow(df), 100000)
train_df = df[train_idx,]
test_df = df[-train_idx,]
fit = randomForest(truth~modcounts + unmodcounts + misscounts, data=train_df, importance=TRUE, ntree=1000)
probabilities = predict(fit, test_df, type="prob")
prediction = ifelse(predict(fit, test_df, type="response") == 1, 1, -1)
test_df$methylation = probabilities[,2]
test_df$truth = as.numeric(as.character(test_df$truth))
test_df$correct = ifelse(as.character(prediction) == as.character(test_df$truth), 1, -1)
#test_df$phred=sapply(1:nrow(test_df), FUN = function(x) {
#  truth = ifelse(test_df$truth[x]==1, 2, 1)
#  probability = 1 - probabilities[x, truth]
#  -10*log10(probability)
#})
phred = -10*log10(sum(test_df$correct==-1)/nrow(test_df))
threshold=0.5
pdf(paste(out_filename, "rf_fit", "pdf", sep="."))
ggplot(data=test_df, aes(x=methylation)) + 
  geom_histogram(position = "dodge", aes(fill=factor(truth)), alpha=0.4) +
  geom_vline(xintercept = threshold) +
  labs(title="Random Forest",
       x="Methylation",
       y="Density",
       caption=paste0("Phred = ",phred),
       fill="Truth") +
  scale_fill_manual(labels=c("Unmethyl","Methyl"), values=c("blue","red")) + 
  scale_y_sqrt()
dev.off()
rf_df = test_df

#####################################################
# check if above fits improve dependence between correctness and misses
#####################################################

# check dependence of correctness on misses
binwidth = 0.02
plot_df = data.frame(misses=seq(from=0, to=1, by=binwidth))
plot_df$naive = get_pass_prop(df, 0, 0.5)
plot_df$linear = get_pass_prop(linear_df, NA, NA)
plot_df$polynomial = get_pass_prop(polynomial_df, NA, NA)
plot_df$svm = get_pass_prop(svm_df, NA, NA)
plot_df$`random forest` = get_pass_prop(rf_df, NA, NA)
plot_df = melt(plot_df, id="misses")
pdf(paste(out_filename, "pass_prop_by_misses_by_method", "pdf", sep="."))
ggplot() +
  geom_line(data=plot_df, aes(x=misses, y=value, colour=variable)) +  
  labs(y="Proportion of sites called correctly",
       x="Proportion of X's",
       colour="Method") +
  coord_cartesian(ylim=c(0.76, 1))
dev.off()

###################################################################
# The T. Speed Attack - Bayesian Approximation
###################################################################

select = function(df, expr) {
  df[which(with(df, eval(parse(text=expr)))),]
}

bins=256
smooth=4
pr_p_m_given_M = density(select(df, "truth==1")$mods, n=bins, bw=smooth/bins, from=0, to=1)
pr_p_m_given_C = density(select(df, "truth==-1")$mods, n=bins, bw=smooth/bins, from=0, to=1)
stopifnot(pr_p_m_given_M$x == pr_p_m_given_C$x)
log_pr_p_m_given_truth = data.frame(idx=pr_p_m_given_M$x, empirical=log(pr_p_m_given_M$y / pr_p_m_given_C$y))
log_pr_p_m_given_truth[which(with(log_pr_p_m_given_truth, empirical == Inf | empirical == -Inf)),] = NA
log_pr_p_m_given_truth_fit = loess(empirical~idx, data=log_pr_p_m_given_truth, control=loess.control(surface="direct"))
log_pr_p_m_given_truth$predicted = predict(log_pr_p_m_given_truth_fit, log_pr_p_m_given_truth$idx)


# plots
#####################################################
p <- plot_ly(x=~pr_p_m_given_M$x, y=~pr_p_m_given_M$y) %>%
  layout(xaxis = list(title="Pr(p_m)"), yaxis=list(title="Density"), title="Pr(p_m | M)")
export(p, file=paste(out_filename, "pr_p_m_given_M", "png", sep="."))
p <- plot_ly(x=~pr_p_m_given_C$x, y=~pr_p_m_given_C$y) %>%
  layout(xaxis = list(title="Pr(p_m)"), yaxis=list(title="Density"), title="Pr(p_m | C)")
export(p, file=paste(out_filename, "pr_p_m_given_C", "png", sep="."))
p <- plot_ly(log_pr_p_m_given_truth, x=~idx, y=~empirical) %>%
  layout(xaxis = list(title="Pr(p_m)"), yaxis=list(title="Log Ratio"), title="Log Pr(p_m | M) / Pr(p_m | C)")
export(p, file=paste(out_filename, "log_pr_p_m_given_truth", "png", sep="."))
df$ll_corrected = predict(log_pr_p_m_given_truth_fit, df$mods)
roc_curve = roc(truth~ll_corrected, data=df)
threshold = coords(roc_curve, "best", ret="threshold")[1]
df$call_corrected = ifelse(df$ll_corrected > threshold, 1, -1)
df$correct_corrected = ifelse(df$call == df$truth, 1, -1)
pdf(paste(out_filename, "pass_prop_by_misses_by_truth_bayes_1d", "pdf", sep="."))
print(plot_pass_by_prop_X(df, TRUE))
dev.off()
#####################################################

bins=128
smooth = 4
illegal_idx = sapply(1:(bins**2), FUN = function(x) { i = (x-1) %/% bins + 1; j = (x-1) %% bins + 1; ifelse(i+j > bins + 1, TRUE, FALSE) })

pr_p_m_p_x_given_M =
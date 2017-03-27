library(ggplot2)
library(jsonlite)
library(reshape2)
library(gplots)
library(ggtern)
library(grid)
library("e1071")
library(gridExtra)
library(pROC)
library(MASS)
library(misc3d)
library(plotly)
library(glmnet)

read_file <- function(in_filename, truth, subsample=1) {
  idx = list()
  if (FALSE) {
    mat = read.table(in_filename, header=TRUE, sep="\t", col.names=c("chromosome", "start", "end", "proportion", "M", "U", "u", "X", "x", "D", "d"))
    idx$start = 2
    idx$end = 3
    idx$prop= 4
    idx$M = 5
    idx$U = 6
    idx$u = 7
    idx$X = 8
    idx$x = 9
    idx$D = 10
    idx$d = 11
  } else {
    mat = read.table(in_filename, header=TRUE, sep="\t", col.names=c("chromosome", "pos", "proportion", "M", "U", "u", "X", "x", "D", "d"))
    idx$start = 2
    idx$end = 2
    idx$prop= 3
    idx$M = 4
    idx$U = 5
    idx$u = 6
    idx$X = 7
    idx$x = 8
    idx$D = 9
    idx$d = 10
  }
  df = data.frame(
    position=(mat[,idx$start]+mat[,idx$end])/2, 
    methylation=mat[,idx$prop],
    M = mat[,idx$M],
    U = mat[,idx$U],
    u = mat[,idx$u],
    X = mat[,idx$X],
    x = mat[,idx$x],
    D = mat[,idx$D],
    d = mat[,idx$d]
  )
  df$coverage = df$M + df$U + df$X + df$D + df$u + df$x + df$d
  df$pM = df$M/df$coverage
  df$pU = df$U/df$coverage
  df$pu = df$u/df$coverage
  df$pX = df$X/df$coverage
  df$px = df$x/df$coverage
  df$pD = df$D/df$coverage
  df$pd = df$d/df$coverage
  df = df[which(df$coverage != 0),]
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
  df$misses = df$pX + df$px + df$pD + df$pd
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
canonical_filename = "detailed/canonical.methylation.txt"
canonical_truth = -1
#title = "dcm-/dam- Mutant"
#modified_filename = "/stornext/Home/data/allstaff/g/gigante.s/python/nanomod/data/lucattini/fasta/modified_rerun.methylation.txt"
modified_filename = "detailed/modified.methylation.txt"
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

# look at the data
totals_meth = melt(colSums(df[which(df$truth==1),c("M","U","u","X","x","D","d")]))
totals_meth$base = rownames(totals_meth)
totals_meth$truth = 1
totals_unmeth = melt(colSums(df[which(df$truth==-1),c("M","U","u","X","x","D","d")]))
totals_unmeth$base = rownames(totals_unmeth)
totals_unmeth$truth = -1
totals = rbind(totals_meth, totals_unmeth)
pdf(paste(out_filename, "6D_totals", "pdf", sep="."))
ggplot(data=totals, aes(x=base, y=value, fill=
                          factor(truth, levels=c(-1, 1), labels=c("Unmethyl","Methyl")))) +
  geom_bar(stat="identity", position="dodge") + 
  labs(fill="Truth", y="Count", x="Base")
dev.off()

density_df = melt(df, id="truth", measure.vars = c("pM","pU","pu","pX","px","pD","pd"))
pdf(paste(out_filename, "6D_densities", "pdf", sep="."))
ggplot(data=density_df, aes(x=value)) +
  facet_grid(factor(truth, levels=c(-1, 1), labels=c("Unmethyl", "Methyl"))~variable, scales="free") +
  geom_density(fill="grey50", alpha=0.4) + 
  coord_cartesian(ylim=c(0,2))
dev.off()

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

binwidth = 0.02
plot_df = data.frame(misses=seq(from=0, to=1, by=binwidth))
plot_df$correctness_all = get_pass_prop(df, 0, 0.5)
plot_df$correctness_corrected = sapply(plot_df$misses, FUN = function(x) { s1 = x - binwidth/2; s2 = x+binwidth/2; sum(with(df, correct_corrected == 1 & misses >= s1 & misses < s2))/sum(with(df, misses >= s1 & misses < s2)) })
plot_df$correctness_M = sapply(plot_df$misses, FUN = function(x) { 
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
plot_df$correctness_corrected_M = sapply(plot_df$misses, FUN = function(x) { 
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
plot_df$correctness_C = sapply(plot_df$misses, FUN = function(x) { s1 = x - binwidth/2; s2 = x+binwidth/2; sum(with(df, truth == -1 & correct == 1 & misses >= s1 & misses < s2))/sum(with(df, truth == -1 & misses >= s1 & misses < s2)) })
plot_df$correctness_corrected_C = sapply(plot_df$misses, FUN = function(x) { s1 = x - binwidth/2; s2 = x+binwidth/2; sum(with(df, truth == -1 & correct_corrected == 1 & misses >= s1 & misses < s2))/sum(with(df, truth == -1 & misses >= s1 & misses < s2)) })
#plot_df$call_methylated = sapply(plot_df$misses, FUN = function(x) { s1 = x - binwidth/2; s2 = x+binwidth/2; sum(with(df, call == 1 & misses >= s1 & misses < s2))/sum(with(df, misses >= s1 & misses < s2)) })
#plot_df$call_corrected_methylated = sapply(plot_df$misses, FUN = function(x) { s1 = x - binwidth/2; s2 = x+binwidth/2; sum(with(df, call_corrected == 1 & misses >= s1 & misses < s2))/sum(with(df, misses >= s1 & misses < s2)) })
#plot_df$true_methylated = sapply(plot_df$misses, FUN = function(x) { s1 = x - binwidth/2; s2 = x+binwidth/2; sum(with(df, truth == 1 & misses >= s1 & misses < s2))/sum(with(df, misses >= s1 & misses < s2)) })
plot_df = melt(plot_df, id="misses")
pdf(paste(out_filename, "pass_prop_by_misses_by_truth", "pdf", sep="."))
ggplot() +
  geom_line(data=plot_df, aes(x=misses, y=value, colour=variable)) +  
  labs(y="Proportion of sites called correctly",
       x="Proportion of X's")
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
# build methylation on polynomial fit
###########################################
#train_idx = sample(nrow(df), 200000)
df$truth = factor(as.numeric(as.character(df$truth)), levels=c(-1, 1), labels=c(-1, 1))
train_df = df[train_idx,]
test_df = df[-train_idx,]
fit = glm(truth~(pM + pU + pu + pX + px + pD + pd)^6, data=train_df, family=binomial())
train_df$methylation = predict(fit, train_df)
roc_curve = roc(truth~methylation, data=train_df)
threshold = coords(roc_curve, "best", ret="threshold")
test_df$methylation = predict(fit, test_df)
test_df$correct = ifelse((test_df$methylation > threshold & test_df$truth == 1) |
                           (test_df$methylation < threshold & test_df$truth == -1), 1, -1)
glm_phred=-10*log10(sum(test_df$correct==-1)/nrow(test_df))
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
       caption=paste0("Phred = ",glm_phred),
       fill="Truth") +
  #scale_x_continuous(limits=c(0,1)) + 
  scale_fill_manual(labels=c("Unmethyl","Methyl"), values=c("blue","red"))
dev.off()
polynomial_df = test_df

###########################################
# build methylation on GLMNET
###########################################
#train_idx = sample(nrow(df), 200000)
df$truth = factor(as.numeric(as.character(df$truth)), levels=c(-1, 1), labels=c(-1, 1))
train_df = df[train_idx,]
test_df = df[-train_idx,]
glmnet_cv = cv.glmnet(as.matrix(train_df[,c("pM", "pU", "pu", "pX", "px", "pD", "pd")]), train_df$truth, family="binomial")
lambda = glmnet_cv$lambda.1se
glmnet_fit = glmnet_cv$glmnet.fit
train_df$methylation = predict(glmnet_fit, as.matrix(train_df[,c("pM", "pU", "pu", "pX", "px", "pD", "pd")]), s=lambda, type="response")
roc_curve = roc(truth~methylation, data=train_df)
threshold = coords(roc_curve, "best", ret="threshold")
test_df$methylation = predict(glmnet_fit, as.matrix(test_df[,c("pM", "pU", "pu", "pX", "px", "pD", "pd")]), s=lambda, type="response")
test_df$correct = ifelse((test_df$methylation > threshold & test_df$truth == 1) |
                           (test_df$methylation < threshold & test_df$truth == -1), 1, -1)
glmnet_phred=-10*log10(sum(test_df$correct==-1)/nrow(test_df))
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
       caption=paste0("Phred = ",glmnet_phred),
       fill="Truth") +
  #scale_x_continuous(limits=c(0,1)) + 
  scale_fill_manual(labels=c("Unmethyl","Methyl"), values=c("blue","red"))
dev.off()
polynomial_df = test_df

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
       colour="Method")
dev.off()

#####################################################
# checking time lags
#####################################################
ccf((df_forward$methylation-mean(df_forward$methylation)), (df_forward$methylation-mean(df_forward$methylation)), lag.max=40)
ccf((df_reverse$methylation-mean(df_reverse$methylation)), (df_reverse$methylation-mean(df_reverse$methylation)), lag.max=40)

#pdf(paste(out_filename, "correctness_by_position","pdf", sep="."))
rang=108500:109500
#colScale = rainbow(100, start=0.4, end=0.7)
colScale = gray(1:100 / 100)
colCovScale = gray(1:101 / 100)
rowSideCols = colScale[round(df_forward$methylation[rang], digits=2)*100]
colSideCols = colCovScale[round(df_forward$coverage[rang]/max(df_forward$coverage[rang]),digits=2)*100+1]
methbypos = abs(outer(df_forward$methylation[rang], df_forward$methylation[rang], `-`))
colnames(methbypos) = rownames(methbypos) = df_forward$position[rang]
heatmap.2(as.matrix(methbypos), Colv=FALSE, Rowv=FALSE, dendrogram="none", trace="none", main="forward", labRow=NA, labCol=NA,RowSideColors = rowSideCols,ColSideColors = colSideCols)

rowSideCols = colScale[round(df_reverse$methylation[rang], digits=2)*100]
methbypos = abs(outer(df_reverse$methylation[rang], df_reverse$methylation[rang], `-`))
colnames(methbypos) = rownames(methbypos) = df_reverse$position[rang]
heatmap.2(as.matrix(methbypos), Colv=FALSE, Rowv=FALSE, dendrogram="none", trace="none", main="reverse", labRow=NA, labCol=NA,RowSideColors = df_reverse$methylation)

pdf(paste(out_filename, "methylation_to_misses", "pdf", sep="."))
ggplot(df, aes(x=misses, y=methylation)) +
  stat_density2d(geom = "raster", aes(alpha=..density..), fill = "dodgerblue", contour = FALSE) +
  stat_density2d() + 
  geom_smooth(method='lm', formula=y~x) + 
  labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                     "Y-Int =",signif(fit$coef[[1]],5 ),
                     " Slope =",signif(fit$coef[[2]], 5),
                     " P =",signif(summary(fit)$coef[2,4], 5))) + 
  ylab("Methylation percentage") + 
  xlab("Proportion of reads not called")
dev.off()

pdf(paste(out_filename, "misses_cdf", "pdf", sep="."))
ggplot(cdf, aes(x=z)) + 
  geom_line(aes(y=pass, col="Correct call")) + 
  geom_line(aes(y=fail, col="Incorrect call")) +
  ylab("CDF") +
  xlab("Proportion of sites called as neither mC or C")
dev.off()

pdf(paste(out_filename, "misses_pdf", "pdf", sep="."))
ggplot(cdf, aes(x=z)) + 
  geom_line(aes(y=pass_pdf, col="Correct call")) + 
  geom_line(aes(y=fail_pdf, col="Incorrect call")) +
  ylab("PDF") +
  xlab("Proportion of sites called as neither mC or C")
dev.off()

pdf(paste(out_filename, "correctness_to_misses", "pdf", sep="."))
ggplot(df, aes(x=misses)) + geom_bar(group=df$correct)
dev.off()

pdf(paste(out_filename, "scatter", "pdf", sep="."))
ggplot(data=df, aes(x=position, y=ll)) + 
  geom_point(alpha=0.05) + 
  xlab("Genomic Position") +
  scale_y_continuous(limits=c(-4,4)) + 
  ylab("Methylation Log Likelihood Ratio") +
  geom_hline(yintercept = truth[1]*0.7, colour="red") +
  ggtitle(paste("CpG Methylation by Position on E.Coli K12 MG1655", title, sep=" "))
dev.off()

pdf(paste(out_filename, "density", "pdf", sep="."))
ggplot(data=df, aes(x=position, y=ll)) + 
  stat_density2d(geom = "raster", aes(alpha=..density..), fill = "dodgerblue", contour = FALSE) +
  geom_smooth(method='gam', alpha=0.5, colour="dodgerblue4", fill="dodgerblue3", level=0.99) +
  xlab("Genomic Position") +
  scale_y_continuous(limits=c(-4,4)) + 
  ylab("Methylation Log Likelihood Ratio") +
  geom_hline(yintercept = truth[1]*0.7, colour="red") +
  ggtitle(paste("CpG Methylation by Position on E.Coli K12 MG1655", title, sep=" "))
dev.off()

pdf(paste(out_filename, "error_rate", "pdf", sep="."))
ggplot(data=threshold_df, aes(x=thresholds, y=error_rate)) + 
  geom_line() + 
  xlab("Log likelihood ratio threshold") +
  scale_y_continuous(limits=c(0,3)) + 
  ylab("Methylation call error rate (%)") +
  ggtitle(paste("CpG methylation error rate vs threshold\nE.Coli K12 MG1655", title, sep=" "))
dev.off()

pdf(paste(out_filename, "num_calls", "pdf", sep="."))
ggplot(data=threshold_df, aes(x=thresholds, y=num_calls)) + 
  geom_line() + 
  xlab("Log likelihood ratio threshold") +
  scale_y_continuous(limits=c(0,num_sites)) + 
  ylab("Number of calls") +
  ggtitle(paste("CpG methylation number of calls vs threshold\nE.Coli K12 MG1655", title, sep=" "))
dev.off()

pdf(paste(out_filename, "threshold", "pdf", sep="."))
ggplot(data=threshold_df, aes(x=prop_called, y=error_rate)) + 
  geom_line() + 
  xlab("Proportion of CpG sites called") +
  scale_y_continuous(limits=c(0,3)) + 
  scale_x_continuous(limits=c(0,1)) +
  ylab("Methylation call error rate (%)") +
  ggtitle(paste("CpG methylation error rate vs proportion of sites called\nE.Coli K12 MG1655", title, sep=" "))
dev.off()

pdf(paste(out_filename, "threshold", "roc", "pdf", sep="."))
ggplot(data=threshold_df, aes(x=prop_not_called, y=correct_rate)) + 
  geom_line() + 
  xlab("Proportion of CpG sites not called") +
  scale_y_continuous(limits=c(97,100)) + 
  scale_x_continuous(limits=c(0,1)) +
  ylab("Methylation call correctness rate (%)") +
  ggtitle(paste("CpG methylation error rate vs proportion of sites called\nE.Coli K12 MG1655", title, sep=" "))
dev.off()

peakfinder <- function(d, sites=3, binwidth=1000){
  dh <- hist(d,breaks=num_sites%/%binwidth)
  ins <- dh[["counts"]]
  nbins <- length(ins)
  ss <- which(ins%in%sort(ins)[(length(ins)-sites):(length(ins)-1)]) ## pick the top 3 intensities
  dh[["mids"]][ss]
}


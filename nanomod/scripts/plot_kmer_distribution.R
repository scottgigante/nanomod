# show kmer distributions as a function of current
library(rhdf5)
library(ggplot2)
library(reshape2)

#data_dir = "~/python/nanomod/data/lucattini_no_phage"
data_dir = "~/python/nanomod/data/lucattini_no_phage_random_errors"
num_files = 30
original_group = "Analyses/Basecall_1D_000/BaseCalled_template/Events"
new_group = "Analyses/AlignToRef/CurrentSpaceMapped_template/Events"

# subset dataset
dataset = list.files(path=data_dir, pattern = "\\.fast5$", full.names=TRUE)
dataset = dataset[sample(length(dataset), num_files)]

# read file
output_per_file = sapply(dataset, FUN = function(f) {
  original_df = h5read(f, original_group)
  new_df = h5read(f, new_group)
  
  # align data frames
  original_df$alignment = apply(original_df, 1, FUN = function(x) {
    idx = which(signif(as.numeric(new_df$start), digits=3) == signif(as.numeric(x[2]), digits=3) & 
                  signif(as.numeric(new_df$length), digits=3) == signif(as.numeric(x[4]), digits=3) & 
                  as.character(new_df$model_state) == as.character(x[5]) & 
                  signif(as.numeric(new_df$weights), digits=3) == signif(as.numeric(x[7]), digits=3) & 
                  signif(as.numeric(new_df$p_A), digits=3) == signif(as.numeric(x[11]), digits=3) & 
                  signif(as.numeric(new_df$p_C), digits=3) == signif(as.numeric(x[12]), digits=3) & 
                  signif(as.numeric(new_df$p_G), digits=3) == signif(as.numeric(x[13]), digits=3))
    if (length(idx) > 1) print("multiple match, panic")
    ifelse(length(idx) == 1, idx, NA)
  })
  original_df = original_df[which(!is.na(original_df$alignment)),]
  if (all(sort(unique(original_df$alignment)) == original_df$alignment)) {
    aligned_df = new_df[original_df$alignment,c("mean", "stdv", "kmer", "length")]
    aligned_df$original_mean = original_df$mean
    aligned_df$original_stdv = original_df$stdv
    raw = do.call("c",apply(aligned_df, 1, function(x) {
      length = as.numeric(as.character(x[4]))
      original_mean = as.numeric(as.character(x[5]))
      rep(original_mean, length*4000)
    }))
    med = median(raw)
    MAD = median(abs(med - raw))
    aligned_df$new_mean = (aligned_df$original_mean - med)/MAD
    aligned_df$filename = f
  } else {
    aligned_df = matrix(ncol=0, nrow=0)
  }
  t(aligned_df)
})

# combine dataframes
output_per_file = sapply(output_per_file[-which(sapply(output_per_file, is.null))], FUN = t)
output = as.data.frame(do.call("rbind", output_per_file))
output$mean = as.numeric(as.character(output$mean))
output$stdv = as.numeric(as.character(output$stdv))
output$original_mean = as.numeric(as.character(output$original_mean))
output$original_stdv = as.numeric(as.character(output$original_stdv))
output$new_mean = as.numeric(as.character(output$new_mean))
output$kmer = as.character(output$kmer)

means_df = melt(output, id.vars=c("kmer", "filename"), measure.vars=c("original_mean", "new_mean", "mean"))
means_df$idx = NA
ordering = dcast(means_df[which(means_df$variable=="original_mean"),], kmer~variable, median)
ordering$idx = order(ordering$original_mean)
means_df[which(means_df$variable=="original_mean"),]$idx = sapply(means_df[which(means_df$variable=="original_mean"),"kmer"], function(k) {  which(ordering$idx==which(ordering$kmer==k)) })
ordering = dcast(means_df[which(means_df$variable=="mean"),], kmer~variable, median)
ordering$idx = order(ordering$mean)
means_df[which(means_df$variable=="mean"),]$idx = sapply(means_df[which(means_df$variable=="mean"),"kmer"], function(k) { which(ordering$idx==which(ordering$kmer==k)) })
ordering = dcast(means_df[which(means_df$variable=="new_mean"),], kmer~variable, median)
ordering$idx = order(ordering$new_mean)
means_df[which(means_df$variable=="new_mean"),]$idx = sapply(means_df[which(means_df$variable=="new_mean"),"kmer"], function(k) { which(ordering$idx==which(ordering$kmer==k)) })


ggplot(data=means_df, aes(x=idx, y=value, color=filename)) +
  facet_grid(variable~., scales="free") +
  geom_point(size=0, alpha=0.1, show.legend=FALSE) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

stdvs_df = melt(output, id.vars=c("kmer", "filename"), measure.vars=c("original_stdv", "stdv"))
stdvs_df$idx = NA
ordering = dcast(stdvs_df[which(stdvs_df$variable=="original_stdv"),], kmer~variable, median)
ordering$idx = order(ordering$original_stdv)
stdvs_df[which(stdvs_df$variable=="original_stdv"),]$idx = sapply(stdvs_df[which(stdvs_df$variable=="original_stdv"),"kmer"], function(k) {  which(ordering$idx==which(ordering$kmer==k)) })
ordering = dcast(stdvs_df[which(stdvs_df$variable=="stdv"),], kmer~variable, median)
ordering$idx = order(ordering$stdv)
stdvs_df[which(stdvs_df$variable=="stdv"),]$idx = sapply(stdvs_df[which(stdvs_df$variable=="stdv"),"kmer"], function(k) { which(ordering$idx==which(ordering$kmer==k)) })

ggplot(data=stdvs_df, aes(x=idx, y=value, color=filename)) +
  facet_grid(variable~., scales="free") +
  geom_point(size=0, alpha=0.1, show.legend=FALSE) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

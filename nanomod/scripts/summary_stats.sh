#!/bin/sh
# reads stdin and prints summary statistics
# total, count, mean, median, min and max
# you can pipe this through scripts or redirect input <
# modified from http://unix.stackexchange.com/a/13779/27194

# Usage: python seq_length.py in.fasta | cut -f 2 | bash summary_stats.sh
# Usage: samtools view 116.sorted.bam | gawk 'match($6, /([0-9]*)S(.*)/, a) match($6, /.*[A-Z]([0-9]*)S/, b) {print length($10) - a[1] - b[1]}' | bash summary_stats.sh
sort -n | awk '
  $1 ~ /^[0-9]*(\.[0-9]*)?$/ {
    a[c++] = $1;
    sum += $1;
  }
  END {
    ave = sum / c;
    if( (c % 2) == 1 ) {
      median = a[ int(c/2) ];
    } else {
      median = ( a[c/2] + a[c/2-1] ) / 2;
    }
    OFS="\t";
    { printf ("Total:\t""%'"'"'d\n", sum)}
    { printf ("Count:\t""%'"'"'d\n", c)}
    { printf ("Mean:\t""%'"'"'d\n", ave)}
    { printf ("Median:\t""%'"'"'d\n", median)}
    { printf ("Min:\t""%'"'"'d\n", a[0])}
    { printf ("Max:\t""%'"'"'d\n", a[c-1])}
  }
'

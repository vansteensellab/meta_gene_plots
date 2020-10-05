#!/usr/bin/env Rscript

library(ggplot2)
library(data.table)
library(argparse)


parser <- ArgumentParser(description='Process count matrix into normalized running means')

parser$add_argument('--signal', '-s', help='Matrix with counts of signal from deeptools')
parser$add_argument('--binsize', '-bs', type='integer', help='binsize used in the matrix')
parser$add_argument('-a', type='integer', help='number of bp after POI')
parser$add_argument('-b', type='integer', help='number of bp before POI')
parser$add_argument('--body', type='integer', help='body length')
parser$add_argument('--scaled', type='logical', help='scaled region around gene body instead of POI',
                    default=F)
parser$add_argument('--unscaled3', '-u3', type='integer', help="unscaled region on 3'")
parser$add_argument('--unscaled5', '-u5', type='integer', help="unscaled region on 5'")
parser$add_argument('--bed', help="bed-file with POI's or regions")
parser$add_argument('--output', '-o', help="output file")
parser$add_argument('-k', type='integer', help='number of bins around center', default=20)
parser$add_argument('--min_count', '-m', type='integer',
                    help='minimum read count in run sum', default=10)
parser$add_argument('--step', type='integer', help='number of bins each step', default=2)
parser$add_argument('--exp', help='experiment name')
parser$add_argument('--ctrl', help='control name')

argv <- parser$parse_args()

save(argv, file='test.Rdata')

## assuming k is even, we can get a "running" sum for k bins around bin i.
## example with k = 10:
##
## [ ][ ][ ][#][#][#][#][#]|[i][#][#][#][#][ ][ ]
##
## center: |
## bin i: [i]
## selected bin: [#]
##
runsum_bins <- function(dt, k, step){
    runsum_list = lapply(seq(1,ncol(dt),step), function(i){
        if (i - k/2 < 1){
            start = 1
            end = k
        } else if (i + k/2 > ncol(dt)){
            start = ncol(dt) - k + 1
            end = ncol(dt)
        } else {
            start = i - k/2
            end = i + k/2 - 1
        }
        rowSums(dt[, start:end], na.rm=T)
    })
    data.table(do.call(cbind, runsum_list))
}


# P = fread(argv$bed, stringsAsFactors=F)
signal_dt = fread(cmd=paste('zcat', argv$signal), stringsAsFactors=F, skip=1)
signal_runsum = runsum_bins(signal_dt[,-c(1:6)], argv$k, argv$step)



nuc_step = argv$binsize * argv$step
if (argv$scaled){
    tss_pos_vec = ((-argv$b / nuc_step) : (argv$unscaled5 / nuc_step - 1)) * nuc_step
    tes_pos_vec = ((-argv$unscaled3 / nuc_step) : (argv$a / nuc_step - 1)) * nuc_step
    mid_pos_vec = 0 : (argv$body / nuc_step - 1) * nuc_step
    pos_vec = c(paste0('TSS', ifelse(tss_pos_vec < 0, '', '+'), tss_pos_vec),
                paste0('MID', '+', mid_pos_vec),
                paste0('TES', ifelse(tes_pos_vec < 0, '', '+'), tes_pos_vec))
} else {
    pos_vec = as.character(((-argv$b / nuc_step) : (argv$a / nuc_step - 1)) * nuc_step
                           + nuc_step/2)
}

colnames(signal_runsum) = pos_vec

norm_dt = cbind(signal_dt[,4], signal_runsum)
colnames(norm_dt)[1] = 'transcript_ID'

fwrite(norm_dt, argv$out, sep='\t')

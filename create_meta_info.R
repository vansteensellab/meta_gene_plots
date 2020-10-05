
library(data.table)
library(stringr)
dam_list = list.files(paste0('/DATA/usr/t.v.schaik/proj/sManzo_pADamID/',
                             'ts190515_pADamID_RPE_Top1_DRB/results/counts/bin-gatc'),
                      full.names=T)

meta_dt = data.table(str_match(dam_list, '.*/pADamID-(.*)_(.*)-gatc.counts.*'))
colnames(meta_dt) = c('file', 'exp_id','target')

control_dt = meta_dt[target=='Dam', ]
exp_dt = meta_dt[target!='Dam', ]

exp_dt[, name:=paste0(exp_id, '_', target)]
exp_dt[, signal:=file]

setkey(control_dt, exp_id)

exp_dt[, control:=control_dt[exp_id, 'file']]
exp_dt[, type:='dam']


fwrite(exp_dt[,c('name', 'type', 'signal', 'control')], 'dam_files_TOP_DRB.txt',
       sep='\t')

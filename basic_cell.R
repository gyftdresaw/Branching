

# Some plots/investigation of basic cell simulation
#  -- arbitrary waiting time distribution for division
#  -- produce n progeny per division

setwd('Documents/Dinner/Branching/')

# extract exponential data file
exp1_p2.df = read.table('results/basic_ncell_exp1_p2.txt',sep='\t')
exp1_p2.df = as.data.frame(t(exp1_p2.df))
exp1_p2.df = exp1_p2.df[-dim(exp1_p2.df)[1],] # remove last row of NA

# other exponential mean data
exp075_p2.df = read.table('results/basic_ncell_exp075_p2.txt',sep='\t')
exp075_p2.df = as.data.frame(t(exp075_p2.df))
exp075_p2.df = exp075_p2.df[-dim(exp075_p2.df)[1],] # remove last row of NA

# other exponential mean data
exp045_p3.df = read.table('results/basic_ncell_exp045_p3.txt',sep='\t')
exp045_p3.df = as.data.frame(t(exp045_p3.df))
exp045_p3.df = exp045_p3.df[-dim(exp045_p3.df)[1],] # remove last row of NA

## some gamma distribution data
gam5_02_p2.df = read.table('results/basic_ncell_gam5_02_p2.txt',sep='\t')
gam5_02_p2.df = as.data.frame(t(gam5_02_p2.df))
gam5_02_p2.df = gam5_02_p2.df[-dim(gam5_02_p2.df)[1],] # remove last row of NA

## other gamma distribution data
gam6_02_p3.df = read.table('results/basic_ncell_gam6_02_p3.txt',sep='\t')
gam6_02_p3.df = as.data.frame(t(gam6_02_p3.df))
gam6_02_p3.df = gam6_02_p3.df[-dim(gam6_02_p3.df)[1],] # remove last row of NA

require(ggplot2)
require(reshape2)

## exponential mean 1 ##
d = melt(exp1_p2.df,id.vars='V1')
# plot all the growth curves on one plot
g = ggplot(d,aes(V1,value,col=variable))
g = g + geom_line()
g = g + stat_function(fun=function(x) {x/log(10)},geom='line',colour='black',size=2,linetype='longdash')
g = g + scale_x_continuous(limits=c(0,14)) + scale_y_log10(limits=c(1,1e7))
g = g + xlab('Time') + ylab('Num Cells')
g = g + ggtitle(expression(paste('exponential ',lambda,'=1.0 ',nu,'=2',sep='')))
g = g + theme(legend.position='none')
g

ggsave('plots/ncells_exp1_p2.png')

d = melt(exp075_p2.df,id.vars='V1')
# plot all the growth curves on one plot
g = ggplot(d,aes(V1,value,col=variable))
g = g + geom_line()
g = g + stat_function(fun=function(x) {0.75*x/log(10)},geom='line',colour='black',size=2,linetype='longdash')
g = g + scale_x_continuous(limits=c(0,14)) + scale_y_log10(limits=c(1,1e7))
g = g + xlab('Time') + ylab('Num Cells')
g = g + ggtitle(expression(paste('exponential ',lambda,'=0.75 ',nu,'=2',sep='')))
g = g + theme(legend.position='none')
g

ggsave('plots/ncells_exp075_p2.png')

d = melt(exp045_p3.df,id.vars='V1')
# plot all the growth curves on one plot
g = ggplot(d,aes(V1,value,col=variable))
g = g + geom_line()
g = g + stat_function(fun=function(x) {2*0.45*x/log(10)},geom='line',colour='black',size=2,linetype='longdash')
g = g + scale_x_continuous(limits=c(0,14)) + scale_y_log10(limits=c(1,1e7))
g = g + xlab('Time') + ylab('Num Cells')
g = g + ggtitle(expression(paste('exponential ',lambda,'=0.45 ',nu,'=3',sep='')))
g = g + theme(legend.position='none')
g

ggsave('plots/ncells_exp045_p3.png')

## gamma distribution plots ##

k_gam = (2^(1/5) - 1)/0.2
d = melt(gam5_02_p2.df,id.vars='V1')
# plot all the growth curves on one plot
g = ggplot(d,aes(V1,value,col=variable))
g = g + geom_line()
g = g + stat_function(fun=function(x) {k_gam*x/log(10)},geom='line',colour='black',size=2,linetype='longdash')
g = g + scale_x_continuous(limits=c(0,14)) + scale_y_log10(limits=c(1,1e7))
g = g + xlab('Time') + ylab('Num Cells')
g = g + ggtitle(expression(paste('gamma ',alpha,'=5 ',beta,'=0.2 ',nu,'=2',sep='')))
g = g + theme(legend.position='none')
g

ggsave('plots/ncells_gam5_02_p2.png')

k_gam = (3^(1/6) - 1)/(1/5)
d = melt(gam6_02_p3.df,id.vars='V1')
# plot all the growth curves on one plot
g = ggplot(d,aes(V1,value,col=variable))
g = g + geom_line()
g = g + stat_function(fun=function(x) {k_gam*x/log(10)},geom='line',colour='black',size=2,linetype='longdash')
g = g + scale_x_continuous(limits=c(0,14)) + scale_y_log10(limits=c(1,1e7))
g = g + xlab('Time') + ylab('Num Cells')
g = g + ggtitle(expression(paste('gamma ',alpha,'=6 ',beta,'=0.2 ',nu,'=3',sep='')))
g = g + theme(legend.position='none')
g

ggsave('plots/ncells_gam6_02_p3.png')





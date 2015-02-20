

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







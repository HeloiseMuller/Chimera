#####################################
#### Read arguments and files #######
######################################

args = commandArgs(trailingOnly=FALSE)
scriptPath <- dirname(sub("--file=","", args[grep("--file", args)]))
print(scriptPath)

args = commandArgs(trailingOnly=TRUE)
for(i in 2:length(args))
{
    eval(parse(text=sub("--", "", args[i])))
}

help = FALSE
if(!exists("config_file")){help = TRUE; print("config_file is missing")}
if(!exists("depth_dir")){help = TRUE; print("depth_dir is missing")}

if(help==TRUE)
{
    print("1-Coverage.R --args --config_file=\"config.Rdata\" --depth_dir=\"directory/\"")
    quit("no")
}

require(data.table)
require(dplyr)
require(stringr)
library(data.table)

depth=TRUE #since true, need to fill the variable `Coverage` in the config file
IPS=FALSE
source(config_file)
source(paste0(scriptPath, "/functions.Rdata"))

#### Inputs ###

cov = lapply(paste0(depth_dir, samples, "_CtBV.cov"),fread)


### Main ###
names = c("Caterpillar", "Adult-1", "Adult-6", "PoolF1-A", "PoolF1-B", "PoolF1-C", "PoolF1-D", "PoolF1-E", expression(italic("N. typhae")), expression(italic("G. sparganii")), expression(italic("C. phragmitella")), "Control")
svg(paste0(dir,"/Figures/Figure2.svg"),  width = 8, height = 5)    
par(xpd=NA) #to be allowed to draw outside figure
layout(matrix(c(1,1,2), 1, 3, byrow=TRUE))
barplot(t(as.matrix(Coverage[,3:5])), beside=T, names.arg=names,
    cex.names = 0.7,
    legend.text= c(species[1], species[2], "CtBV"),
    args.legend = list(bty = "n", x = "top", ncol=3, inset=0.1), 
    col = c("darkgreen", "gold", "red3"), 
    ylim = c(0, ceiling(max(Coverage$depth1, Coverage$depth2, Coverage$depthBV)*1.2)), 
    ylab = "Depth (X)",
    las=2,
    main = "(a)")
segments(2,-45,32,-45, lwd=0.5)
text(15,-52, labels=expression(italic("S. nonagrioides")), cex=0.8)

barplot(t(as.matrix(Coverage[,6:7])), beside=T, names.arg=names, legend.text=F, col = c("darkgreen", "gold"), ylim = c(0,1), ylab = "Proportion of the genome covered", cex.names=0.7, las=2, main="(b)")
segments(1,-0.145,24,-0.145, lwd=0.5)
text(12,-0.165, labels=expression(italic("S. nonagrioides")), cex=0.8)

dev.off()

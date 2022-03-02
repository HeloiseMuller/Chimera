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

png(paste0(dir,"/Figures/Coverage.png"),  width = 8, height = 5, units = 'in', res = 600)    
layout(matrix(c(1,1,2), 1, 3, byrow=TRUE))
barplot(t(as.matrix(Coverage[,3:5])), beside=T, names.arg=Coverage$Names,
    legend.text= c(species[1], species[2], "Proviral segments"), 
    col = c("darkgreen", "gold", "red3"), 
    ylim = c(0, ceiling(max(Coverage$depth1, Coverage$depth2, Coverage$depthBV)*1.3)), 
    ylab = "Depth (X)",
    las=2)

barplot(t(as.matrix(Coverage[,6:7])), beside=T, names.arg=Coverage$Names, legend.text=F, col = c("darkgreen", "gold"), ylim = c(0,1), ylab = "Proportion of the genome covered", cex.names=0.7, las=2)
dev.off()

##Coverage for each segments

#Edit colnames 
cov = lapply(cov, function(x) {
   colnames(x) <- c("Contig","Start","End", "Segment", "Depth")
   return(x)
  })

#### Coverage of each segments proportionally to coverage on C. typhae

#Arrange segments in order
cov = lapply(cov, fSeg)
for (i in seq_along(cov)){
    cov[[i]] = arrange(cov[[i]], -Seg)
    cov[[i]]$Sample = samples[i]
    }

#Relatice depth. Here, one table with all depths
cov_seg = mutate(cov[[1]], relativeDepth = Depth-Coverage[1,]$depth2) %>% select(., Segment, relativeDepth, Seg, Sample) 
for (i in 2:length(cov)){
    cov_seg = rbind(cov_seg, mutate(cov[[i]], relativeDepth = Depth-Coverage[i,]$depth2) %>% select(., Segment, relativeDepth, Seg, Sample) )
    }
cov_seg = mutate(cov_seg, relativeDepth = ifelse(relativeDepth<0, 0, relativeDepth))
print(filter(cov_seg, Sample=="NH212C"))

png(paste0(dir,"/Figures/Coverage_segments.png"),  width = 10, height = 5, units = 'in', res = 600)     
par(mar=c(4, 0, 1, 6), xpd=TRUE)
dotchart(cov[[1]]$Depth-Coverage[1,]$depth2, 
    labels= cov[[1]]$Segment, cex = 0.8,
    bg = colors[1], pt.cex = 0.8, 
    xlim = c(0, ceiling(max(cov_seg$relativeDepth))*1.5), 
    xlab = expression(paste("Segment depths minus the depth on ", italic('C. typhae'))))
for (i in 2:length(cov)){
    points(cov[[i]]$Depth-Coverage[i,]$depth2, 1:nrow(cov[[1]]), col = colors[i], pch = 19, cex = 0.8)
}
legend("topright", inset=c(-0.15, .2), legend = names, col = colors, fill = colors, pch=19, cex = 0.8)
dev.off()


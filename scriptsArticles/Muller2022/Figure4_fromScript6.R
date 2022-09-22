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
if(!exists("meta")){help = TRUE; print("meta is missing")}

if(help==TRUE)
{
    print("6-NbIEvsDepth.R --args --config_file=\"config.Rdata\" --depth_dir=\"directory/\" --meta=\"Table_CtBV.txt\"")
    quit("no")
}

require(data.table)
require(dplyr)
require(stringr)
require(ggplot2)

depth=TRUE
IPS=FALSE
source(config_file)
source(paste0(scriptPath, "/functions.Rdata"))

#### Inputs ###
print("Reading depths...")
cov = lapply(paste0(depth_dir, samples, "_CtBV.cov"),fread)

cov = lapply(cov, function(x) {
   colnames(x) <- c("Contig","Start","End", "Segment", "Depth")
   return(x)
  })

print("Reading chimera..")
#We also want the chimeria from PCR dup, since we have them in depth
chimera_alongSegments = lapply(paste0(dir, "/Chimera/", samples, "_all_chimera_alongSegments.txt"), fread, header = T, sep = "\t")
#ALso reads the file without PCR dup that count the number of IE
chimera_alongSegments_IE = lapply(paste0(dir, "/Chimera/", samples, "_all_chimera_alongSegments_IE.txt"), fread, header = T, sep = "\t")

print("Reading metadata...")
#Read metadata (mostly to know who has HIM)
meta=fread(meta)

#Summarie infos
summarize_chimera = lapply(chimera_alongSegments_IE, function(x){
    i = which(samples==unique(x$sample))
    x = x %>% group_by(sample, Segment) %>% summarise(NbChimera_noPCRdup = sum(NbChimera), NbIE = n())
    x = chimera_alongSegments[[i]]  %>% group_by(sample, Segment) %>% summarise(NbChimera=n()) %>% left_join(., x, by=c("sample", "Segment"))
    })
#### Functions ###
  
#Function to merge tables coverage and summarize_chimera with the colomns of interest
fMerge = function(x,y){
    dt = select(x, c("Depth", "Seg", "Segment")) %>% left_join(., select(y, -c("sample", "Segment")), by="Seg")
    return(dt)
    }
  
##Correlation between coverage segments and NbIE ?

#Remove line corresponding to duplicated segments
cov = lapply(cov, function(x) {
   x <- x[!str_detect(x$Segment, "Duplication"),]
   return(x)
  })

#Make a colomn in both table where segment names have the same format
cov = lapply(cov, fSeg)
summarize_chimera = lapply(summarize_chimera, fSeg)
meta = fSeg(meta)

#Merge tables cov and summarize on chimera
summarize_chimera_depth = vector(mode = "list", length = length(samples))
for (i in 1:length(samples)){
    summarize_chimera_depth[[i]] = fMerge(cov[[i]], summarize_chimera[[i]])
    summarize_chimera_depth[[i]]$sample = samples[i]
    }

#Panels for article
summarize_chimera_depth = summarize_chimera_depth[-1]
panels = c("(a)", "(b)", "(c)", "(d)", "(e)")
namesN = c("Adult-1", "Adult-6", expression(italic(" N. typhae")), expression(italic("G. sparganii")), expression(italic("C. phragmitella")))
svg(paste0(dir, "../FiguresArticle/Figure4.svg"),  width = 10, height = 6) 
#layout(matrix(c(1,1,2,2,3,3,4,4,0,5,5,0), nrow=3, byrow=TRUE))
layout(matrix(c(0,1,1,2,2,0,3,3,4,4,5,5), nrow=2, byrow=TRUE))
par(xpd=T, mai=c(0.4,0.5,0.5,0)) #to be able to draw segments under the graph
lapply(seq_along(summarize_chimera_depth), function(i){
	       x = summarize_chimera_depth[[i]]
	       x[is.na(x$NbChimera),]$NbChimera=0
	       pt=plot(x$NbChimera~x$Depth, ylab=ifelse(i==2 | i==4 | i==5, "", "Number of chimeric reads"), 
		       xlab="", #set it with title() to bring it closer axis 
		       pch=ifelse(i==2 & (x$Seg==20 | x$Seg==23), 18, 20),
		       col=ifelse(x$Seg %in% filter(meta, !is.na(HIM_start))$Seg, "blue","red"),
		       cex=1.1,
		       xlim=c(0, max(x$Depth)),
		       main =  panels[i])
	       	mtext(namesN[i], side=3)
	       title(xlab="Depth", line=2, cex.lab=1.2)
if(i==2){
    segments(filter(summarize_chimera_depth[[2]], Seg == 20)$Depth, 0, filter(summarize_chimera_depth[[2]], Seg == 20)$Depth, -10, lty=2, lwd=0.8)
    text(filter(summarize_chimera_depth[[2]], Seg == 20)$Depth, -13, "S20/33", cex=0.75, col="red", font=2)
    segments(filter(summarize_chimera_depth[[2]], Seg == 23)$Depth, 0, 15, -10, lty=2, lwd=0.8)
    text(15, -13, "S23", cex=0.75, col="red", font=2)
}
})

dev.off()

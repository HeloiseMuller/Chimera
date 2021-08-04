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
    print("6-NbIEvsDepth.R --args --config_file=\"config.Rdata\" --depth_dir=\"directory/\"")
    quit("no")
}

require(data.table)
require(dplyr)
require(stringr)
library(data.table)

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

print("Reading IEs..")
chimera_alongSegments_IE = lapply(paste0(dir, "/Chimera/", samples, "_all_chimera_alongSegments_IE.txt"), fread, header = T, sep = "\t")
IE = lapply(chimera_alongSegments_IE, function(x){
    i = which(samples==unique(x$sample))
    x = x %>% group_by(sample, Segment) %>% summarise(NbChimera = sum(NbChimera), NbIE = n())
    })
    
#### Functions ###
  
#Function to merge tables coverage and NbIE with the colomns of interest
fMerge = function(x,y){
    dt = select(x, c("Depth", "Seg")) %>% left_join(.,y[, c("NbChimera", "Seg")], by="Seg")
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
IE = lapply(IE, fSeg)

#Merge tables cov and NbIE
dt = vector(mode = "list", length = length(samples))
for (i in 1:length(samples)){
    dt[[i]] = fMerge(cov[[i]], IE[[i]])
    dt[[i]]$Sample = samples[i]
    }


#Check if correlation among lines without NA
dt_cor = lapply(dt, function(x) {
   x <- x[!is.na(x$NbChimera),]
   print(cor.test(x$NbChimera, x$Depth, method = "spearman"))
   return(x)
  })
 
#Pool samples
dt_tot = bind_rows(dt) %>%  group_by(Seg) %>% summarise(Depth = sum(Depth), NbChimera = sum(NbChimera))


#Keep only segments with IE for correlation
dt_cor_tot = dt_tot[!is.na(dt_tot$NbChimera),]
cor.test(dt_cor_tot$NbChimera, dt_cor_tot$Depth, method = "spearman")


pdf(paste0(dir, "/Figures/CoverageSegments_NbReads.pdf"))
dt_tot[is.na(dt_tot$NbChimera),]$NbChimera=0 #Replace NA by 0 to see it on plot
plot(dt_tot$NbChimera~dt_tot$Depth, ylab="Number of chimeric reads", xlab="Depth", pch=20,
    col=ifelse(dt_tot$NbChimera==0, "red","blue"),
    main = "Sum up all samples")
text(dt_tot$Depth, dt_tot$NbChimera+2, labels = dt_tot$Seg, cex=0.5)
text(500,1900, paste0("rho = ", round(cor.test(dt_cor_tot$NbChimera, dt_cor_tot$Depth, method = "spearman")$estimate, 3)), cex=1)
#abline(v=sum(Coverage[1:length(samples),]$depth2), lwd=1, lty=2, col="gold")


lapply(dt, function(x){
    sample = unique(x$Sample)
    print(sample)
    x[is.na(x$NbChimera),]$NbChimera=0
    plot(x$NbChimera~x$Depth, ylab="Number of chimeric reads", xlab="Depth", pch=20, 
    col=ifelse(x$NbChimera==0, "red","blue"), xlim=c(0, max(x$Depth)), main = Coverage[which(Coverage$Sample==sample),]$Names)
    text(x$Depth, x$NbChimera+max(x$NbChimera)/37, labels = x$Seg, cex=0.5)
  #  abline(v=Coverage[which(Coverage$Sample==sample),]$depth2, lwd=1, lty=2, col="gold")
})


dev.off()




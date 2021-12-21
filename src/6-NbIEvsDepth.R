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

print("Reading IEs..")
chimera_alongSegments_IE = lapply(paste0(dir, "/Chimera/", samples, "_all_chimera_alongHIMs_IE.txt"), fread, header = T, sep = "\t")
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


##Barplot to compare NbRead and Cov in each sample 

#Calcutae ratio between cov and NbIE
fratio = function(x){
    x$ratio = x$Depth/x$NbChimera
    x= setorder(x, NbChimera)
    return(x)
}

dt_cor = lapply(dt_cor, fratio)

png(paste0(dir, "/Figures/Mirror_Depth_NbRead.png"),  width = 10, height = 10, units = 'in', res = 600)    
#Figure in mirror Cov & NbReads with ratio on top of each bars

#function to round to the upper accuracy
round_any = function(x, accuracy, f=round){f(x/ accuracy) * accuracy}


layout(matrix(c(1,1,1,1,2,2,2,2,3,3,3,3,0,0,4,4,4,4,5,5,5,5,0,0), 2, 3*4, byrow=TRUE))
lapply(seq_along(dt_cor), function(i) {
    x = setorder(dt_cor[[i]], NbChimera)
    pt = barplot(x$NbChimera, names.arg=x$Seg, ylim=c(-round_any(max(x$Depth), ifelse(max(x$Depth)/100>1, 100, 10), f=ceiling),round_any(max(x$NbChimera),  ifelse(max(x$Depth)/100>1, 100, 10), f=ceiling)+30), cex.names=0.5, las=2,  col="gray", border=F, yaxt='n', main = names[i])
    barplot(-x$Depth, names.arg=NA, ylim=c(-round_any(max(x$Depth),  ifelse(max(x$Depth)/100>1, 100, 10), f=ceiling),round_any(max(x$NbChimera),  ifelse(max(x$Depth)/100>1, 100, 10), f=ceiling)+30), , add=T, las=2, col="grey43", border=F, yaxt='n')
    text(3, max(x$NbChimera)*0.9, paste0("rho = ", round(cor.test(x$NbChimera, x$Depth, method = "spearman")$estimate, 3)), cex=0.9)
    text(pt, x$NbChimera+15 , round(fratio(x)$ratio,1), cex=0.6)
    axis(side=2, at=seq(-round_any(max(x$Depth),  ifelse(max(x$Depth)/100>1, 100, 10), f=ceiling), round_any(max(x$NbChimera),  ifelse(max(x$Depth)/100>1, 100, 10), f=ceiling), ifelse(max(x$Depth)/100>1, 100, 10)), labels=gsub('-','',seq(-round_any(max(x$Depth),  ifelse(max(x$Depth)/100>1, 100, 10), f=ceiling), round_any(max(x$NbChimera),  ifelse(max(x$Depth)/100>1, 100, 10), f=ceiling), ifelse(max(x$Depth)/100>1, 100, 10))), las=2, cex.axis=0.5)
    return(x)
})

dev.off()  

##Boxplot
dt_box = bind_rows(dt_cor)
dt_box$Sample <- factor(dt_box$Sample, levels = samples, ordered = T)

png(paste0(dir, "/Figures/Boxplot_Depth_NbRead.png"),  width = 7, height = 5, units = 'in', res = 600)    
dt_box %>% ggplot(aes(x=Sample, y=ratio, fill=Sample)) +
    geom_boxplot(fill = colors, outlier.shape = NA) +
    geom_jitter(size=0.7) +  
    scale_x_discrete(labels = names) +
    theme(legend.position = "none") +
    xlab("") +  ylab("Ratio depth/chimera") +
    geom_hline(yintercept = 0.63, linetype = "dashed")
dev.off()





test = lapply(dt_cor, function(x){
    i = which(samples==unique(x$Sample))
    #d = filter(Coverage, Sample==unique(x$Sample))$depth2
    d = filter(dt[[i]], Seg == 15)$Depth
    x$DepthCorrected = x$Depth-d+0.43 #depth on control
    x$ratioCorrected = x$DepthCorrected / x$NbChimera
    return(x)
    })
    
t = lapply(seq_along(test), function(i){
    colnames(test[[i]]) = c("Seg", paste0("ratioCorrected_",samples[i])) 
    return(test[[i]])
    })
t = Reduce(function(x, y) merge(x, y, by = "Seg", all = TRUE), t)    
    
test = bind_rows(test)
test$Sample <- factor(test$Sample, levels = samples, ordered = T)
  
png(paste0(dir, "/Figures/Boxplot_Depth_NbRead_Corrected.png"),  width = 7, height = 5, units = 'in', res = 600)    
test %>% ggplot(aes(x=Sample, y=ratioCorrected, fill=Sample)) +
    geom_boxplot(fill = colors, outlier.shape = NA) +
    geom_jitter(size=0.7) +  
    scale_x_discrete(labels = names) +
    theme(legend.position = "none") +
    xlab("") +  ylab("Ratio depth/chimera") +
    geom_hline(yintercept = 0.63, linetype = "dashed")
dev.off()




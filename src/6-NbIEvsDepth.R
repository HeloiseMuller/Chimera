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
    dt = select(x, c("Depth", "Seg")) %>% left_join(., select(y, -c("sample")), by="Seg")
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
    summarize_chimera_depth[[i]]$Sample = samples[i]
    }

#Check if correlation for HIM-containing segments
print("Corelation between number of chimera (including PCR dup) of HIM-containing BV and sequencing depth")
dt_cor = lapply(summarize_chimera_depth, function(x) {
   x = filter(x, Seg %in% filter(meta, !is.na(HIM_start))$Seg)
   print(cor.test(x$NbChimera, x$Depth, method = "spearman"))
   return(x)
  })
 
#Pool samples
summarize_chimera_depth_sumSamples = bind_rows(summarize_chimera_depth) %>%  group_by(Seg) %>%
    summarise(Depth = sum(Depth), NbChimera = sum(NbChimera, na.rm=T), NbChimera_noPCRdup = sum(NbChimera_noPCRdup, na.rm=T), NbIE = sum(NbIE, na.rm=T))

#Keep only segments with HIM for correlation
dt_cor_sumSamples = filter(summarize_chimera_depth_sumSamples, Seg %in% filter(meta, !is.na(HIM_start))$Seg)
print("Corelation when sum up all samples")
cor.test(dt_cor_sumSamples$NbChimera, dt_cor_sumSamples$Depth, method = "spearman")

#Plot the number of chimera (with PCR dups as a function of the sequencing depth)
pdf(paste0(dir, "/Figures/CoverageSegments_NbReads.pdf"))

summarize_chimera_depth_sumSamples[is.na(summarize_chimera_depth_sumSamples$NbChimera),]$NbChimera=0 #Replace NA by 0 to see it on plot
plot(summarize_chimera_depth_sumSamples$NbChimera~summarize_chimera_depth_sumSamples$Depth, ylab="Number of chimeric reads", xlab="Depth", pch=20,
    col=ifelse(summarize_chimera_depth_sumSamples$Seg %in% filter(meta, !is.na(HIM_start))$Seg, "blue","red"),
    main = "Sum up all samples")
text(summarize_chimera_depth_sumSamples$Depth, summarize_chimera_depth_sumSamples$NbChimera+2, labels = summarize_chimera_depth_sumSamples$Seg, cex=0.5)
text(500,1900, paste0("rho = ", round(cor.test(dt_cor_sumSamples$NbChimera, dt_cor_sumSamples$Depth, method = "spearman")$estimate, 3)), cex=1)
#abline(v=sum(Coverage[1:length(samples),]$depth2), lwd=1, lty=2, col="gold")


lapply(summarize_chimera_depth, function(x){
    sample = unique(x$Sample)
    x[is.na(x$NbChimera),]$NbChimera=0
    plot(x$NbChimera~x$Depth, ylab="Number of chimeric reads", xlab="Depth", pch=20, 
    col=ifelse(x$Seg %in% filter(meta, !is.na(HIM_start))$Seg, "blue","red"), 
    xlim=c(0, max(x$Depth)), main = Coverage[which(Coverage$Sample==sample),]$Names)
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


layout(matrix(c(1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4,5,5,5,5,6,6,6,6), 2, 3*4, byrow=TRUE))
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
    xlab("") +  ylab("Ratio depth/chimera")
dev.off()


#Correct for sequencing depth on proviral segmets

dt_cor_corrected = lapply(dt_cor, function(x){
    i = which(samples==unique(x$Sample))
    d = filter(Coverage, Sample==unique(x$Sample))$depth2
    #d = filter(summarize_chimera_depth[[i]], Seg == 15)$Depth
    x$DepthCorrected = x$Depth-d+0.43 #depth on control
    x$ratioCorrected = x$DepthCorrected / x$NbChimera
    return(x)
    })
 
    
dt_cor_corrected = bind_rows(dt_cor_corrected)
dt_cor_corrected$Sample <- factor(dt_cor_corrected$Sample, levels = samples, ordered = T)
  
png(paste0(dir, "/Figures/Boxplot_Depth_NbRead_Corrected.png"),  width = 7, height = 5, units = 'in', res = 600)    
dt_cor_corrected %>% ggplot(aes(x=Sample, y=ratioCorrected, fill=Sample)) +
    geom_boxplot(fill = colors, outlier.shape = NA) +
    geom_jitter(size=0.7) +  
    scale_x_discrete(labels = names) +
    theme(legend.position = "none") +
    xlab("") +  ylab("Ratio depth/chimera") +
    geom_hline(yintercept = 0.63, linetype = "dashed")
dev.off()

names(colors) = samples
png(paste0(dir, "/Figures/Ratio.png"),  width = 7, height = 5, units = 'in', res = 600)
names(names) = samples
dt_cor_corrected %>% ggplot(aes(x=as.character(Seg), y=ratioCorrected)) +
    geom_point() +
    scale_color_manual(values=colors) +
    facet_grid(sample~., labeller = labeller(sample=names)) +
    theme(legend.position = "none") +
    xlab("Segments") +  ylab("Ratio depth/chimera") +
    geom_hline(yintercept = 0.63, linetype = "dashed")
dev.off()

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

if(help==TRUE)
{
    print("5-NbIE_alongSegments --args --config_file=\"config.Rdata\" ")
    quit("no")
}

require(data.table)
require(dplyr)
require(stringr)
library(data.table)

depth=TRUE
IPMH=TRUE
source(config_file)
source(paste0(scriptPath, "/functions.Rdata"))

print("Reading IEs..")
#Read file10 = IE
chimera_alongSegments_IE = lapply(paste0(dir, "/Chimera/", samples, "_all_chimera_alongSegments_IE.txt"), fread, header = T, sep = "\t")



#Liste Segment with integrations to give order
HIMSegments = c("S1", "S4", "S7",
    "S10", "S11", "S12", "S14", "S16", "S17", "S18", 
    "S24", "S26", "S27", "S28", 
    "S32", "S35") 

IE = lapply(chimera_alongSegments_IE, function(x){
    i = which(samples==unique(x$sample))
    x = x %>% group_by(sample, Segment) %>% summarise(NbChimera = sum(NbChimera), NbIE = n()) %>%
    mutate(., Tot_ReadsInSample = Tot_ReadsOnSesamia[i]) %>%
    mutate(., IPMH = NbIE*1000000/Tot_ReadsInSample)
    print(x)
    #For the figures, we are intereted only by segments with HIM
    y = filter(x, Segment %in% HIMSegments) 
    #Check if some HIM-segments are not integrated
    missing = HIMSegments[! HIMSegments %in% x$Segment]
    #If there are HIM-segments with no integration, add an empty line
    if (length(missing)>0) {
        y = rbind(y, data.frame(sample = unique(x$sample), Segment = missing, NbChimera =0, NbIE=0, Tot_ReadsInSample = Tot_ReadsOnSesamia[i], IPMH = 0))
        }
    return(y)
    
    #Print interesting info
    cat("\n")
    print(samples[i])
    print(y)
    })
cat("\n")


###### Do following figures with IPMS

#Prepare data for stacked barplot
data_stacked = matrix(nrow = length(samples), ncol = length(HIMSegments)+1) #Need data like that to make the stacked plot
colnames(data_stacked) = c(HIMSegments, "All segments")
rownames(data_stacked) = names

for (i in HIMSegments){
    for (j in 1:length(samples)){
        data_stacked[names[j],i] <- subset(IE[[j]], Segment==i)$IPMH 
    }#png("../All_samples/Graph_NbIE_RPM_Segments.png",  res = 600, width = 7, height = 7, units = "in")

}

data_stacked[,"All segments"] = rowSums(data_stacked, na.rm = T)

cat("\n")
print("Number of IPMS by HIM-containing segments in all samples:")
print(data_stacked)

#### Figures comparing the relative numbers of chimeric reads, with abolute values on top of the bars

#Absoulte values of Number of IE and chimeric reads:
absolute_IE_seg = bind_rows(IE) %>% group_by(Segment) %>% summarise(TotNbChimera = sum(NbChimera), TotNbIE=sum(NbIE))
absolute_IE_seg = absolute_IE_seg[match(HIMSegments, absolute_IE_seg$Segment),] #Segments in same order as on graph

png(paste0(dir, "/Figures/IPMH_Segments.png"),  res = 600, width = 7, height = 7, units = "in")
par(mar=c(5,4.5,1,1))
x = barplot(data_stacked[,1:length(HIMSegments)], names.arg = gsub('egment_','',HIMSegments), 
        cex.names = 1, cex.axis = 1.5, cex.lab = 1.5,  las = 2,
        main = "", 
        ylab = "Number of IE per million reads of host", 
        border = NA,
        col = colors,
        legend.text = T,
        ylim = c(0,ceiling(max(colSums(data_stacked[,1:length(HIMSegments)]))))) 
#Add value abosulte IE on top bars:
text(x, colSums(data_stacked)+0.07, labels =  as.character(as.matrix(absolute_IE_seg[,3])), cex = 0.9)
#Add value abosulte chimeric reads on top bars:
text(x, colSums(data_stacked)+0.02, labels =  paste0('(',as.character(as.matrix(absolute_IE_seg[,2])),')'), cex = 0.8)

dev.off()

absolute_IE_sample = bind_rows(lapply(IE, function(x){
 x = data.frame(sample = unique(x$sample), TotNbChimera = sum(x$NbChimera), TotNbIE = sum(x$NbIE))
 }))
 
#png("../All_samples/Graph_NbIE_RPM_samples.png",  res = 600, width = 7, height = 7, units = "in")
png(paste0(dir, "/Figures/IPMH_Samples.png"),  res = 600, width = 7, height = 7, units = "in")
par(mar=c(5.8,4.5,1,1))
x = barplot(data_stacked[,"All segments"], 
        cex.names = 1.2, cex.axis = 1.5, cex.lab = 1.5,  las =2,
        main = "", 
        ylab = "Number of IE per million reads of host",   
        border = NA,
        col = colors,
        legend.text = F,
        ylim = c(0,max(ceiling(data_stacked[,"All segments"]))))
text(x, data_stacked[,"All segments"]+0.11, labels =  as.character(as.matrix(absolute_IE_sample[,3])), cex = 0.9) 
text(x, data_stacked[,"All segments"]+0.03, labels =  paste0('(',as.character(as.matrix(absolute_IE_sample[,2])),')'), cex = 0.8) 
dev.off()

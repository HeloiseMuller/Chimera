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
    print("5-NbIE_alongSegments --args --config_file=\"config.Rdata\"") 
    quit("no")
}

require(data.table)
require(dplyr)
require(stringr)

depth=TRUE
IPMH=TRUE
source(config_file)
source(paste0(scriptPath, "/functions.Rdata"))

print("Reading IEs..")
#Read file10 = IE
chimera_alongSegments_IE = lapply(paste0(dir, "/Chimera/", samples, "_all_chimera_alongSegments_IE.txt"), fread, header = T, sep = "\t")
chimera_alongHIM = lapply(paste0(dir, "/Chimera/", samples, "_all_chimera_alongHIMs_IE.txt"), fread, header = T, sep = "\t")


#Liste Segments with integrations to give order
HIMSegments = c("S1", "S4", "S7",
    "S10", "S11", "S12", "S14", "S16", "S17", "S18", 
    "S24", "S26", "S27", "S28", 
    "S32", "S35") 
    
###### Calculate IPMS in HIM only ######

HIM_IE = lapply(chimera_alongHIM, function(x){
    i = which(samples==unique(x$sample))
    y = x %>% group_by(sample, Segment) %>% summarise(NbChimera = sum(NbChimera), NbIE = n()) %>%
    mutate(., Tot_ReadsInSample = Tot_ReadsOnSesamia[i]) %>%
    mutate(., IPMH = NbIE*1000000/Tot_ReadsInSample)
    #For the figures, we are intereted only by segments with HIM
    #y = filter(y, Segment %in% HIMSegments) 
    #Check if some HIM-segments are not integrated
    missing = HIMSegments[! HIMSegments %in% y$Segment]
    #If there are HIM-segments with no integration, add an empty line
    if (length(missing)>0) {
        y = rbind(y, data.frame(sample = unique(y$sample), Segment = missing, NbChimera=0, NbIE=0, Tot_ReadsInSample = Tot_ReadsOnSesamia[i], IPMH = 0))
        }
    return(y)
    })
    
###### Calculate all IPMH along segments ######

Seg_IE = lapply(chimera_alongSegments_IE, function(x){
    i = which(samples==unique(x$sample))
    y = x %>% group_by(sample, Segment) %>% summarise(NbChimera = sum(NbChimera), NbIE = n()) %>%
    mutate(., Tot_ReadsInSample = Tot_ReadsOnSesamia[i]) %>%
    mutate(., IPMH = NbIE*1000000/Tot_ReadsInSample)
    #Check if some HIM-segments are not integrated
    missing = HIMSegments[! HIMSegments %in% y$Segment]
    #If there are HIM-segments with no integration, add an empty line
    if (length(missing)>0) {
        y = rbind(y, data.frame(sample = unique(y$sample), Segment = missing, NbChimera=0, NbIE=0, Tot_ReadsInSample = Tot_ReadsOnSesamia[i], IPMH = 0))
        }
    
    #Add a colomn specifing if it is a segment with a HIM
    y = mutate(y,  HIM = ifelse(Segment %in% HIMSegments, "yes", "no"))
    
    #Add a colomn giving the percentage of IE in HIM (0% if no HIM)
    tmp = merge(y, HIM_IE[[i]], by = "Segment", suffixes = c(".Seg",".HIM")) %>%
        mutate(., PercIEinHIM = `NbIE.HIM`*100/`NbIE.Seg`) %>%
        select(., c("Segment", "PercIEinHIM"))
    y = left_join(y, tmp, by = "Segment")
    y[is.na(y$PercIEinHIM),]$PercIEinHIM = 0
    
    #Print interesting info
    cat("\n")
    print(samples[i])
    print(y, n=Inf) #Inf to print all rows
    
    return(y)
    })
cat("\n")




###### Do following figures with IPMH

#Prepare data for stacked barplot
data_stacked = matrix(nrow = length(samples), ncol = length(HIMSegments)+1) #Need data like that to make the stacked plot
colnames(data_stacked) = c(HIMSegments, "All segments")
rownames(data_stacked) = names

for (i in HIMSegments){
    for (j in 1:length(samples)){
        data_stacked[names[j],i] <- subset(HIM_IE[[j]], Segment==i)$IPMH 
    }
}

data_stacked[,"All segments"] = rowSums(data_stacked, na.rm = T)

cat("\n")
print("Number of IPMH in HIM in all samples:")
print(data_stacked)

#### Figures comparing the relative numbers of chimeric reads, with abolute values on top of the bars

#Absoulte values of Number of IE and chimeric reads:
absolute_IE_seg = bind_rows(HIM_IE) %>% group_by(Segment) %>% summarise(TotNbChimera = sum(NbChimera), TotNbIE=sum(NbIE))
absolute_IE_seg = absolute_IE_seg[match(HIMSegments, absolute_IE_seg$Segment),] #Segments in same order as on graph

svg(paste0(dir, "../FiguresArticle/Figure3.svg"), width = 10, height = 6)
namesM  = c("Caterpillar", expression(italic("N. typhae")), "Adult-1",  expression(italic("G. sparganii")), "Adult-6", expression(italic("C. phragmitella")))
colorM = c("red", "khaki1", "blue", "lightgoldenrod2","darkolivegreen4", "lightgoldenrod3")

layout(matrix(c(1,1,1,1,2,2,2), nrow=1,  byrow=TRUE))
par(mar=c(5,4.5,1,1))

x = barplot(data_stacked[,1:length(HIMSegments)], names.arg = gsub('egment_','',HIMSegments), 
        cex.names = 1, cex.axis = 1.5, cex.lab = 1.3,  las = 2,
        main = "(a)", 
        ylab = "Number of HIM-mediated IPMH", 
        border = NA,
        col = colors,
        legend.text = T, args.legend=list(legend=namesM, fill=colorM, bty = "n", x = "top", ncol=3, inset=0.1),
        ylim = c(0,4)) 
#Add value abosulte IE on top bars:
text(x, colSums(data_stacked)+0.1, labels =  as.character(as.matrix(absolute_IE_seg[,3])), cex = 0.9)
#Add value abosulte chimeric reads on top bars:
text(x, colSums(data_stacked)+0.03, labels =  paste0('(',as.character(as.matrix(absolute_IE_seg[,2])),')'), cex = 0.8)

absolute_IE_sample = bind_rows(lapply(HIM_IE, function(x){
 x = data.frame(sample = unique(x$sample), TotNbChimera = sum(x$NbChimera), TotNbIE = sum(x$NbIE))
 }))
 
par(mar=c(5.8,4.5,1,1), xpd=NA)
namesN  = c("Caterpillar", "Adult-1", "Adult-6", expression(italic("N. typhae")), expression(italic("G. sparganii")), expression(italic("C. phragmitella")))
x = barplot(data_stacked[,"All segments"], names.arg = namesN, 
        cex.names = 0.85, cex.axis = 1.5, cex.lab = 1.3, las = 2,
        main = "(b)", 
        ylab = "Number of HIM-mediated IPMH",   
        border = NA,
        col = colors,
        legend.text = F,
        ylim = c(0,max(ceiling(data_stacked[,"All segments"]))))
text(x, data_stacked[,"All segments"]+0.12, labels =  as.character(as.matrix(absolute_IE_sample[,3])), cex = 0.9) 
text(x, data_stacked[,"All segments"]+0.04, labels =  paste0('(',as.character(as.matrix(absolute_IE_sample[,2])),')'), cex = 0.8) 

segments(0.3,-0.46,3.6,-0.46, lwd=0.5)
text(1.8,-0.51, labels=expression(italic("S. nonagrioides")), cex=0.8)

dev.off()

### Script inspired from script 5- ###

library(data.table)
library(dplyr)
library(ggplot2)

depth=TRUE
IPMH=TRUE

source("Chimera/src/functions.Rdata")
source("~config.Rdata")

HIMSegments = c("S1", "S4", "S7",
    "S10", "S11", "S12", "S14", "S16", "S17", "S18", 
    "S24", "S26", "S27", "S28", 
    "S32", "S35") 
    
chimera_alongSegments_IE = lapply(paste0(dir, "/Chimera/", samples, "_all_chimera_alongSegments_IE.txt"), fread, header = T, sep = "\t")
chimera_alongHIM = lapply(paste0(dir, "/Chimera/", samples, "_all_chimera_alongHIMs_IE.txt"), fread, header = T, sep = "\t")


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
    
    #Add a colomn giving the percentage of IE in HIM (0% if no HIM) and NbIEinHIM
    tmp = merge(y, HIM_IE[[i]], by = "Segment", suffixes = c(".Seg",".HIM")) %>%
        mutate(., NbIEinHIM = `NbIE.HIM`) %>%
        mutate(PropIEinHIM = `NbIE.HIM`/`NbIE.Seg`) %>%
        select(., c("Segment", "NbIEinHIM", "PropIEinHIM"))
    y = left_join(y, tmp, by = "Segment")
    y[is.na(y$PropIEinHIM),]$PropIEinHIM = 0
    y[is.na(y$NbIEinHIM),]$NbIEinHIM = 0
    
    #Print interesting info
    cat("\n")
    print(samples[i])
    print(y, n=Inf) #Inf to print all rows
    
    return(y)
    })

dt_HIMcontaining = bind_rows(Seg_IE) %>% filter(.,  Segment %in% HIMSegments)

#Sum up value of the samples corresponding to S. nonagrioides
dt_HIMcontainingHost = filter(dt_HIMcontaining, sample == samples[1] | sample == samples[2] | sample == samples[3])  %>%
    select(., -c("Tot_ReadsInSample", "IPMH")) %>%
    group_by(Segment) %>% summarise(across(where(is.numeric), sum)) %>%
    mutate(., PropIEinHIM = NbIEinHIM/NbIE, sample = "Host") 

#Merge S. nonagrioides and the non-target species
dt_HIMcontaining2 = bind_rows(dt_HIMcontainingHost, select(filter(dt_HIMcontaining, sample %in% samples[4:6]), -c("Tot_ReadsInSample", "IPMH", "HIM")))

     
### DO stat per sample and per segment, comparing to Host
#SHould do prop.test but we need at least 5 in each category
#For host, <5 for most segments outside HIM --> do fisher test

test <- function(Host, NT){
    contingence = rbind(Host, NT) %>% mutate(NbIEhorsHIM = NbIE-NbIEinHIM) %>% select(., c('NbIEhorsHIM', 'NbIEinHIM')) 
    test = fisher.test(contingence)
    pvalue = test$p.value
    return(pvalue)
}

allpvalue = data.frame(Segment=NA, sample=NA, pvalue=NA)

i=1
for (seg in unique(dt_HIMcontaining2$Segment)) {
    HostSeg = filter(dt_HIMcontaining2, sample=="Host" & Segment==seg)
    for (NT in samples[4:6]){
        allpvalue[i,]$Segment = seg
        allpvalue[i,]$sample = NT
        allpvalue[i,]$pvalue = test(HostSeg, filter(dt_HIMcontaining2, sample == NT & Segment==seg))   
        i=i+1        
    }
}

#Add stars according to pvalue
dt_HIMcontaining2_stat = full_join(dt_HIMcontaining2, allpvalue)
annotation = mutate(dt_HIMcontaining2_stat,
    significant = case_when(
    sample=="Host" | pvalue > 0.05 ~ "",
    pvalue<0.05 & pvalue > 0.01 ~ "*",
    pvalue<0.01 & pvalue > 0.001 ~ "**",
    pvalue<0.001  ~ "***"),
    xstar = case_when(
    sample==samples[4] ~ 2,
    sample==samples[5] ~ 3,
    sample==samples[6] ~ 4
    ),
    ystar = PropIEinHIM*100+1) %>% 
    select(., c(Segment, xstar, ystar, significant, sample)) %>%
    filter(., !is.na(xstar))

#Order samples
dt_HIMcontaining2_stat$sample = factor(dt_HIMcontaining2_stat$sample, level = c("Host", samples[4:6]))

names_fig = c(expression(italic("S. nonagrioides")), expression(italic("N. typhae")), expression(italic("G. sparganii")), expression(italic("C. phragmitella")))
svg("Figure5.svg",  width = 12, height = 7)
ggplot(dt_HIMcontaining2_stat, aes(fill=sample, y=PropIEinHIM*100, x=sample)) + 
    geom_bar(position="dodge", stat="identity") +
    facet_grid(~factor(Segment, levels=HIMSegments)) +
    scale_x_discrete(limit = c("Host", samples[4:6]),
                     labels = NULL) + #names_fig
    labs( x = "", y = "% IE in HIM") +
    theme(legend.position = c(0.5,0), legend.direction="horizontal",
        axis.text.x = element_text(angle=90, size=7), #plus besoin si blanck
        axis.ticks.x = element_blank(),
        panel.background = element_rect(fill='transparent')) +
    scale_fill_manual(values= colors[c(1,4:6)], name="",  labels = names_fig) +
    geom_text(data=annotation, mapping= aes(x = xstar, y = ystar, label = significant), size=3) 
dev.off()






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
if(!exists("bed")){help = TRUE; print("bed is missing")}

if(help==TRUE)
{
    print("3-FindChimera_alongSegments.R --args --config_file=\"config.Rdata\" --bed=\"CtBV.bed\" --meta=\"Table_CtBV\"")
    quit("no")
}

require(data.table)
require(Hmisc)
require(dplyr)
require(stringr)

source(config_file)
source(paste0(scriptPath, "/functions.Rdata"))

print("Reading argument files")
bed = fread(bed)
colnames(bed) = c("contig", "begin","end","Segment")
meta = fread(meta)

print("Reading files containing all chimera...")
#Read file2 = all chimera but organized
chimera_Ctyphae_vs_Snonagrioides = lapply(paste0(dir, "/Chimera/", samples, "_chimera_", spWasp, "_vs_", spHost, "_ordered.txt"), fread, header = T, sep = "\t")

#Read file 8 = chimera that mapped on the re-integrated segments mapped after on proviral segments
SReintegrated = lapply(seq_along(samples), function(i){
    file = paste0(dir, "/Chimera/", samples[i], "_chimera_SReintegrated_vs_", spHost, ".txt")
    if(file.exists(file)){
    print(paste0("Reading file containing chimera on re-integrated segments of ", samples[i], "..."))
    fread(file) 
    }  else { 
     print(paste0("No file containing chimera on re-integrated segments for ", samples[i]))
    data.frame()}
    })

##############################################################
############## Process re-integrated segments ################
##############################################################

#ATTENTION coordinates in Cotesia begin at the beggining of the segment (instead of scf)
#Put all obs with Sesamia in colomn sp and Ctyphae in colomn sp.sz

SReintegrated_rightCoord = lapply(seq_along(samples), function(i){
    if(nrow(SReintegrated[[i]])!=0){
    print("Proccessing file containing chimera on re-integrated segments")
        x = chimera_organized(SReintegrated[[i]])
        x$Begin_segReI = as.numeric(gsub(".*[:]([^.]+)[-].*", "\\1", x$subject.s)) 
        y = x %>% mutate(chimericPoint = chimericPoint+Begin_segReI) %>%
        mutate(sStart.s = sStart.s+Begin_segReI) %>%
        mutate(sEnd.s = sEnd.s+Begin_segReI) %>%
        subset(., select=-Begin_segReI) #Don't need the colomn Begin_segReI anymore
        #Correct the names of the subject.s 
        y$subject.s =  gsub('Segment_.*::','',y$subject.s) %>% gsub(':.*','',.)
        return(y)
    } else {
        return(data.frame())
    }
})

##############################################################
############## Chimeric reads along segments #################
##############################################################
print("Merging tables containing all reads and the ones on re-integrated segments...")
chimera_alongSegments = lapply(seq_along(samples), function(i) {
    #In file containing all chimera, keep only chimera falling in segments & add a colomns specifying this segment
    x = bind_rows(apply(bed, 1, function(x) chimera_segment(file_chimera=chimera_Ctyphae_vs_Snonagrioides[[i]], bed_segment=x)))
    #If there were chimera on re-integrated segment, add a colomns specifying the name of the segment & bind table with the one above
    if(nrow(SReintegrated[[i]])!=0){
    x = x %>% rbind(., bind_rows(apply(bed, 1, function(x) chimera_segment(file_chimera=SReintegrated_rightCoord[[i]], bed_segment=x))))
    }
    #Return the table and add a colomn containing sample name
    return(x %>%  mutate(., sample = samples[i]))   
    })
    
##Extract chimera in HIM
chimera_alongHIM = lapply(seq_along(samples), function(i) {
    #Coordinates of the HIMs:
    bed = meta[, c(1,2,7,8)]
    colnames(bed) = c("Segment","contig", "begin", "end")    
    bind_rows(apply(bed, 1, function(x) chimera_segment(file_chimera=chimera_alongSegments[[i]], bed_segment=x)))
    })


########################################################
############# Calculate exact number of IE #############
########################################################
print('Calculating the number of IE on segments...')
chimera_alongSegments_IE = lapply(chimera_alongSegments, function(x){
    x = x %>% count(sample, Segment, sp.s, subject.s, chimericPoint, sp, subject, insertionCoord)
    names(x)[which(names(x) == "n")] <- "NbChimera"
    return(x)
    })

### Object containing only IE in HIM ###
chimera_alongHIM_IE = lapply(seq_along(samples), function(i) {
    #Coordinates of the HIMs:
    bed = meta[, c(1,2,7,8)]
    colnames(bed) = c("Segment","contig", "begin", "end")    
    bind_rows(apply(bed, 1, function(x) chimera_segment(file_chimera=chimera_alongSegments_IE[[i]], bed_segment=x)))
    })

########################################################
####################### Outputs #######################
########################################################
#For each sample, print the number of chimeric reads and the number of IE
resume = data.frame(NbChimera=NA, NbIE_Tot=NA)
cat("\n")

for (i in 1:length(samples)){
    resume$NbChimera = nrow(chimera_alongSegments[[i]])
    resume$NbIE_Tot = nrow(chimera_alongSegments_IE[[i]])
    resume$NbIE_HIMcontainingSegments = nrow(filter(chimera_alongSegments_IE[[i]], Segment %in%  meta[!is.na(meta$HIM_start),]$Segment))
    resume$NbIE_inHIM = nrow(chimera_alongHIM_IE[[i]])
    print(names[i])
    print(resume)
    cat("\n")
}



print("If there are chimera along segments, saving outputs...")
for (i in 1:length(samples)){
    #save only non empty outputs
    if (nrow(chimera_alongSegments[[i]])!=0){
        #Save file 9 = file of chimeric Reads along all segments I used for all my downstream analysis
        write.table(chimera_alongSegments[[i]],
        paste0(dir, "/Chimera/", samples[i], "_all_chimera_alongSegments.txt"),
        sep = "\t", row.names = F, col.names = T, quote = F)
        #Save file 10 = file of IE along all segments I used for all my downstream analysis too
        write.table(chimera_alongSegments_IE[[i]],
        paste0(dir, "/Chimera/", samples[i], "_all_chimera_alongSegments_IE.txt"),
        sep = "\t", row.names = F, col.names = T, quote = F)
        #Save file 11 = file of chimeric Reads along HIM only
        write.table(chimera_alongHIM[[i]],
        paste0(dir, "/Chimera/", samples[i], "_all_chimera_alongHIMs.txt"),
        sep = "\t", row.names = F, col.names = T, quote = F)
        #Save file 12 = file of IE along HIMs
        write.table(chimera_alongHIM_IE[[i]],
        paste0(dir, "/Chimera/", samples[i], "_all_chimera_alongHIMs_IE.txt"),
        sep = "\t", row.names = F, col.names = T, quote = F)
    } 
}



#####################################
#### Read arguments and files #######
######################################

args = commandArgs(trailingOnly=FALSE)
scriptPath <- dirname(sub("--file=","",args[grep("--file",args)]))
print(scriptPath)

args = commandArgs(trailingOnly=TRUE)
for(i in 2:length(args))
{
    eval(parse(text=sub("--", "", args[i])))
}

help = FALSE
if(!exists("config_file")){help = TRUE; print("config_file")}
if(!exists("ReInserted_Segments")){help = TRUE; print("ReInserted_Segments")}

if(help==TRUE)
{
    print("2-FindChimera_alongReIntegratedSegments.R --args --config_file=\"config.Rdata\" --ReInserted_Segments=\"ReInserted_CtBV.bed\" ")
    quit("no")
}

require(data.table)
require(Hmisc)
require(dplyr)
require(stringr)

source(config_file)
source(paste0(scriptPath, "/functions.Rdata"))

cat("\n")

#Load coordintes Re-inserted segments
print("Reading argument files")
ReInserted_Segments = fread(ReInserted_Segments)
colnames(ReInserted_Segments) = c("contig", "begin","end","Segment")

#Raw output FindChimericReads = file1
print("Reading samples")
chimera_Ctyphae_vs_Snonagrioides = lapply(paste0(dir, "/Chimera/", samples, "_chimera_", spWasp, "_vs_", spHost, ".txt"), fread, header = T, sep = "\t")

#####################################
#### Process inputs #######
######################################

#Put all obs with Sesamia in colomn sp and Ctyphae in colomn sp.s & Remove PCR duplicates if PCRdup=FALSE
print("Organizing chimera")
chimera = lapply(chimera_Ctyphae_vs_Snonagrioides, chimera_organized)


#Save file2 
for (i in 1:length(samples)){
    write.table(chimera[[i]], paste0(dir, "/Chimera/", samples[i], "_chimera_", spWasp, "_vs_", spHost, "_ordered.txt"),
    sep = "\t", row.names = F, col.names = T, quote = F)
    }

#Make bed file of all chimera (start need to be before end)
chimeraCt_bed = lapply(chimera, function(x){
    x = data.frame(scf = x$subject.s, 
        start = ifelse(x$sStart.s<x$sEnd.s, x$sStart.s, x$sEnd.s),
        end = ifelse(x$sStart.s<x$sEnd.s, x$sEnd.s, x$sStart.s))
    })

#Save file4
for (i in 1:length(samples)){
    write.table(chimeraCt_bed[[i]],
    paste0(dir, "/Chimera/", samples[i], "_chimera_", spWasp, "_vs_", spHost, "_ordered.bed_", spWasp),
    sep = "\t", row.names = F, col.names = T, quote = F)
    }
    
#Nombre reads chimerique:
print("Number of chimera in each sample, with PCR duplications:")
lapply(chimera, nrow) 

##############################################################
#### For Segments where re-insertion, sum up the reads #######
##############################################################

##### Extract the reads from the real segment and its re-insertion
#Indeed, re-integrated segments cannot form circle so integration in S. nonagrioides impossible. So if we find some chimeric reads on them, it means they are actually on the proviral segments

SReintegrated_cat = lapply(seq_along(chimera), function(i) {
    #Get data corresponding to re-inserted segments
    bind_rows(apply(ReInserted_Segments, 1, function(x) chimera_segment(file_chimera=chimera[[i]], bed_segment=x)))  
    })
    
##Save all these reads = file5
for (i in 1:length(samples)){
    if (nrow(SReintegrated_cat[[i]])!=0){
    write.table(select(SReintegrated_cat[[i]], readName),
    paste0(dir, "/Chimera/", samples[i], "_Reads_SReintegrated_cat.lst"),
    sep = "\t", row.names = F, col.names = T, quote = F)
    }
}
    
#Information on the procedure that follow
cat("\n")
samplesWithSReI = c()
for (i in 1:length(samples)){
    if (nrow(SReintegrated_cat[[i]])!=0) {
       samplesWithSReI = append(samplesWithSReI, samples[i])
    }
}

if (length(samplesWithSReI)==0){
    print("No sample have chimera on re-integrated segments. Continue with srcipt 3-")
} else {
    print("The following samples have chimera on re-integrated segments, so do not forget to run the script 2-bis before continuing the pipeline with the scrit 3-: ")
    print(samplesWithSReI)
}

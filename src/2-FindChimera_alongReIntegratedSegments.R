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
    print("2-FindChimera_alongReIntegratedSegments.R --args --config_file=\"config.Rdata\" --ReInserted_Segments=\"ReInserted_CtBV.bed\"")
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
colnames(ReInserted_Segments) = c("contig", "begin","end","segment")

#Raw output FindChimericReads = file1
print("Reading samples")
chimera_Ctyphae_vs_Snonagrioides = lapply(paste0(dir, "/Chimera/", samples, "_chimera_", spWasp, "_vs_", spHost, ".txt"), fread, header = T, sep = "\t")

#####################################
#### Process inputs #######
######################################

#Remove PCR duplicates and Put all obs with Sesamia in colomn sp and Ctyphae in colomn sp.s
print("Organizing chimera")
chimera = lapply(chimera_Ctyphae_vs_Snonagrioides, chimera_organized)

#Save file2 (like file1 but organized)
for (i in 1:length(samples)){
    write.table(chimera[[i]],
    paste0(dir, "/Chimera/", samples[i], "_chimera_", spWasp, "_vs_", spHost, "_ordered.txt"),
    sep = "\t", row.names = F, col.names = T, quote = F)
    }

#Make bed file of all chimera (start need to be before end)
chimeraCt_bed = lapply(chimera, function(x){
    x = data.frame(scf = x$subject.s, 
        start = ifelse(x$sStart.s<x$sEnd.s, x$sStart.s, x$sEnd.s),
        end = ifelse(x$sStart.s<x$sEnd.s, x$sEnd.s, x$sStart.s))
    })
    
for (i in 1:length(samples)){
    write.table(chimeraCt_bed[[i]],
    paste0(dir, "/Chimera/", samples[i], "_chimera_", spWasp, "_vs_", spHost, "_ordered.bed_", spWasp),
    sep = "\t", row.names = F, col.names = T, quote = F)
    }
    
#Nombre reads chimerique:
print("Number of chimera in each sample:")
lapply(chimera, nrow) #15200


##############################################################
#### For Segments where re-insertion, sum up the reads #######
##############################################################

##### Extract the reads from the real segment and its re-insertion
#Indeed, re-integrated segments cannot form circle so integration in S. nonagrioides impossible. So if we find some chimeric reads on them, it means they are actually on the proviral segments

SReInserted = c()
for (i in ReInserted_Segments$segment){
     SReInserted = append(gsub("_Ri.*", "", i), SReInserted)    
}   
SReInserted = unique(SReInserted)

SReintegrated_cat = lapply(seq_along(chimera), function(i) {
    #Get data corresponding to re-inserted segments
    bind_rows(apply(ReInserted_Segments, 1, function(x) chimera_segment(file_chimera=chimera[[i]], bed_segment=x)))  
    })
    
##Save all these reads = file5
for (i in 1:length(samples)){
    write.table(select(SReintegrated_cat[[i]], readName),
    paste0(dir, "/Chimera/", samples[i], "_Reads_SReintegrated_cat.lst"),
    sep = "\t", row.names = F, col.names = T, quote = F)
    }
    
##Message on the procedure that follow
cat("\n")
print("If the file '*_Reads_SReintegrated_cat' has been created, do not forget to run the script 2-bis before doing anything else")

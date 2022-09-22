require(data.table)
setwd("~/working directory")

#import tabular blastn output [outputblastn_CtBV_versus_775LepidopteranGenomes.txt or outputblastn_HdIV_versus_775LepidopteranGenomes.txt]
#these two blastn outputs were obtained using blastn with the following options: blastn -task blastn -query PDV_segments.fasta -db Lepidopteran_genomes  -outfmt '6 qaccver qlen saccver stitle slen pident length mismatch gapopen qstart qend sstart send evalue bitscore' -num_threads 20 -out outputblastn
#PDV_segments used in Heisserer et al. are those provided in Supplementary Dataset 6 and NCBI accession numbers of lepidopteran genomes are provided in Supplementary Table 6
dt=fread("outputblastn_CtBV_versus_775LepidopteranGenomes.txt")

#filter out hits shorter than 200 bp and/or with an evalue higher than 0.0001
dt2=dt[dt$V14<0.0001 & dt$V7>300]

#import coordinates of PDV HIM boundaries [can be found in Supplementary Dataset 7 of Heisserer et al.]
HIM=fread("HIMCoordinates.txt")

#filter out all blastn hits involving HdIV segments devoid of HIM
dt_seg_HIM=dt2[dt2$V1 %in% HIM$SegmentName]

#number of blastn hits involving PDV devoid of HIM
dt_seg_no_HIM=dt2[!(dt2$V1 %in% HIM$SegmentName)]
setorder(dt_seg_no_HIM, -V7)

#we determine the direction in which the target is aligned
dt_seg_HIM$strand=ifelse(dt_seg_HIM$V12 < dt_seg_HIM$V13, "+", "-")

#add column with HIM start coordinate
dt_seg_HIM$HIMstart=HIM$HIMstart[match(dt_seg_HIM$V1, HIM$SegmentName)]

#add column with HIM end coordinate
dt_seg_HIM$HIMend=HIM$HIMend[match(dt_seg_HIM$V1, HIM$SegmentName)]

#create data table containing all blastn hits covering both HIM boundaries (we expect very few if integration is mediated by HIM)
dtHIMstart_End=dt_seg_HIM[dt_seg_HIM$V10 > dt_seg_HIM$HIMend & dt_seg_HIM$V11 < dt_seg_HIM$HIMstart | dt_seg_HIM$V11 > dt_seg_HIM$HIMend & dt_seg_HIM$V10 < dt_seg_HIM$HIMstart]

#create data table containing all blastn hits for which start or end coordinate fall between HIM start and HIM end coordinates (as expected if integration occurs through breaks within HIM)
dtHIM=dt_seg_HIM[dt_seg_HIM$V10 < dt_seg_HIM$HIMend & dt_seg_HIM$V10 > dt_seg_HIM$HIMstart | dt_seg_HIM$V11 < dt_seg_HIM$HIMend & dt_seg_HIM$V11 > dt_seg_HIM$HIMstart]

#determine which coordinate in HdIV falls in HIM
dtHIM$coordHIM=ifelse(dtHIM$V10<dtHIM$HIMend & dtHIM$V10>dtHIM$HIMstart, dtHIM$V10, dtHIM$V11)

#we determine to which junction corresponds each blastn hit based on coordinates in HdIV
dtHIM$junction=ifelse(dtHIM$V11 > dtHIM$coordHIM, "J2", "J1")

#we determine if lepidopteran contig is aligned in sense or antisense orientation
dtHIM$rev=dtHIM$V12 > dtHIM$V13
dtHIM$stran=ifelse(dtHIM$rev==FALSE, "+", "-")

#new column to get species name
dtHIM$V3_sp=gsub(pattern = "TPA_asm_", replacement = "", x = dtHIM$V3)

#new dt to get species name
species=as.data.table(tstrsplit(dtHIM$V3_sp, "_", fixed=TRUE))
species$name=paste(species$V2, species$V3, sep = " ")
setnames(species, "V2", "GenusName")
setnames(species, "V3", "SpeciesName")

#add column with species name
dtHIM$species=species$name

#determine start and end coordinates of HdIV in lepidopteran contig
dtHIM$start=ifelse(dtHIM$V12 < dtHIM$V13, dtHIM$V12, dtHIM$V13)
dtHIM$end=ifelse(dtHIM$V12 < dtHIM$V13, dtHIM$V13, dtHIM$V12)

#duplicate column with lepidopteran contig name
dtHIM$V3_bis=dtHIM$V3

#remove dot from lepidopteran name to be able to work with seqtk later
dtHIM$V3_bis=gsub("\\.", "_", dtHIM$V3_bis)

#add a number to each blastn hit
dtHIM$num=c(1:nrow(dtHIM))
setorder(dtHIM, V3_bis, start)

#create BED with name of lepidopteran contig, start and end coordinate of HdIV fragment to merge them
#this is to remove redundancy caused by the fact that some segments are similar to each other and thus some lepidopteran PDV fragments align to both of them (e.g. Hd12 and Hd16 or BV10 and BV17)
toMergeBV=as.data.table(cbind(dtHIM$V3_bis, dtHIM$start, dtHIM$end, dtHIM$num, dtHIM$V15, dtHIM$stran))

write.table(toMergeBV, "toMergecat_outputblastnALL775_IV2_300bp.txt", row.names = F, col.names = F, quote = F, sep = "\t")

#Coordinates are merged using the following commands: dos2unix toMerge.txt and bedtools merge -s -c 4 -o distinct -i toMerge.txt  > toMerge.txt.merged

#import merged coordinates of blastn hits
merged=fread("toMerge.txt.merged")

#create a dt containing number of each blastn hit merged separately in two columns, in order to assign merged coordinates to initial dt of blastn hits (dtHIM)
ref=as.data.table(tstrsplit(merged$V4, ",", fixed = TRUE))

#create column containing number of blastn hit, in order to assign merged coordinates to initial dt of blastn hits (dtHIM)
merged$ref=ref$V1

#add merged start coordinate to initial dt (dtHIM)
dtHIM$FinalStart=merged$V2[match(dtHIM$num, merged$ref)]

#add merged end coordinate to initial dt (dtHIM)
dtHIM$FinalEnd=merged$V3[match(dtHIM$num, merged$ref)]

#remove all blastn hits that do not have merged blastn coordinates (those that have been merged with other blastn hits)
dtHIM=dtHIM[!is.na(dtHIM$FinalEnd)]

#import taxonomy of lepidopteran species for which a genome is available in Genbank as of December 2021, obtained from "assembly_summary_genbank.txt" file
#This file is available here: https://ftp.ncbi.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt

#taxonomy was assigned to each lepidopteran family using the taxonomizr package, which relies on taxid provided in "assembly_summary_genbank.txt"
taxo=fread("taxoLepido.txt")

#add column indicating the family of lepidopteran
dtHIM$family=taxo$V28[match(dtHIM$species, taxo$V30)]

#create dt containing PDV coordinates in lepidopteran genomes in order to extract them with seqtk
PDV_coord = as.data.table(cbind(dtHIM$V3_bis, ifelse(dtHIM$V13 > dtHIM$V12, dtHIM$V12, dtHIM$V13), ifelse(dtHIM$V13 > dtHIM$V12, dtHIM$V13, dtHIM$V12)))
write.table(PDV_coord, "PDVlongerthan300bp_coord.txt", row.names = F, col.names = F, quote = F, sep = "\t")

#add column to be able to add PDV sequence based on correspondance in seq names in dtHIM and in PDV sequences
options(scipen = 100)
dtHIM$nameCor=paste(dtHIM$V3_bis, ifelse(dtHIM$V13 > dtHIM$V12, dtHIM$V12 + 1, dtHIM$V13 + 1), sep = ":")
dtHIM$nameCor2=paste(dtHIM$nameCor, ifelse(dtHIM$V13 > dtHIM$V12, dtHIM$V13, dtHIM$V12), sep = "-")

#extract PDV sequences from lepidopteran genomes with seqtk subseq Lepidopteran_genomes PDVlongerthan300bp_coord.txt > PDVlongerThan300bp.fasta
#PDVlongerThan300bp.fasta can be made using all sequenced provided in Supplementary Tables 1 & 2 of Heisserer et al. 

#import BV sequences and add them to new column in dtHIM
library(Biostrings)
BV_HIM=readDNAStringSet("PDVlongerThan300bp.fasta")
n=as.data.table(names(BV_HIM))
n$seq=as.data.table(BV_HIM)
dtHIM$seq=n$seq[match(dtHIM$nameCor2, n$V1)]

#generate summary tables of numbers of junctions (example is given for HdIV-like fragments)
summaryJ1=dcast(data = dtHIM[dtHIM$junction=="J1"], formula = species~V1, length)
summaryJ2=dcast(data = dtHIM[dtHIM$junction=="J2"], formula = species~V1, length)
summary=merge(summaryJ1, summaryJ2, by = "species", all = T)
taxo=fread("taxoLepido.txt")
setnames(taxo, "V30", "species")
summary$family=taxo$V28[match(summary$species, taxo$species)]
setnames(summary, c("species", "H12_J1",	"H16_J1",	"Hd27_J1",	"Hd47_J1",	"H12_J2", "H16_J2", "Hd27_J2",	"Hd44.2_J2",	"Hd47_J2", "family"))
setcolorder(summary, c("species", "family", "H12_J1",	"H12_J2", "H16_J1",	"H16_J2", "Hd27_J1",	"Hd27_J2",	"Hd47_J1",	"Hd47_J2", "Hd44.2_J2"))
summary[is.na(summary)] <- 0
summary$sum=summary$H12_J1 +	summary$H12_J2 +	summary$H16_J1 + summary$H16_J2 + summary$Hd13_J1 + summary$Hd27_J1 + summary$Hd27_J2 + summary$Hd47_J1 + summary$Hd47_J2 + summary$Hd44.2_J2
setorder(summary, -sum)
summary$family=ifelse(summary$family != "0", summary$family, "Hesperiidae")
stat=fread("stat.assembly.txt")
summary$AssembSize=stat$V12[match(summary$species, stat$V13)]
summary$N50=stat$V11[match(summary$species, stat$V13)]
summary$Access=stat$V10[match(summary$species, stat$V13)]
write.table(summary, "summarycat_outputblastnALL775_IV2_300_0.0001.txt", row.names = F, col.names = T, quote = F, sep = "\t")

family=aggregate(summary[ , 3:12], by = list(summary$family), FUN = sum)
write.table(family, "summaryFamilycat_outputblastnALL775_IV2_300_0.0001.txt", row.names = F, col.names = T, quote = F, sep = "\t")












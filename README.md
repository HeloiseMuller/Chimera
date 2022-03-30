# Pipeline Chimera

## FindChimericReads.R

This script look for chimeric reads between two genomes.
It ouputs a file named  *_chimera_spWasp_vs_spHost.txt

The other scripts process chimeric reads that correspond to polyDNAviruses.

## 1-Coverage.R

This script make a plot on the sequencing depth and coverage on species1, species2, and polyDNAviruses.

## 2-FindChimera_alongReIntegratedSegments.R

This script read *_chimera_spWasp_vs_spHost.txt and outputs several files:
-  *_chimera_spWasp_vs_spHost_ordered.txt: the only difference with the input is that the table is ordered, ie that for each chimeric reads, the information of the wasp has the suffix .s and the information of host has no suffix.
-  *_chimera_spWasp_vs_spHost_ordered.bedSpWasp: bed file of the coordinates on the wasp where the chimeric reads map
-  *_Reads_SReintegrated_cat.lst: list on the reads that map on polydnavirus taht have to pe processed differently. Indeed, Re-integrated segments cannot integrate in the host genome so these chimeric reads probably correspond to the parental proviral segments.

## 2bis-FindChimera_alongReIntegratedSegments.sh

This bash script check whether chimeric reads from re-integrated segments map on the parental proviral segments. It outputs *_Reads_SReintegrated_cat_vs_CtBV.txt

## 3-FindChimera_alongSegments.R




# All the chimeric reads that fall on re-integrated segments
file5 = *_Reads_SReintegrated_cat.lst

#### Outputs specific to re-integrated segments ####

# seqtk subseq on file5 to get fasta reads 
file6 = *_Reads_SReintegrated_cat.fasta

# Blastn of file6 on proviral segments
file7 = *_Reads_SReintegrated_cat_vs_CtBV.txt

# Run the script FindChimercReads on file7 (keep the same blast on S. nonagrioides)
file8 = *_chimera_SReintegrated_vs_S.nonagrioides.txt


##### Final outputs of my script FindChimera_AlongSegments.R ####

file9 = *_chimera_alongSegments.txt
This a the file of chimeric Reads along segments I used for all my downstream analysis
Here, chimeric reads in HIMs but also outside (really anything that falls in proviral segments, without any filtre)
Here, the chimeric reads that fell in re-integrated segments are sum up with the ones in the proviral segments

file10 = *_all_chimera_alongSegments_IE.txt
This is a file cointing IE. It is like file 9 but we keep only one line per IE



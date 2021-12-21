#### Output of the script FindChimercReads ####
Here, we have all the chimera inter-species, not matter their positions
We still have PCR duplicat
file1 = *_chimera_C.typhae_vs_S.nonagrioides.txt

#### Outputs of my script FindChimera_AlongSegments ####

# Like file1 but file2 is organiszed (all Ctyphae in the same colomn)
# No PCR duplicat
file2 = *_chimera_C.typhae_vs_S.nonagrioides_ordered.txt

# Bed files of file2, with either Ctyphae coordinates or S.nonagrioides
file3 = *_chimera_C.typhae_vs_S.nonagrioides_ordered.bed_Ctyphae
file4 = *_chimera_C.typhae_vs_S.nonagrioides_ordered.bed_Snonagrioides

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



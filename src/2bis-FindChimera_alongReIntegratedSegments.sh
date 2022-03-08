samples="SR2 SR32 SRA SRB SRD SRE"
dir="/mnt/65To/Heloise/ResistantCaterpillars"

for i in $samples
do
    echo $i
	#Extract reads saved in the R script=file6:
	echo "Making fasta file"
	seqtk subseq $dir/blastn/${i}_trimmed_vs_Ctyphae.fasta $dir/Chimera/${i}_Reads_SReintegrated_cat.lst > $dir/Chimera/${i}_Reads_SReintegrated_cat.fasta
	echo "Running blastn"
	#Run blastn on sample_Reads_SReintegrated_cat.fasta, against the proviral segments=file7: 
    blastn -query $dir/Chimera/${i}_Reads_SReintegrated_cat.fasta -subject $dir/../ref/CtBV.fasta  -outfmt 6 -max_target_seqs 2 -out $dir/blastn/${i}_Reads_SReintegrated_cat_vs_CtBV.txt  
	#Don't need to blastn again on S. nonagrioides, we can keep the initial output
	#Run the Rscript FindChimericReads.R with the new output on C. typhae's side  = file8
done

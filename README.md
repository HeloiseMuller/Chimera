# Pipeline Chimera

This pipeline allows to process data the same way as in Muller et al. (2021) https://doi.org/10.1128/JVI.00684-21. It was written for the article Muller et al. (2022) https://doi.org/10.1111/mec.16685. It is composed of several scripts, located in `src/`, to run in the same order as listed here. 
The numeroted scripts process chimeric reads specifically for studies on polyDNAviruses.

In `examples/`, one can find examples of files needed to run the scripts. Examples are given for chimeric reads between CtBV (bracoviruses of _Cotesia typhae_) and the leiptoptera _Sesamia nonagrioides_.

In `scriptsArticles/`, one can find the scripts used in Muller et al. (2022) and Heisserer et al. (2022). 

## Run blastn
Blastn has to be run on the two reference genomes between which one is looking for chimeric reads. One can use the WorkflowBowBlast (https://github.com/HeloiseMuller/WorkflowBowBlast) to process the data until this step.

## FindChimericReads.R

This script look for chimeric reads between two genomes.
It reads the outputs of the two blastn and it ouputs a file named  *_chimera_spWasp_vs_spHost.txt
One can directly modify the head of the script:
```
Path and variables to set
species = c("spWasp","spHost")                      
path="/home/user/dir/"
sample="sample"
print(sample)
TrimmedReads=T #if false, all reads have the same length, if true, each read has a different length
print("Read length of the reads")
readLength=fread(paste0(path,"blastn/sample_trimmed_vs_spWasp.length"), header = F)   #number if TrimmedReads==F (e.g. 150) ; table if TrimmedReads==T (length 1st, readname 2nd)
```
Watch out to have the exact same read names in readLength and the blastn outputs.

The script will automatically read the 2 outputs of blastn located in the subdirectory blastn of `path`:
`path/blastn/sample_trimmed_vs_spWasp.txt` and `path/blastn/sample_trimmed_vs_spWasp_vs_spHost.txt`
If the locations or the file names are different, one can directly give the name of the files in the script with the variables `blastOnGenome1` and `blastOnGenome2`.

The script will output the chimeric reads found between spWasp and spHost in the subdirectory Chimera of `path`, under the name sample_chimera_spWasp_vs_spHost.txt. One can change this automatic name by changing the file name at the very bottom of the script in write.table().

Run the script with `Rscript FindChimericReads.R`

## 1-Coverage.R

This script makes a plot of the sequencing depth and coverage on species1, species2, and polyDNAviruses.

## 2-FindChimera_alongReIntegratedSegments.R

This script reads the file *_chimera_spWasp_vs_spHost.txt and outputs several files:
-  *_chimera_spWasp_vs_spHost_ordered.txt: the only difference with the input is that the table is ordered, i.e. that for each chimeric read, the information on the wasp has the suffix .s and the information on the host has no suffix.
-  *_chimera_spWasp_vs_spHost_ordered.bedSpWasp: bed file of the coordinates on the wasp where the chimeric reads map
-  *_Reads_SReintegrated_cat.lst: list of the reads that map on polydnavirus that have to be processed differently. Indeed, Re-integrated segments cannot integrate in the host genome so these chimeric reads probably correspond to the parental proviral segments.

Run the scrip with:
```
Rscript 2-FindChimera_alongReIntegratedSegments.R --args --config_file=\"config.Rdata\" --ReInserted_Segments=\"ReInserted_CtBV.bed\"`
```

## 2bis-FindChimera_alongReIntegratedSegments.sh

This bash script checks whether chimeric reads from re-integrated segments map on the parental proviral segments. It outputs *_Reads_SReintegrated_cat_vs_CtBV.txt. One can use this output to run FindChimericReads.R again, to validate these chimeric reads. For this, replace the output of the blastn against spWasp by *_Reads_SReintegrated_cat_vs_CtBV.txt, keep the output of the blastn againt spHost, and change the name of the output at the very bottom of the scipt FindChimericReads.R.

## 3-FindChimera_alongSegments.R

This script identifies chimeric reads and estimate the number of integration events.

Inputs: 
- a configuration file
- the bed file of proviral segments
- metadata about the proviral segments, with their positions, the presence of HIM, and their orientations
- *_Reads_SReintegrated_cat_vs_CtBV.txt. The script automaticaly looks whether this file exists in the directory. If it does, it sums up their number of chimeric reads. If the file does not exist, it will simply pass this step. 

Outputs:
- *_all_chimera_alongSegments.txt: file that contains all chimeric reads between the host and the proviral segments from the bed file.
- *_all_chimera_alongSegments_noPCRdup.txt: same as above, but without the PCR duplicates.
-  *_all_chimera_alongHIMs.txt: same as above (no PCR duplicates either), but the file contains only chimeric reads that map in HIMs.
- *_all_chimera_alongSegments_IE.txt: file that contais one line per integration events (IE). It specifies how many chimeric reads cover each IE.
- *_all_chimera_alongHIMs_IE.txt: as above, but the file contains only IE mediated by HIMs.

Run the scrip with:
```
Rscript 3-FindChimera_alongSegments.R --args --config_file=\"config.Rdata\" --bed=\"CtBV.bed\" --meta=\"Table_CtBV.txt\"
```

## 4-FigureChimericReads_alongSegments.R

This R script outputs a pdf that shows where chimeric reads map for all segments with more than 2 chimeric reads. On the right, it plots a zoom on HIM, if the segment has HIM.

One has to run the script a first time to visualy set some graphical parametes, before running it a second time.
```
Rscript 4-FigureChimericReads_alongSegments.R --args --config_file=\"config.Rdata\" --Annotation_segments=\"Table_CtBV.txt\" --initialize=1
```
Depending on the height of the bars, fill the file `ParametersFigure`. The first colomn is the name of the segment, the second set the height of the y axis of the left pannel (min, max, tick every) and the third set it for the right panels. Set NA in the third colomn for segments which are devoid of HIM.

Then, the user can run the script a second time with:
```
Rscript 4-FigureChimericReads_alongSegments.R --args --config_file=\"config.Rdata\" --Annotation_segments=\"Table_CtBV.txt\" --initialize=0 --ParametersFigure=\"ParametersFigure.txt\" --init_bas=400 --init_ytext=200 --init_size_seg=650 --init_yleg=850
```
The user might want to modify the four last parameters to get a cleaner graph. \
 `--init_bas`: set the width of the rectangle representing the proviral segment \
 `--init_ytext`: set the y position of the name of the proviral segment \
 `--init_size_seg`: set the size of the segments that shows the begining and end of each segment \
 `--init_yleg`: set the y position of the number under the segments that shows the begining and end of each segment 
 
 ## 5-NbIE_alongSegments.R
 
 This script process the number of IE along segments, without PCR duplications.
 
  ```
 5-NbIE_alongSegments --args --config_file=\"config.Rdata\"
 ```
 It reads two files:*_all_chimera_alongSegments_IE.txt and *_all_chimera_alongHIMs_IE.txt
 
 It outputs 2 figures: both show the number of HIM-mediated IPMH (number of Integration event Per Million reads of Host sequenced), one for each sample, and one for each proviral segment.
 It also outputs summaized tables giving the number of IPMH for each sample and segment, and also the percentage of IE that falls in HIM.
 
 ## 6-NbIEvsDepth.R
 
This script compares the number of chimeric reads to the sequencing depth for each segment.
Since PCR duplicates were not filtered out to calculate the sequencing depth, we keep them too in the number of chimeric reads.
```
6-NbIEvsDepth.R --args --config_file=\"config.Rdata\" --depth_dir=\"directory/\" --meta=\"Table_CtBV.txt\"
```
`depth_dir` is the direcotry where is located one file per sample, each containing the sequencing depth on each segment. These files are composed of 5 colomns: the first 3 colomns correspond to a bed file for each segment, the 4th is the segment name, and the 5th is the sequencing depth on this segment. This file can be generated with CovWindows (https://github.com/HeloiseMuller/CovWindows).

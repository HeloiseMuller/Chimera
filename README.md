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
-  *_Reads_SReintegrated_cat.lst: list of the reads that map on polydnavirus taht have to pe processed differently. Indeed, Re-integrated segments cannot integrate in the host genome so these chimeric reads probably correspond to the parental proviral segments.

## 2bis-FindChimera_alongReIntegratedSegments.sh

This bash script check whether chimeric reads from re-integrated segments map on the parental proviral segments. It outputs *_Reads_SReintegrated_cat_vs_CtBV.txt

## 3-FindChimera_alongSegments.R

This script identified chimeric reads and estimate the number of integration events.

Inputs: 
- a configuration file
- the bed file of proviral segments
- metadata about the proviral segments, with their positions, the presence of HIM, and their orientation
- *_Reads_SReintegrated_cat_vs_CtBV.txt if the file exist in the directory, it sum up their number of chimeric reads

Outputs:
- *_all_chimera_alongSegments.txt: file that contains all chimeric reads between the host and the proviral segments from the bed fime
- *_all_chimera_alongSegments_noPCRdup.txt: same as above, but without the PCR duplications
-  *_all_chimera_alongHIMs.txt: same as above (no PCR dup either), but the file contains only chimeric reads that map in HIMs.
- *_all_chimera_alongSegments_IE.txt: file that contais one line per integration events (IE). Specified how many chimeric reads cover each IE.
- *_all_chimera_alongHIMs_IE.txt: as above, but the file contains only IE mediated by HIMs

Run the scrip with:
```
Rscript 3-FindChimera_alongSegments.R --args --config_file=\"config.Rdata\" --bed=\"CtBV.bed\" --meta=\"Table_CtBV.txt\"
```

## 4-FigureChimericReads_alongSegments.R

This R script output a pdf that shows where chimeric reads map for all segments with more than 2 chimeric reads. On the right, it plots a zoom on HIM, if the segment has HIM

One has to run the script a first time to visualy set some graphical parametes, before running it a second time.
```
Rscript 4-FigureChimericReads_alongSegments.R --args --config_file=\"config.Rdata\" --Annotation_segments=\"Table_CtBV.txt\" --initialize=1
```
Depending on the height of the bars, fill the file `ParametersFigure`. The first colomn is the name of the segment, the second set the left pannels and the third set the right panels. One has to give 3 numbers: the begining of the x axis, the end, show value on the axis.

Then, the user can run the script a second time with:
```
Rscript 4-FigureChimericReads_alongSegments.R --args --config_file=\"config.Rdata\" --Annotation_segments=\"Table_CtBV.txt\" --initialize=0 --ParametersFigure=\"ParametersFigur.txt\" --init_bas=400 --init_ytext=200 --init_size_seg=650 --init_yleg=850
```
The user might want to modify the four last parameters to get a cleaner graph.
`--init_bas`: set the width of the rectangle representing the proviral segment
 `--init_ytext`: set the y position of the name of the proviral segment
 `--init_size_seg`: set the size of the segments that shows the begining and end of each segment 
 `--init_yleg`: set the y position of the number under the segments that shows the begining and end of each segment 
 
 ## 5-NbIE_alongSegments.R
 
 This script process the number of IE along segments, without PCR duplications.
 
  ```
 5-NbIE_alongSegments --args --config_file=\"config.Rdata\"
 ```
 It reads two files:*_all_chimera_alongSegments_IE.txt and *_all_chimera_alongHIMs_IE.txt
 
 It outputs 2 figures: both show the the number of IMPH (number of integration event per million reads of host sequenced), one for each sample, and one for each proviral segment.
 It also outputs summaized tables giving the number of IPMH for each sample and segment, and also the percentage of IE that falls in HIM
 
 ## 6-NbIEvsDepth.R
 

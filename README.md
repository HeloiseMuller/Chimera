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

Inputs: 
- a configuration file
- the bed file of proviral segments
- metadata about the proviral segments, with their positions, the presence of HIM, and their orientation
- *_Reads_SReintegrated_cat_vs_CtBV.txt if the file exist in the directory

Outputs:
- *_all_chimera_alongSegments.txt: file that contains all chimeric reads between the host and the proviral segments from the bed fime
- *_all_chimera_alongSegments_noPCRdup.txt: same as above, but without the PCR duplications
-  *_all_chimera_alongHIMs.txt: same as above (no PCR dup either), but the file contains only chimeric reads that map in HIMs.
- *_all_chimera_alongSegments_IE.txt: file that contais one line per integration events (IE). Specified how many chimeric reads cover each IE.
- *_all_chimera_alongHIMs_IE.txt: as above, but the file contains only IE mediated by HIMs

Run the scrip with
```
Rscript 3-FindChimera_alongSegments.R --args --config_file=\"config.Rdata\" --bed=\"CtBV.bed\" --meta=\"Table_CtBV.txt\"
```

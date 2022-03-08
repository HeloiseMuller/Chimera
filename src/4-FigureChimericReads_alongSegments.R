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
if(!exists("init_bas")){help = TRUE; print("init_bas")}
if(!exists("init_ytext")){help = TRUE; print("init_ytext")}
if(!exists("init_size_seg")){help = TRUE; print("init_size_seg")}
if(!exists("init_yleg")){help = TRUE; print("init_yleg")}

if(help==TRUE)
{
    print("4-FigureChimericReads_alongSegments.R --args --config_file=\"config.Rdata\" --Annotation_segments=\"Table_CtBV.txt\" --init_bas=400 --init_ytext=200 --init_size_seg=650 --init_yleg=850")
    quit("no")
}
  
require(data.table)
require(Hmisc)
require(dplyr)
require(stringr)
library(data.table)

source(config_file)
source(paste0(scriptPath, "/functions.Rdata"))

print("Read the file cointing all chimera along segments, without PCR duplications")
chimera_C.typhae_vs_S.nonagrioides = bind_rows(lapply(paste0(dir,"/Chimera/", samples, "_all_chimera_alongSegments_noPCRdup.txt"), fread))

print("Read annotation table")
annotation = fread(Annotation_segments)

print("Number of chimera on each segment")
table(chimera_C.typhae_vs_S.nonagrioides$Segment)
annotation_segHIM = annotation[annotation$Segment %in% unique(chimera_C.typhae_vs_S.nonagrioides$Segment),]

##############################################################
######################## FUNCTION ###########################
##############################################################

#Function to plot chimeric reads on segments
plot_chim_segments = function(chimera_segment, annotation, PosJ1, PosJ2, quality, yaxis, yaxisHIM){

#quality = "forestgreen" or "darkorange" or "red"
#yaxis: set y axis of the left plot c(0, max, interval)
#ytextHIM: set y axis of the right plot c(0, max, interval)
#PosJ1 and PosJ2: Positions I took into account to count IE (parameter not needed in this version of the script)

#Parameters in the file annotation
#Contig: contig where the segment is located
#Start: Begining of the segment
#End: End of the segment
#HIM_start: Beginning of the HIM of the segment
#HIM_end: End of the HIM of the segment
#Orientation : si -1, séquence 3'-> 5' donc il faut inverser

#Parameters to set graphical constants 
breaks1=1000
breaks2=100
cex_axis=0.5
lwd_axis=0.5
hadj_axis=0.5
tck_axis=-0.1
tick_ratio=0.3
lwd_minor=0.5
cex_text=0.8
lwd_seg=0.5


#Set all these parameters by hand for the first segment
    #--> use ratios for the othets
bas = yaxis[2]/1500*init_bas         #width of the segment on the first plot (here, bas=400 for the first segment)
ytext = yaxis[2]/1500*init_ytext       #position to write the name of the segment (e.g S1)
size_seg = yaxis[2]/1500*init_size_seg   #taille du segment indiquand les positions de début et fin dans le 1e plot
yleg = yaxis[2]/1500*850 *init_yleg       #position en y de la légende du 1e plot (here, bas=850 for the first segment)

#Don't need to set those ones by hand
size_segHIM = yaxisHIM[2] * size_seg / yaxis[2]     #taille du segment indiquand les positions de début et fin dans le 2e plot
ytextHIM = yaxisHIM[2] * yleg / yaxis[2]            #position en y de la légende du 2e plot
basHIM = yaxisHIM[2] * bas / yaxis[2]               #width of the segment on the second plot


annotation_seg = annotation[which(annotation$Segment == unique(chimera_segment$segment)),]

if (annotation_seg$Orientation==-1){
  xlim_C=c(annotation_seg$End*annotation_seg$Orientation, (annotation_seg$End-length_max)*annotation_seg$Orientation)
  text_C_debut=annotation_seg$End-annotation_seg$Start
  text_C_fin="1"
  xlim_HIM = c(annotation_seg$HIM_end*annotation_seg$Orientation, (annotation_seg$HIM_end-max_HIM)*annotation_seg$Orientation)
  text_HIM_debut=annotation_seg$End-annotation_seg$HIM_start 
  text_HIM_fin=annotation_seg$End-annotation_seg$HIM_end
} else {
  xlim_C = c(annotation_seg$Start, annotation_seg$Start+length_max)
  text_C_debut="1"
  text_C_fin=annotation_seg$End-annotation_seg$Start
  xlim_HIM = c(annotation_seg$HIM_start,annotation_seg$HIM_start+max_HIM)
  text_HIM_debut=annotation_seg$HIM_start-annotation_seg$Start
  text_HIM_fin=annotation_seg$HIM_end-annotation_seg$Start
}
hist(chimera_segment$chimericPoint*annotation_seg$Orientation, breaks = breaks1, xlab = "", ylab = "",main = "",
     xaxt='n', yaxt='n',
     ylim=c(yaxis[1], yaxis[2]),
     xlim = xlim_C)
axis(side=2, at=seq(from=yaxis[1], to=yaxis[2], by=yaxis[3]),labels=seq(from=yaxis[1], to=yaxis[2], by=yaxis[3]),
     lty = 1, col='black', las=2, cex.axis=cex_axis, lwd = lwd_axis, hadj = hadj_axis, tck=tck_axis)
minor.tick(ny=2, nx=1, tick.ratio = tick_ratio, y.args=list(lwd=lwd_minor))
rect(annotation_seg$Start*annotation_seg$Orientation, -bas, annotation_seg$End*annotation_seg$Orientation, 0, col = quality, border = NA)
text(x=(annotation_seg$Start+annotation_seg$End)/(2*annotation_seg$Orientation),y=-ytext, deparse(substitute(chimera_segment)), cex = cex_text)
rect((mean(c(annotation_seg$HIM_start, annotation_seg$HIM_end))-150)*annotation_seg$Orientation, 0, (mean(c(annotation_seg$HIM_start, annotation_seg$HIM_end))+150)*annotation_seg$Orientation, -bas, col = "white", border = NA)
segments(annotation_seg$End*annotation_seg$Orientation, -bas, annotation_seg$End*annotation_seg$Orientation, -size_seg, lwd = lwd_seg)
segments(annotation_seg$Start*annotation_seg$Orientation, -bas, annotation_seg$Start*annotation_seg$Orientation, -size_seg, lwd = lwd_seg)
text(x=annotation_seg$Start*annotation_seg$Orientation, y=-yleg, text_C_debut, cex = cex_text)
text(x=annotation_seg$End*annotation_seg$Orientation, y=-yleg, text_C_fin, cex = cex_text)

#Did i took orientation into account there?
chimera_aroundHIM = chimera_segment[which(chimera_segment$chimericPoint>=annotation_seg$HIM_start & chimera_segment$chimericPoint<=annotation_seg$HIM_end),]
hist(chimera_aroundHIM$chimericPoint*annotation_seg$Orientation, breaks = breaks2, xaxt='n', yaxt='n',
    xlab = "", ylab = "", main = "",
    col = "black",
    border = F,
    xlim = xlim_HIM, 
    ylim = c(0,yaxisHIM[2]))
axis(side=2, at=seq(from=yaxisHIM[1], to=yaxisHIM[2], by=yaxisHIM[3]),labels=seq(from=yaxisHIM[1], to=yaxisHIM[2], by=yaxisHIM[3]), 
     lty = 1, col='black', las=2, cex.axis=cex_axis, lwd = lwd_axis, hadj = hadj_axis, tck=tck_axis)
minor.tick(ny=2, nx=1, tick.ratio = tick_ratio, y.args=list(lwd=lwd_minor))
segments(annotation_seg$HIM_start*annotation_seg$Orientation, 0, annotation_seg$HIM_start*annotation_seg$Orientation, -size_segHIM, lwd = lwd_seg)
text(x=annotation_seg$HIM_start*annotation_seg$Orientation, y=-ytextHIM, text_HIM_debut, cex = cex_text)
segments(annotation_seg$HIM_end*annotation_seg$Orientation, 0, annotation_seg$HIM_end*annotation_seg$Orientation, -size_segHIM, lwd = lwd_seg)
text(x=annotation_seg$HIM_end*annotation_seg$Orientation, y=-ytextHIM, text_HIM_fin, cex = cex_text)
rect(annotation_seg$HIM_start*annotation_seg$Orientation, -basHIM, annotation_seg$HIM_end*annotation_seg$Orientation, 0, col = "white", border = quality)

}

##############################################################
############## Chimeric reads along segments #################
##############################################################

#Parameters to have proportional figures
length_max = max(annotation_segHIM$End - annotation_segHIM$Start) #S28 --> permet d'avoir la même echelle en x sur le graph colonne gauche
max_HIM = max(annotation_segHIM$HIM_end - annotation_segHIM$HIM_start, na.rm=T)

pdf(paste0(dir,"/Figures/ChimericReads_alongSegments_oriented.pdf"))
par(mfrow = c(8,2) , mar = c(4.1,4.1,0,3), oma=c(0.2,0.2,2,2), xpd=T)
#mar: marges interieur (pour chaque plot)
#oam: marges exérieurs
#xpd=T: considéré tout ce qui est autout figure (axe, légende, ect) comme étant à lintérieur de la région du plot

#PosJ1 and PosJ2 are useless. It's the positions I took into account for IE

###########
##Circle 1
S1 = chimera_C.typhae_vs_S.nonagrioides[chimera_C.typhae_vs_S.nonagrioides$segment=="Segment_1",]
plot_chim_segments(chimera_segment=S1, annotation, PosJ1=10780347, PosJ2=10780388, quality = "forestgreen", yaxis=c(0,200,100),  yaxisHIM=c(0,100,50))

###########
##Circle 4
S4 = chimera_C.typhae_vs_S.nonagrioides[chimera_C.typhae_vs_S.nonagrioides$segment=="Segment_4",]
plot_chim_segments(chimera_segment=S4, annotation, PosJ1=998158, PosJ2=998200, quality = "forestgreen", yaxis=c(0,20,10), yaxisHIM=c(0,20,10))

###########
##Circle 7
S7 = chimera_C.typhae_vs_S.nonagrioides[chimera_C.typhae_vs_S.nonagrioides$segment=="Segment_7",]
plot_chim_segments(chimera_segment=S7, annotation, PosJ1=355800, PosJ2=355852, quality = "forestgreen", yaxis=c(0,50,25), yaxisHIM=c(0,20,10))

########### 
##Circle 10 
S10 = chimera_C.typhae_vs_S.nonagrioides[chimera_C.typhae_vs_S.nonagrioides$segment=="Segment_10",]
plot_chim_segments(chimera_segment=S10, annotation, PosJ1=1075760, PosJ2=1075802, quality = "forestgreen", yaxis=c(0,10,5), yaxisHIM=c(0,10,5))


###########
##Circle 11
S11 = chimera_C.typhae_vs_S.nonagrioides[chimera_C.typhae_vs_S.nonagrioides$segment=="Segment_11",]
plot_chim_segments(chimera_segment=S11,annotation, PosJ1=1347890, PosJ2=1347940, quality = "forestgreen", yaxis=c(0,20,10),  yaxisHIM=c(0,10,5))


###########
##Circle 12
S12 = chimera_C.typhae_vs_S.nonagrioides[chimera_C.typhae_vs_S.nonagrioides$segment=="Segment_12",]
plot_chim_segments(chimera_segment=S12, annotation, PosJ1=1359908, PosJ2=1359969, quality = "forestgreen", yaxis=c(0,20,10), yaxisHIM=c(0,20,10))


###########
##Circle 14  #Oups! Probleme avec position Start ou End. Check annotation
S14 = chimera_C.typhae_vs_S.nonagrioides[chimera_C.typhae_vs_S.nonagrioides$segment=="Segment_14",]
plot_chim_segments(chimera_segment=S14, annotation, PosJ1=1317495, PosJ2=NA,  quality = "red", yaxis=c(0,20,10), yaxisHIM=c(0,20,10))

########### 
##Circle 16 
S16 = chimera_C.typhae_vs_S.nonagrioides[chimera_C.typhae_vs_S.nonagrioides$segment=="Segment_16",]
plot_chim_segments(chimera_segment=S16, annotation, PosJ1=132515, PosJ2=132555,  quality = "forestgreen", yaxis=c(0,20,10), yaxisHIM=c(0,20,10))

###########
##Circle 17
S17 = chimera_C.typhae_vs_S.nonagrioides[chimera_C.typhae_vs_S.nonagrioides$segment=="Segment_17",]
plot_chim_segments(chimera_segment=S17, annotation, PosJ1=1076654, PosJ2=1076696, quality = "forestgreen", yaxis=c(0,20,10), yaxisHIM=c(0,20,10))


###########  
##Circle 18
S18 = chimera_C.typhae_vs_S.nonagrioides[chimera_C.typhae_vs_S.nonagrioides$segment=="Segment_18",]
plot_chim_segments(chimera_segment=S18, annotation, PosJ1=122785, PosJ2=122826, quality = "forestgreen", yaxis=c(0,20,10), yaxisHIM=c(0,20,10))

###########  
##Circle 24
S24 = chimera_C.typhae_vs_S.nonagrioides[chimera_C.typhae_vs_S.nonagrioides$segment=="Segment_24",]
plot_chim_segments(chimera_segment=S24 ,annotation, PosJ1=85737, PosJ2=85789, quality = "red", yaxis=c(0,20,10), yaxisHIM=c(0,20,10))

###########
##Circle 26 
S26 = chimera_C.typhae_vs_S.nonagrioides[chimera_C.typhae_vs_S.nonagrioides$segment=="Segment_26",]
plot_chim_segments(chimera_segment=S26, annotation, PosJ1=11087799, PosJ2=NA, quality = "forestgreen", yaxis=c(0,50,25), yaxisHIM=c(0,50,25))

###########
##Circle 27
S27 = chimera_C.typhae_vs_S.nonagrioides[chimera_C.typhae_vs_S.nonagrioides$segment=="Segment_27",]
plot_chim_segments(chimera_segment=S27, annotation, PosJ1=30034, PosJ2=30086, quality = "darkorange", yaxis=c(0,30,15),  yaxisHIM=c(0,30,15))


###########
##Circle 28 
S28 = chimera_C.typhae_vs_S.nonagrioides[chimera_C.typhae_vs_S.nonagrioides$segment=="Segment_28",]
plot_chim_segments(chimera_segment=S28, annotation, PosJ1=41330, PosJ2=41403, quality = "forestgreen", yaxis=c(0,30,15),  yaxisHIM=c(0,30,15))


###########
##Circle 32
S32 = chimera_C.typhae_vs_S.nonagrioides[chimera_C.typhae_vs_S.nonagrioides$segment=="Segment_32",]
plot_chim_segments(chimera_segment=S32, annotation, PosJ1=63682, PosJ2=63751, quality = "forestgreen", yaxis=c(0,30,15), yaxisHIM=c(0,30,15))


###########
##Circle 35
S35 = chimera_C.typhae_vs_S.nonagrioides[chimera_C.typhae_vs_S.nonagrioides$segment=="Segment_35",]
plot_chim_segments(chimera_segment=S35, annotation, PosJ1=90753, PosJ2=90806, quality = "forestgreen", yaxis=c(0,20,10), yaxisHIM=c(0,20,10))

dev.off()



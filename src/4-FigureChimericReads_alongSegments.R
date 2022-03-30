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
if(!exists("initialize") & !exists("ParametersFigure")){help = TRUE; print("Either ParametersFigure.txt is missing, or you need to set initialize to true")}
if(initialize==1 & exists("ParametersFigure")){help = TRUE; print("initialize and ParametersFigure are not compatible")}
if(initialize!=1 & !exists("ParametersFigure")){help = TRUE; print("ParametersFigure missing")}
if(initialize!=1 & !exists("init_bas")){help = TRUE; print("init_bas missing")}
if(initialize!=1 & !exists("init_ytext")){help = TRUE; print("init_ytext missing")}
if(initialize!=1 & !exists("init_size_seg")){help = TRUE; print("init_size_seg missing")}
if(initialize!=1 & !exists("init_yleg")){help = TRUE; print("init_yleg missing")}

if(help==TRUE)
{
    print("4-FigureChimericReads_alongSegments.R --args --config_file=\"config.Rdata\" --Annotation_segments=\"Table_CtBV.txt\" --initialize=0 --ParametersFigure=\"ParametersFigureAlongSegment.txt\" --init_bas=400 --init_ytext=200 --init_size_seg=650 --init_yleg=850")
    quit("yes")
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

if(exists("ParametersFigure")){
    print("Read parameters for the figure")
    ParametersFigure = fread(ParametersFigure, quote=F)
    }

print("Number of chimera on each segment")
table(chimera_C.typhae_vs_S.nonagrioides$Segment)
annotation_segHIM = filter(annotation, !is.na(HIM_start))

##############################################################
######################## FUNCTION ###########################
##############################################################

#Function to plot chimeric reads on segments
plot_chim_segments = function(chimera_segment, annotation, PosJ1, PosJ2, quality, yaxis, yaxisHIM, HIM){

    #quality = "forestgreen" or "darkorange" or "red"
    #yaxis: set y axis of the left plot c(0, max, interval)
    #ytextHIM: set y axis of the right plot c(0, max, interval)

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
    
    #Segment to draw:
    annotation_seg = annotation[which(annotation$Segment == unique(chimera_segment$Segment)),]

    #To initialize, let plot set parameters automatically. The user can then set the parameters to get a clean figure
    if(initialize==1){
        hist(chimera_segment$chimericPoint*annotation_seg$Orientation, breaks = breaks1, xlab = annotation_seg$Segment, ylab = "", main = "",   xaxt='n')
        #If there are HIMs, make a zoom on this area
        if(HIM==TRUE){
            #Did i took orientation into account there?
            chimera_aroundHIM = chimera_segment[which(chimera_segment$chimericPoint>=annotation_seg$HIM_start & chimera_segment$chimericPoint<=annotation_seg$HIM_end),]
            hist(chimera_aroundHIM$chimericPoint*annotation_seg$Orientation, breaks = breaks2, xaxt='n', xlab="", ylab = "", main = "", col = "black", border = F)
        } else {
            plot(NULL, xlim=c(0,1), ylim=c(0,1), ylab="", xlab="", xaxt='n', yaxt='n')
        }

    #once the user chose saw the size of the bars with initialize, they can set yaxis and yaxisHIM
    } else {
        
        #Set all the following parameters as ratio of yaxis.The user might want to modify
        bas = yaxis[2]/1500*init_bas         #width of the segment on the first plot (here, bas=400 for the first segment)
        ytext = yaxis[2]/1500*init_ytext       #position to write the name of the segment (e.g S1)
        size_seg = yaxis[2]/1500*init_size_seg   #taille du segment indicant les positions de début et fin dans le 1e plot
        yleg = yaxis[2]/1500*init_yleg       #position en y de la légende du 1e plot (here, bas=850 for the first segment)

        #Don't need to set those ones by hand
        size_segHIM = yaxisHIM[2] * size_seg / yaxis[2]     #taille du segment indiquand les positions de début et fin dans le 2e plot
        ytextHIM = yaxisHIM[2] * yleg / yaxis[2]            #position en y de la légende du 2e plot
        basHIM = yaxisHIM[2] * bas / yaxis[2]               #width of the segment on the second plot

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
        hist(chimera_segment$chimericPoint*annotation_seg$Orientation, breaks = breaks1, xlab = "", ylab = "",main = "", xaxt='n', yaxt='n', ylim=c(yaxis[1], yaxis[2]), xlim = xlim_C)
        axis(side=2, at=seq(from=yaxis[1], to=yaxis[2], by=yaxis[3]),labels=seq(from=yaxis[1], to=yaxis[2],
        by=yaxis[3]), lty = 1, col='black', las=2, cex.axis=cex_axis, lwd = lwd_axis, hadj = hadj_axis, tck=tck_axis)
        minor.tick(ny=2, nx=1, tick.ratio = tick_ratio, y.args=list(lwd=lwd_minor))
        rect(annotation_seg$Start*annotation_seg$Orientation, -bas, annotation_seg$End*annotation_seg$Orientation, 0, col = quality, border = NA)
        text(x=(annotation_seg$Start+annotation_seg$End)/(2*annotation_seg$Orientation),y=-ytext, annotation_seg$Segment, cex = cex_text)
        rect((mean(c(annotation_seg$HIM_start, annotation_seg$HIM_end))-150)*annotation_seg$Orientation, 0, (mean(c(annotation_seg$HIM_start, annotation_seg$HIM_end))+150)*annotation_seg$Orientation, -bas, col = "white", border = NA)
        segments(annotation_seg$End*annotation_seg$Orientation, -bas, annotation_seg$End*annotation_seg$Orientation, -size_seg, lwd = lwd_seg)
        segments(annotation_seg$Start*annotation_seg$Orientation, -bas, annotation_seg$Start*annotation_seg$Orientation, -size_seg, lwd = lwd_seg)
        text(x=annotation_seg$Start*annotation_seg$Orientation, y=-yleg, text_C_debut, cex = cex_text)
        text(x=annotation_seg$End*annotation_seg$Orientation, y=-yleg, text_C_fin, cex = cex_text)
        
        #If there are HIMs, make a zoom on this area
        if(HIM==TRUE){
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
        #Otherwise, leave an empty pannel
        } else {
            plot(NULL, xlim=c(0,1), ylim=c(0,1), ylab="", xlab="", xaxt='n', yaxt='n')
        }
    }
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


###########

if(initialize==1){
    for (i in annotation$Segment){
        print(i)
        seg = chimera_C.typhae_vs_S.nonagrioides[chimera_C.typhae_vs_S.nonagrioides$Segment==i,]
        if(nrow(seg)<2){
            print(paste0("Number of chimeric reads <2: ", i))
        } else if(i %in% annotation_segHIM$Segment){
            plot_chim_segments(chimera_segment=seg, annotation, quality = "forestgreen", yaxis=NA,  yaxisHIM=NA, HIM=TRUE)
        } else {
            plot_chim_segments(chimera_segment=seg, annotation, quality = "forestgreen", yaxis=NA,  yaxisHIM=NA, HIM=FALSE)
        }
    }
} else {
    for (i in annotation$Segment){
        print(i)
        seg = chimera_C.typhae_vs_S.nonagrioides[chimera_C.typhae_vs_S.nonagrioides$Segment==i,]
         #Get yaxis of this segment
        yaxis=filter(ParametersFigure,Segment==i)$yaxis
        yaxis = as.integer(unlist(strsplit(yaxis, ",")))
        if(nrow(seg)<2){
            print(paste0("Number of chimeric reads <2 for segment: ", i))
        } else if(i %in% annotation_segHIM$Segment) {
            #Get yaxisHIM of this semgent
            yaxisHIM=filter(ParametersFigure,Segment==i)$yaxisHIM
            yaxisHIM = as.integer(unlist(strsplit(yaxisHIM, ",")))  
            plot_chim_segments(chimera_segment=seg, annotation, quality = "forestgreen", yaxis=yaxis,  yaxisHIM=yaxisHIM, HIM=TRUE)
        } else {
            plot_chim_segments(chimera_segment=seg, annotation, quality = "forestgreen", yaxis=yaxis,  yaxisHIM=yaxisHIM, HIM=FALSE)
        }
    }
}

dev.off()

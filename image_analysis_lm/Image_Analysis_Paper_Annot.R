#Crop HE images from spatial transcriptomics, to get spot-associated cell size
setwd("~/working_directory")
#load packages
library(dplyr)
packages <- c("imager", "magick", "ggplot2","data.table")
lapply(packages, library, character.only=TRUE)
#read data for baseline samples
samples <- c("42", "44", "46", "48", "49", "50", "51", "52", "54", "55")
#combine information for x and y coordinates on the images with cluster assignments
for (i in 1:length(samples)) {
coordinates <- read.table(paste("spatial_spot_coordinates/S",samples[i],"_tissue_positions_list.csv",sep=""), sep=",", col.names = c("barcode", "Selection", "X", "Y", "Xpixel", "Ypixel"))
subtype <- read.table(paste("seu_clusters-baseline_all-selected/seu_S",samples[i],"_clusters.csv",sep=""), sep=",", header= TRUE)
combined <- merge(coordinates, subtype, by="barcode")
write.table(combined, file = paste("Size_Analysis/Combined_Sample", samples[i],".csv",sep=""), row.names = F, col.names = T, sep=",")
}
#print coordinates for each subject, to see if orientation is the same as HE image
setwd("spatial_spot_coordinates/")
temp <- list.files(pattern = "*.csv")
read.all <- lapply(temp, read.csv, col.names = c("barcode", "Selection", "X", "Y", "Xpixel", "Ypixel"))
w <- 0
for (i in read.all) {
  coord <- ggplot(data = as.data.frame(i), aes(x=X, y=Y, color=Selection)) + geom_point()
  w <- w + 1
  ggsave(paste(temp[w],".jpeg",sep=""), device = "jpg")
  coord <- ggplot(data = as.data.frame(i), aes(x=as.integer(Y), y=as.integer(-X), color=Selection)) + geom_point()
  ggsave(paste("rotated",temp[w],".jpeg",sep=""), device = "jpg")
}
#read in combined files with coordinated and spot assignment
setwd("~/05_Data/Spatial_Transcriptomics/Size_Analysis/")
temp1 <- list.files(pattern = "*.csv")
read.all <- lapply(temp1, read.csv)
dir.create("Results") 
#set parameters
d <- 200 
#split into subgroups (LEP, ADIPOQ=PLIN, SAA), create directory, process images
dir.create("Results/Barcode")
dir.create("Results/Barcode/LEP")
dir.create("Results/Barcode/ADIPOQ")
dir.create("Results/Barcode/SAA")
for (i in 1:length(read.all)) {
  setwd("~/05_Data/Spatial_Transcriptomics/Size_Analysis/Results/Barcode/")
  HE_img <- image_read(paste("~/05_Data/Spatial_Transcriptomics/HE_staining/HE_images/S",samples[i],".jpg", sep=""))
  LEP <- read.all[[i]][which(read.all[[i]]$seu_clusters == 2),]
  #based on prior images, Y and X coordinates are actually flipped
  LEP$correctX <- (LEP$Ypixel)
  LEP$correctY <- (LEP$Xpixel)
  #magickR starts counting coordinates in the top-left corner, so to get an image with 200 px in each direction of spot center, subtract 200 px
  LEP$Xpxrange <- as.character(LEP$correctX - d)
  LEP$Ypxrange <- as.character(LEP$correctY - d)
  write.table(x = LEP, file = paste("./LEP/",samples[i],".csv",sep = ""), row.names = F, col.names = T, sep = ",")
  for (j in 1:nrow(LEP)) {
    x <- LEP$Xpxrange[j]
    y <- LEP$Ypxrange[j]
    temp_img <- image_crop(HE_img, paste("400x400+",x,"+",y,sep=""))
    barcode <- as.character(LEP$barcode[j])
    image_write(temp_img, path= paste("./LEP/LEP_S",samples[i],"_",barcode,".tiff",sep=""), format = "tiff")
  }
  ADIPOQ <- read.all[[i]][which(read.all[[i]]$seu_clusters == 3),]
  ADIPOQ$correctX <- (ADIPOQ$Ypixel)
  ADIPOQ$correctY <- (ADIPOQ$Xpixel)
  ADIPOQ$Xpxrange <- as.character(ADIPOQ$correctX - d)
  ADIPOQ$Ypxrange <- as.character(ADIPOQ$correctY - d)
  write.table(x = ADIPOQ, file = paste("./ADIPOQ/",samples[i],".csv",sep = ""), row.names = F, col.names = T, sep = ",")
  for (jj in 1:nrow(ADIPOQ)) {
    x <- ADIPOQ$Xpxrange[jj]
    y <- ADIPOQ$Ypxrange[jj]
    temp_img <- image_crop(HE_img, paste("400x400+",x,"+",y,sep=""))
    barcode <- as.character(ADIPOQ$barcode[jj])
    image_write(temp_img, path= paste("./ADIPOQ/ADIPOQ_S",samples[i],"_",barcode,".tiff",sep=""), format = "tiff")
  }
  SAA <- read.all[[i]][which(read.all[[i]]$seu_clusters == 5),]
  SAA$correctX <- (SAA$Ypixel)
  SAA$correctY <- (SAA$Xpixel)
  SAA$Xpxrange <- as.character(SAA$correctX - d)
  SAA$Ypxrange <- as.character(SAA$correctY - d)
  write.table(x = SAA, file = paste("./SAA/",samples[i],".csv",sep = ""), row.names = F, col.names = T, sep = ",")
  for (jjj in 1:nrow(SAA)) {
    x <- SAA$Xpxrange[jjj]
    y <- SAA$Ypxrange[jjj]
    temp_img <- image_crop(HE_img, paste("400x400+",x,"+",y,sep=""))
    barcode <- as.character(SAA$barcode[jjj])
    image_write(temp_img, path= paste("./SAA/SAA_S",samples[i],"_",barcode,".tiff",sep=""), format = "tiff")
  }
}
#determine adipocyte size using an imagej macro. Get area in pixel, x and y coordinates of center of mass
#
#copy to editor/notepad and save as .ijm
#run("Set Measurements...", "area center display redirect=None decimal=0");
#run("Duplicate...", "title=copy");
#run("8-bit");
#run("Subtract Background...", "rolling=50 light");
#run("Median...", "radius=1");
#setAutoThreshold("Huang dark");
#run("Analyze Particles...", "size=2000-Infinity circularity=0.20-1.00 exclude add");
#selectWindow("copy");
#close();
#roiManager("Measure");
#
#load tables per cluster
size.ADIPOQ <- read.table("./../Barcode/ADIPOQ/Results_ADIPOQ_Barcode.csv", sep=",", header = T)
size.LEP <- read.table("./../Barcode/LEP/Results_LEP_Barcode.csv", sep=",", header = T)
size.SAA <- read.table("./../Barcode/SAA/Results_SAA_Barcode.csv", sep=",", header = T)
#ADIPOQ/PLIN cluster
#calculate distance to center of image= center of spot (at x=200px, y=200px)
size.ADIPOQ$delta <- abs(200 - size.ADIPOQ$XM) + abs(200 - size.ADIPOQ$YM) 
#some formatting to get barcode and sample number
size.ADIPOQ$sample <- substr(size.ADIPOQ$Label, 8,10)
size.ADIPOQ$barcode <- gsub("ADIPOQ_S\\d+_","", size.ADIPOQ$Label)
size.ADIPOQ$barcode <- gsub("\\.tiff\\:\\d+-\\d+","", size.ADIPOQ$barcode)
size.ADIPOQ$sample <- gsub("ADIPOQ_","", size.ADIPOQ$Label)
size.ADIPOQ$sample <- gsub("_\\w+-1\\.tiff\\:\\d+-\\d+","", size.ADIPOQ$sample)
size.ADIPOQ$spot <- paste(size.ADIPOQ$sample,"_",size.ADIPOQ$barcode,sep="")
size.ADIPOQ <- data.table(size.ADIPOQ)
#only select the adipocyte with lowest delta per spot
min.size.ADIPOQ <- size.ADIPOQ[ , .SD[which.min(delta)], by = spot]
#sanity check and get some values to determine threshold
median(size.ADIPOQ$Area) #3954 px
avg.radius <- sqrt(3954/pi) #35 px; assuming circle
#filter for spots with max delta < 120 to reduce likelyhood of adipocytes not overlapping with the spot
min.size.ADIPOQ.120 <- min.size.ADIPOQ[min.size.ADIPOQ$delta < 120,]
#calculate for LEP and SAA with same parameters
size.LEP$delta <- abs(200 - size.LEP$XM) + abs(200 - size.LEP$YM)
size.LEP$sample <- substr(size.LEP$Label, 8,10)
size.LEP$barcode <- gsub("LEP_S\\d+_","", size.LEP$Label)
size.LEP$barcode <- gsub("\\.tiff\\:\\d+-\\d+","", size.LEP$barcode)
size.LEP$sample <- gsub("LEP_","", size.LEP$Label)
size.LEP$sample <- gsub("_\\w+-1\\.tiff\\:\\d+-\\d+","", size.LEP$sample)
size.LEP$spot <- paste(size.LEP$sample,"_",size.LEP$barcode,sep="")
size.LEP <- data.table(size.LEP)
min.size.LEP <- size.LEP[ , .SD[which.min(delta)], by = spot]
min.size.LEP.120 <- min.size.LEP[min.size.LEP$delta < 120,]
size.SAA$delta <- abs(200 - size.SAA$XM) + abs(200 - size.SAA$YM)
size.SAA$sample <- substr(size.SAA$Label, 8,10)
size.SAA$barcode <- gsub("SAA_S\\d+_","", size.SAA$Label)
size.SAA$barcode <- gsub("\\.tiff\\:\\d+-\\d+","", size.SAA$barcode)
size.SAA$sample <- gsub("SAA_","", size.SAA$Label)
size.SAA$sample <- gsub("_\\w+-1\\.tiff\\:\\d+-\\d+","", size.SAA$sample)
size.SAA$spot <- paste(size.SAA$sample,"_",size.SAA$barcode,sep="")
size.SAA <- data.table(size.SAA)
min.size.SAA <- size.SAA[ , .SD[which.min(delta)], by = spot]
min.size.SAA.120 <- min.size.SAA[min.size.SAA$delta < 120,]
#adding further filtering for max size; size over 25000 px unlikely to be adipocyte (might be due to broken cell walls, etc.)
min.size.ADIPOQ.120 <- min.size.ADIPOQ.120[min.size.ADIPOQ.120$Area < 25000,]
min.size.LEP.120 <- min.size.LEP.120[min.size.LEP.120$Area < 25000,]
min.size.SAA.120 <- min.size.SAA.120[min.size.SAA.120$Area < 25000,]
#assign cluster names and combine table
min.size.ADIPOQ.120$cluster <- "ADIPOQ"
min.size.LEP.120$cluster <- "LEP"
min.size.SAA.120$cluster <- "SAA"
size.filtered <- rbind(min.size.ADIPOQ.120, min.size.LEP.120, min.size.SAA.120)
#save results
write.table(size.filtered, file = "./../Barcode/Filtered_Size_WithBarcode.txt", row.names = F, col.names = T, sep="\t")
#calculate mean and median per individual
per.sample <- data.frame(matrix(vector(), 0, 4,
                                     dimnames=list(c(), c("Sample", "Cluster", "Median", "Mean"))),
                              stringsAsFactors=F)
samples2 <- unique(size.filtered$sample)
for (i in 1:10) {
  temp.table <- size.filtered[which(size.filtered$sample == samples2[i]),]
  temp.table.ADIPOQ <- temp.table[which(temp.table$cluster == "ADIPOQ"),]
  temp.table.LEP <- temp.table[which(temp.table$cluster == "LEP"),]
  temp.table.SAA <- temp.table[which(temp.table$cluster == "SAA"),]
  empty <- data.frame(matrix(vector(), 1, 4,
                             dimnames=list(c(), c("Sample", "Cluster", "Median", "Mean"))),
                      stringsAsFactors=F)
  empty$Sample = samples2[[i]]
  empty$Cluster = "ADIPOQ"
  empty$Median = median(temp.table.ADIPOQ$Area)
  empty$Mean = mean(temp.table.ADIPOQ$Area)
  per.sample <- rbind(per.sample, empty)
  empty <- data.frame(matrix(vector(), 1, 4,
                             dimnames=list(c(), c("Sample", "Cluster", "Median", "Mean"))),
                      stringsAsFactors=F)
  empty$Sample = samples2[[i]]
  empty$Cluster = "LEP"
  empty$Median = median(temp.table.LEP$Area)
  empty$Mean = mean(temp.table.LEP$Area)
  per.sample <- rbind(per.sample, empty)
  empty = data.frame(matrix(vector(), 1, 4,
                            dimnames=list(c(), c("Sample", "Cluster", "Median", "Mean"))),
                     stringsAsFactors=F)
  empty$Sample = samples2[[i]]
  empty$Cluster = "SAA"
  empty$Median = median(temp.table.SAA$Area)
  empty$Mean = mean(temp.table.SAA$Area)
  per.sample <- rbind(per.sample, empty)
}
write.table(per.sample, file = "~/05_Data/Spatial_Transcriptomics/Size_Analysis/Results/Size_per_Subject.txt", col.names = T, row.names = F, sep="\t")
#calculate cluster-independent mean and median per sample
median_per_ind <- aggregate(size.filtered[, 4], list(size.filtered$sample), median)
mean_per_ind <- aggregate(size.filtered[, 4], list(size.filtered$sample), mean)
#compare to cell volume data
volume <- read.table("./../Cell_volume.txt", header = T, sep="\t")
volume_area <- merge(median_per_ind, volume, by.x="Group.1", by.y = "Sample")
volume_area <- merge(volume_area, mean_per_ind, by="Group.1")
colnames(volume_area) <- c("Samples", "Median_Area", "Cell_Volume", "Mean_Area")
write.table(volume_area, file = "./../../../Area_Mean_Volume.txt", col.names = T, row.names = F)
ggplot(data = volume_area, aes(x=Area, y=cell_volume)) + geom_point() + geom_smooth(method='lm')
cor.test(volume_area$cell_volume, volume_area$Area, method= "pearson")

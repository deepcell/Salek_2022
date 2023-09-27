# The data used is available in GEO, with accession ID GSE241837
library(umap)
library(scales)

d <- read.csv("featureselection_normalized.csv", as.is=T)
names(d)[1] <- "cell_id"
epcam <- read.csv("epcam.csv", as.is = T, check.names = F)
names(epcam)[1] <- names(d)[1] <- "cell_id"
d <- merge(d, epcam)

um1 <- umap(d[,2:ncol(d)])
um1_layout <- data.frame(um1$layout)
um1_layout$cell_id <- d$cell_id
with(um1_layout, plot(X1, X2, pch = 19, col = alpha("darkgoldenrod",0.5) , cex = 0.1, main = "n_neighbors = 15, min_dist = 0.1"))
write.csv(um1_layout, "umap_nn_15_mindist_p1_epcam.csv")

d2 <- d[!(d$cell_id %in% um1_layout$cell_id[um1_layout$X2 < -10]),]
um2 <- umap(d2[,2:ncol(d2)])
um2_layout <- data.frame(um2$layout)
um2_layout$cell_id <- d2$cell_id
with(um2_layout, plot(X1, X2, pch = 19, col = alpha("darkgoldenrod",0.5) , cex = 0.1, main = "n_neighbors = 15, min_dist = 0.1"))
write.csv(um1_layout, "umap2_nn_15_mindist_p1_epcam.csv")

custom.config <- umap.defaults
custom.config$min_dist <- 0.3
um2 <- umap(d2[,2:ncol(d2)], custom.config)
um2_mindist_p3_layout <- data.frame(um2$layout)
um2_mindist_p3_layout$cell_id <- d2$cell_id
with(um2_mindist_p3_layout, plot(X1, X2, pch = 19, col = alpha("darkgoldenrod",0.5) , cex = 0.1, main = "n_neighbors = 15, min_dist = 0.3"))
write.csv(um2_mindist_p3_layout, "um2_mindist_p3_layout_epcam.csv")

custom.config <- umap.defaults
custom.config$min_dist <- 0.3
um_nn15_mindist_p3 <- umap(d[,2:ncol(d)], custom.config)
um_nn15_mindist_p3_layout <- data.frame(um_nn15_mindist_p3$layout)
um_nn15_mindist_p3_layout$cell_id <- d$cell_id
with(um_nn15_mindist_p3_layout, plot(X1, X2, pch = 19, col = alpha("darkgoldenrod",0.5) , cex = 0.1, main = "n_neighbors = 15, min_dist = 0.3"))
write.csv(um_nn15_mindist_p3_layout, "um_nn15_mindist_p3_epcam.csv")

coul <- colorRampPalette(c("lightgrey","red"))(100)
coul2 <- colorRampPalette(c("lightgrey","blue"))(100)
ul <- um2_mindist_p3_layout
ul <- merge(ul, epcam)
ul <- merge(ul, d2[,c("cell_id","PTPRC")])
ul$Col <- coul[as.numeric(cut(log10(ifelse(ul$EPCAM == 0,1, ul$EPCAM)),breaks = 100))]
ul$Col2 <- coul2[as.numeric(cut(log10(ifelse(ul$PTPRC == 0, 1, ul$PTPRC)),breaks = 100))]
ul1 <- ul[grep("DTC1", ul$cell_id),]
ul2 <- ul[grep("DTC2", ul$cell_id),]
ul3 <- ul[grep("DCT3", ul$cell_id),]
png("DTC_low_epcam.png")
with(ul1, plot(X1, X2, pch = 19, col = "lightgrey", cex = 0.1, main = "", xlab = "umap_1", ylab = "umap_2", asp = 1))
with(ul1[ul1$PTPRC > 0,], points(X1, X2, pch = 19, col = alpha(Col2,0.5), cex = 0.3))
with(ul1[ul1$EPCAM > 0,], points(X1, X2, pch = 19, col = alpha(Col,0.5), cex = 0.3))
dev.off()
png("DTC_intermediate_epcam.png")
with(ul3, plot(X1, X2, pch = 19, col = "lightgrey", cex = 0.1, main = "", xlab = "umap_1", ylab = "umap_2"))
with(ul3[ul3$PTPRC > 0,], points(X1, X2, pch = 19, col = Col2, cex = 0.3))
with(ul3[ul3$EPCAM > 0,], points(X1, X2, pch = 19, col = Col, cex = 0.3))
dev.off()
png("DTC_high_epcam.png")
with(ul2, plot(X1, X2, pch = 19, col = "lightgrey", cex = 0.1, main = "", xlab = "umap_1", ylab = "umap_2"))
with(ul2[ul2$PTPRC > 0,], points(X1, X2, pch = 19, col = Col2, cex = 0.3))
with(ul2[ul2$EPCAM > 0,], points(X1, X2, pch = 19, col = Col, cex = 0.3))
dev.off()

svg("color_bar.svg")
color.bar <- function(lut, min, max=-min, nticks=(max-min)*10+1, ticks=seq(min, max, len=nticks), title='') {
  scale = (length(lut)-1)/(max-min)
  plot(c(0,1), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title, asp = 20)
  axis(2,  las=1)
  for (i in 1:(length(lut)-1)) {
    y = (i-1)/scale + min
    rect(0,y,10,y+1/scale, col=lut[i], border=NA)
  }
}

color.bar(coul, 0, round(max(log10(ul$EPCAM)),1))
color.bar(coul2, 0, round(max(log10(ul$PTPRC)),1))
dev.off()

custom.config <- umap.defaults
custom.config$n_neighbors <- 30
um_nn30_mindist_p1 <- umap(d[,2:ncol(d)], custom.config)
um_nn30_mindist_p1_layout <- data.frame(um_nn30_mindist_p1$layout)
um_nn30_mindist_p1_layout$Cell_ID <- d$cell_id
with(um_nn30_mindist_p1_layout, plot(X1, X2, pch = 19, col = alpha("darkgoldenrod",0.5) , cex = 0.3, main = "n_neighbors = 30, min_dist = 0.1"))
write.csv(um_nn30_mindist_p1_layout, "umap_nn_30_mindist_p1.csv")

custom.config <- umap.defaults
custom.config$n_neighbors <- 60
um_nn60_mindist_p1 <- umap(d[,2:ncol(d)], custom.config)
um_nn60_mindist_p1_layout <- data.frame(um_nn60_mindist_p1$layout)
um_nn60_mindist_p1_layout$Cell_ID <- d$cell_id
with(um_nn60_mindist_p1_layout, plot(X1, X2, pch = 19, col = alpha("darkgoldenrod",0.5) , cex = 0.3, main = "n_neighbors = 60, min_dist = 0.1"))
write.csv(um_nn60_mindist_p1_layout, "um_nn_60_mindist_p1.csv")

custom.config <- umap.defaults
custom.config$min_dist <- 0.3
um_nn15_mindist_p3 <- umap(d[,2:ncol(d)], custom.config)
um_nn15_mindist_p3_layout <- data.frame(um_nn15_mindist_p3$layout)
um_nn15_mindist_p3_layout$Cell_ID <- d$cell_id
with(um_nn15_mindist_p3_layout, plot(X1, X2, pch = 19, col = alpha("darkgoldenrod",0.5) , cex = 0.1, main = "n_neighbors = 15, min_dist = 0.3"))
write.csv(um_nn15_mindist_p3_layout, "um_nn15_mindist_p3.csv")

custom.config <- umap.defaults
custom.config$min_dist <- 0.9
um_nn15_mindist_p9 <- umap(d[,2:ncol(d)], custom.config)
um_nn15_mindist_p9_layout <- data.frame(um_nn15_mindist_p9$layout)
um_nn15_mindist_p9_layout$Cell_ID <- d$cell_id
with(um_nn15_mindist_p9_layout, plot(X1, X2, pch = 19, col = alpha("darkgoldenrod",0.5) , cex = 0.3, main = "n_neighbors = 15, min_dist = 0.9"))
write.csv(um_nn15_mindist_p9_layout, "um_nn15_mindist_p9.csv")

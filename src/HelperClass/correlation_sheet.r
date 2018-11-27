library(cluster)
require(lattice)
library(corrplot)
args <- commandArgs(TRUE)
dat_path <- args[1]
colname_path <- args[2]
pdf_path <- args[3]

x <- scan(colname_path, character(), quote = " ")
print(c(x))

corr <- read.table(dat_path, row.names=1)
if(nrow(corr)==16)
{
	colnames(corr) <- c("kmeans", "kmeans_", "kmedoids", "kmedoids_", "AHC-single", "AHC-single_", "AHC-average", "AHC-average_", "BIRCH", "BIRCH_", 
		"DBSCAN", "DBSCAN_", "OPTICS", "OPTICS_", "SC-kmeans", "SC-kmeans_", "SC-eigen", "SC-eigen_", "AP", "AP_")
}else{
	#colnames(corr) <- c("kmeans", "kmedoids", "AHC-single", "AHC-average", "BIRCH", "DBSCAN", "OPTICS", "SC-kmeans", "SC-eigen", "AP")
	colnames(corr) <- c(x)
}
# rownames(corr) <- 0:(nrow(corr)-1)
# colnames(corr) <- 0:(ncol(corr)-1)
bwr <- c("white", "white", "green")
#bwr <- c("white", "#4080FF", "red")
col3 <- colorRampPalette(bwr, space = "rgb")

pdf(pdf_path, height=25, width=25)
lattice.options(axis.padding=list(factor=0.5))

x.scale <- list(cex=2, alternating=2, col='black',rot=90)
y.scale <- list(cex=2, alternating=1, col='black')

q  = t(as.matrix(corr))
q2 = q[,] 
#corrplot(q2,col=col3(100), tl.col="black", tl.srt=45)
#corrplot(q2,col=col3(100), tl.col="black", tl.srt=45, order = "FPC")
#corrplot(q2,col=col3(100), tl.col="black", tl.srt=45, order = "AOE")
#corrplot(q2,col=col3(100), tl.col="black", tl.srt=45, order = "hclust")

#corrplot(q2,col=col3(100), tl.col="black", tl.srt=45, type = "upper", tl.cex=50/nrow(corr), cl.pos = "n")
corrplot(q2,col=col3(100), tl.col="black", tl.srt=45, tl.cex=50/nrow(corr), cl.pos = "n")

#image(q2, col=col3(100))
dev.off()



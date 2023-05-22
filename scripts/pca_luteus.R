####READ ME####
#this is the pca set up for aciurina mfc data set
####load packages ####
library("adegenet")
library("ade4")
library("vcfR")
library("parallel")
library("viridis")
library("paletteer")
####adegenet’s glPCA####
#This is run in Rstudio (it being an IDE is especially nice for visualizing)
#Read in VCF, one at a time!
test <- read.vcfR("vcf/populations_lc.snps.vcf")

#Convert to genlight format
focal_gl <- vcfR2genlight(test)
#populations map file
popmapbad<-data.frame(id=colnames(test@gt)[2:length(colnames(test@gt))],pop=substr(colnames(test@gt)[2:length(colnames(test@gt))], 7, 8))
#or read it in, one at a time!
popmap<-read.csv("popmaps/popmap_lowcov.csv")

#Add population info, setting up the population array
#popmap5
pop(focal_gl) <- c("NE", "NE", "NE", "NE", "NE", "NE", "NE", "NE", "NW", "NW", "NW", "NW", "NW", "NW", "NW", "NW", "rio", "rio", "rio", "rio", "rio", "rio", "rio", "rio", "sacramento", "sacramento", "sacramento", "sacramento", "wht", "wht", "wht", "wht", "wht", "wht", "wht", "wht", "wht", "rio")
#popmap8
pop(focal_gl) <- c("jhm", "jhm", "jhm", "jhm", "sdc", "sdc", "sdc", "sdc", "sjn", "sjn", "sjn", "sjn", "jmz", "jmz", "jmz", "jmz", "isa", "isa", "isa", "isa", "bda", "bda", "bda", "bda", "sac", "sac", "sac", "sac", "wht", "wht", "wht", "wht", "sjn", "sjn", "wht", "wht", "wht", "wht", "wht", "isa")
#Run PCA
focal_pca <- glPca(focal_gl, n.cores=4, nf=4)
#Not really great, but can be helpful for visualizing potentially problematic samples
scatter(focal_pca, cex=.25)

#My preferred option, I generally make plots here and edit in illustrator
s.class(focal_pca$scores[,c(1,2)], pop(focal_gl), 
        col=magma(10, begin=.8, end=0), clab=, cell=2.5)
#col is the color pallet used, magma is a nice default but try viridis
#cell is effectively a confidence interval (2.5 ~= 95%)

#barplot plots the proportion of the total variance explained by given PCs
#Generally, the more a PC sticks out the more it matters
barplot(focal_pca$eig/sum(focal_pca$eig), main="eigenvalues", 
        col=heat.colors(length(focal_pca$eig)))

## Heat map of genotype
#glPlot(focal_gl)

pca <- glPca(focal_gl, nf = 30)
# Quick plot
scatter(pca, posi = "none")
# Create named color vector for the legend
#legend_colors <- rainbow(10)

legend_colors <- c("turquoise","seagreen","gold", "coral", "hotpink")
legend_colors <- c("turquoise","seagreen","gold", "coral", "hotpink","lightgreen", "navy",  "darkmagenta")

names(legend_colors) <- c("NE", "NW", "rio", "sacramento", "wht")    
names(legend_colors) <- c("bda", "isa", "jhm", "jmz", "sac", "sdc", "sjn", "wht")   
#levels5: NE NW rio sacramento wht -> "NE", "NW", "rio", "sacramento", "wht"
#levels8: bda isa jhm jmz sac sdc sjn wht -> "bda", "isa", "jhm", "jmz", "sac", "sdc", "sjn", "wht"

popmap$pop=as.factor(popmap$pop)

# Plot the first two principal components
plot(x = pca$scores[, 1], y = pca$scores[, 2], col = legend_colors [as.numeric(popmap$pop)], cex = 1, pch = 19)

# Add a legend
legend("bottomleft", legend = names(legend_colors), fill = legend_colors,
       title = "Population", bty = "n", ncol = 1, box.lwd = 0, box.col = "white", cex = 0.8)

# Plot the second and third principal components
plot(x = pca$scores[, 2], y = pca$scores[, 3], col = legend_colors [as.numeric(popmap$pop)], cex = 1, pch = 19)

# Add a legend
legend("bottomleft", legend = names(legend_colors), fill = legend_colors,
       title = "Population", bty = "n", ncol = 1, box.lwd = 0, box.col = "white", cex = 0.8)

# Plot the first and third principal components
plot(x = pca$scores[, 1], y = pca$scores[, 3], col = legend_colors [as.numeric(popmap$pop)], cex = 1, pch = 19)

# Add a legend
legend("topleft", legend = names(legend_colors), fill = legend_colors,
       title = "Population", bty = "n", ncol = 1, box.lwd = 0, box.col = "white", cex = 0.8)

pca_scores <- pca$scores
write.csv(pca_scores, "pca_scores8.csv")

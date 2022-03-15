library(affy)
library(limma)
library(GEOquery)
library(data.table)
library(hgu133plus2.db)
library(RColorBrewer)



#### Reproducing paper having DOI: 10.1155/2019/2478453 

##### Downloading CEL files, Untar and listing for the 1st time

my.gse <- "GSE22255"

my.geo.gse <- getGEO(GEO=my.gse, filename=NULL, destdir="/scratch/bcbf/gourab/projects/stroke/datasets", GSElimits=NULL, GSEMatrix=TRUE, AnnotGPL=TRUE, getGPL=TRUE)

#getGEOSuppFiles(my.gse, makeDirectory=F, baseDir="/scratch/bcbf/gourab/projects/stroke/datasets")   ###### Time consuming in R. Local download and transfer is preferable
untar(paste0("/scratch/bcbf/gourab/projects/stroke/datasets/",my.gse,"_RAW.tar"), exdir=paste0("/scratch/bcbf/gourab/projects/stroke/datasets/",my.gse,"/CEL"))
#my.cels <- list.files(paste0("/scratch/bcbf/gourab/projects/stroke/datasets/",my.gse,"/CEL"), pattern=".CEL")
my.cel <- as.data.frame(list.files(paste0("/scratch/bcbf/gourab/projects/stroke/datasets/",my.gse,"/CEL"), pattern=".CEL"))
colnames(my.cel) <- "CEL_File"
my.cels <- my.cel[grep(".CEL", my.cel$CEL_File), ]

##### Exploring the Data

class(my.geo.gse)
length(my.geo.gse)
names(my.geo.gse)
my.geo.gse <- my.geo.gse[[1]]
class(my.geo.gse)
colnames(pData(my.geo.gse))
pData(my.geo.gse)$data_processing[1]


####### Describing the data

head(exprs(my.geo.gse))
summary(exprs(my.geo.gse))


###### Analyzing the Data

##### Creating pheotype data

my.pdata <- as.data.frame(pData(my.geo.gse), stringsAsFactors=F)
#head(my.pdata)
#dim(my.pdata)

###### Selecting only patient data

my.pdata <- my.pdata[, c("title", "geo_accession", "description", "gender:ch1", "affected status (disease state):ch1")]
my.pdata <- my.pdata[order(rownames(my.pdata)), ]
rownames(my.pdata) <- my.cels
table(rownames(my.pdata) == my.cels)
write.table(my.pdata, file=paste0("/scratch/bcbf/gourab/projects/stroke/datasets/",my.gse,"/CEL/",my.gse,"_SelectPhenoData.txt"), sep="\t", quote=F)

######## DEG analysis

cel.path <- paste0("/scratch/bcbf/gourab/projects/stroke/datasets/",my.gse,"/CEL")

#### To reproduce Zhu et al. result, Move manually the control CEL files

my.affy <- ReadAffy(celfile.path=cel.path, phenoData=paste(cel.path, paste0(my.gse,"_SelectPhenoData.txt"), sep="/"))
pData(my.affy)$sample.levels <- c(rep("Control", 20), rep("Stroke", 20))
pData(my.affy)$sample.labels <- c(paste("Control", 1:20, sep="."), paste("Stroke", 21:40, sep="."))



### Shortening the annotation


#dim(exprs(my.affy))
#pData(my.affy)




##Calculate gene level expression measures
my.rma <- rma(my.affy, normalize=T, background=T)

####### Plotting Density plot of rma values show need for normalization

##make sample.levels a factor

pData(my.rma)$sample.levels <- as.factor(pData(my.rma)$sample.levels)

level.pal <- brewer.pal(3, "Dark2")
level.cols <- level.pal[unname(pData(my.rma)$sample.levels)]


#plotDensities(exprs(my.rma), legend=F, col=level.cols, main="Arrays Not Normalized")
#legend("topright", legend=levels(pData(my.rma)$sample.levels), fill=level.pal)

##Boxplot is another way of showing the same thing.

boxplot(exprs(my.rma), las=2, names=pData(my.rma)$sample.labels, outline=F, col=level.cols, main="Arrays Normalized")


###### Design specification


my.design <- model.matrix(~0 + sample.levels, pData(my.rma))
rownames(my.design) <- pData(my.rma)$sample.labels
colnames(my.design) <- levels(pData(my.rma)$sample.levels)

##determine the average effect (coefficient) for each treatment

my.fit <- lmFit(my.rma, my.design)

##### Making Comparisons (Contrasts) Between Samples

contrast.matrix = makeContrasts(contrasts=c("Control-Stroke"), levels=my.design)
contrast.fits <- contrasts.fit(my.fit, contrast.matrix)

contrast.ebs <- eBayes(contrast.fits, proportion=0.1, trend=FALSE, robust=FALSE)
contrast.tts <- topTable(contrast.ebs, adjust="BH", number=length(contrast.ebs$coefficients), sort.by="none")
contrast.tests <- as.data.frame(decideTests(contrast.ebs, method="separate", adjust.method="BH", p.value=0.05, lfc=0))

########### Plotting DEGs

my.contrast <- "Control-Stroke"
degs <- contrast.tts[contrast.tts$logFC >= 1.58 & contrast.tts$adj.P.Val < 0.05, ]
#updegs <- rownames(contrast.tests[[my.contrast]])[contrast.tests[[my.contrast]][, 1] == 1]
#downdegs <- rownames(contrast.tests[[my.contrast]])[contrast.tests[[my.contrast]][, 1] == -1]

updegs <- rownames(contrast.tests)[contrast.tests$`ISF-ISM` == 1]
downdegs <- rownames(contrast.tests)[contrast.tests$`ISF-ISM` == -1]



######## Annotating with gene ids, symbols etc.

gene.data <- select(hgu133plus2.db, keys=rownames(contrast.tts),
                    keytype = "PROBEID", columns = c("ENTREZID", "GENENAME", "SYMBOL"))

contrast.tts$PROBEID <- rownames(contrast.tts)
rownames(contrast.tts) <- 1:nrow(contrast.tts)


contrast.tests$PROBEID <- rownames(contrast.tests)
rownames(contrast.tests) <- 1:nrow(contrast.tests)
df_tmp <- merge(contrast.tts, contrast.tests)

############ Final Output

df <- merge(gene.data, df_tmp)
write.table(df, file=paste0("/scratch/bcbf/gourab/projects/stroke/datasets/",my.gse,"_DEG_res.txt"), sep="\t", quote=F, row.names = F)

### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###

### DGE analysis for Arabidopsis RNA-Seq samples in Applied Genome Research course ###

library("DESeq2")

# --- loading sampleTable --- #
csvfile <- "/prj/gf-arabseq/project_Ath-Nd1/members/bpucker/GRC/tests/deseq2/clean_sample_table.txt"
sampleTable <- read.csv(csvfile, row.names=1, sep="\t")
sampleTable$genotype <- as.factor( sampleTable$genotype )
summary(sampleTable)

# --- loading the data matrix --- #
count_data_file <- "/prj/gf-arabseq/project_Ath-Nd1/members/bpucker/GRC/tests/deseq2/clean_data_matrix.txt"
countdata <- read.csv(count_data_file,row.names=1, header=T, sep="\t")
summary(countdata)

# --- construction of DESeqDataSet --- #
ddsMat <- DESeqDataSetFromMatrix( countData=countdata, colData=sampleTable, design= ~ genotype )
nrow(ddsMat)

# -- removal of not or low expressed genes --- #
dds <- ddsMat[ rowSums(counts(ddsMat)) > 100, ]
nrow(dds)

# --- plot PCA in R studio --- #
rld <- rlog(dds)
ramp <- 1:2/2
cols <- c( rgb(ramp, 0, 0), rgb(0, ramp, 0), rgb(ramp, 0, ramp), rgb(ramp, 0, ramp) )
print ( plotPCA( rld, intgroup=c( "genotype")  ) )

# --- differential expression analysis --- #
dds <- DESeq(dds)
res <- results(dds)
summary(res)

# --- investigate differentially expressed genes --- #
res.05 <- results( dds, alpha=.05  )
table(res.05$padj < .05  )

small.pvalue.index <- head( order( res$padj ), 20  )
names <- row.names(res)
( sig.gene.names <- names[ small.pvalue.index  ] )

outputfile <- "/prj/gf-arabseq/project_Ath-Nd1/members/bpucker/GRC/tests/deseq2/differentially_expressed_genes.txt"
write( sig.gene.names, outputfile, ncolumns=length( sig.gene.names ), sep="\n"  )

# --- print session information --- #
sessionInfo()

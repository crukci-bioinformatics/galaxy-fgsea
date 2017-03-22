options( show.error.messages=F, error = function () { cat( geterrmessage(), file=stderr() ); q( "no", 1, F ) } )

# we need that to not crash galaxy with an UTF8 error on German LC settings.
loc <- Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")

suppressPackageStartupMessages({
  library("fgsea")
  library("optparse")
  library("ggplot2")
})

option_list <- list(
  make_option(c("-rnk_file", "--rnk_file"), type="character", help="Path to file with differential gene expression result"),
  make_option(c("-out_tab","--out_tab"), type="character", help="Path to output file."),
  make_option(c("-gmt_file", "--gmt_file"), default="h.all.v5.2.symbols.gmt", type="character", help = "Path to Broad gmt file"),
  make_option(c("-min_size", "--min_size"), default=1, help="Minimal size of a gene set to test. All pathways below the threshold are excluded."),
  make_option(c("-max_size", "--max_size"), default=500, help="Maximal size of a gene set to test. All pathways above the threshold are excluded."),
  make_option(c("-n_perm", "--n_perm"),default=1000, help="Number of permutations to do. Minimial possible nominal p-value is about 1/nperm"),
  make_option(c("-top_n","--top_n"), default=10, help="The number of gene sets to produce diagnostics plots for. The top N up-regulated and top N down-regulated sets will be shown"),
  make_option(c("-summary_plot","--summary_plot"), type="character", help="Path to summary plot file."),
  make_option(c("-individual_plot","--individual_plot"), type="character", help="Path to individual plots file."),
  make_option(c("-file_has_header","--file_has_header"),type="logical",default=TRUE,help="If this option is set to TRUE, the tool will assume that the ranked gene-list has a column heading and the gene names commence on the second line")
)

parser <- OptionParser(usage = "%prog [options] file", option_list=option_list)
args = parse_args(parser)

# Vars:
rnk_file = args$rnk_file
#category_file = args$category_file
#length_file = args$length_file
gmt_file = args$gmt_file
out_tab = args$out_tab
min_size = args$min_size
max_size = args$max_size
n_perm = args$n_perm
top_n = args$top_n
summary_plot = args$summary_plot
individual_plot = args$individual_plot
file_has_header = args$file_has_header
### Change to whatever path the gmt files are stored in
path_to_gmt <- "/data/galaxy/git-galaxy-fgsea/"

gmt_file <-paste(path_to_gmt, gmt_file, sep="/")


## If testing locally, change to TRUE and arguments will be set below
run_local <- FALSE


if (run_local) {

  rnk_file <- "testdata/t47d_Treatment_DEA_Prog-vs-Control_all_for_GSEA.rnk"
  gmt_file <- "h.all.v5.2.symbols.gmt"
  path_to_gmt <- "."
  min_size = 5
  max_size = 500
  n_perm=1000
  top_n = 10
  out_tab = "temp.csv"
  summary_plot = "summary.pdf"
  individual_plot = "gene-sets.pdf"
  file_has_header = TRUE
}


## Basically using the steps from the fgsea vignette from now on

rankTab <- read.table(rnk_file,
                      header=file_has_header, colClasses = c("character", "numeric"))
ranks <-rankTab[,2]
names(ranks) <- rankTab[,1]
head(ranks)

## Report an error if gmt_file not found
if(file.exists(gmt_file)) {
  pathways <- gmtPathways(gmt_file)
} else {
  cat(paste("Could not find file ", gmt_file, "in folder", path_to_gmt))
  cat("Printing contents of directory:")
  list.files()
  getwd()
}

fgseaRes <- fgsea(pathways, ranks, minSize=min_size, maxSize=max_size, nperm=n_perm)
dim(fgseaRes)
head(fgseaRes)
cat(paste("Attempting to write to file", out_tab))

## 8th column of the output is a list, so cannot be put easily into a data frame
## http://stackoverflow.com/questions/17291283/outputting-a-dataframe-in-r-to-a-csv

my.df <- data.frame(lapply(fgseaRes, as.character), stringsAsFactors=FALSE)
my.df[,8] <- gsub("c(","",my.df[,8],fixed=TRUE)
my.df[,8] <- gsub(")","",my.df[,8],fixed=TRUE)
my.df[,8] <- gsub("\n", "",my.df[,8],fixed=TRUE)
my.df[,8] <- gsub("\"", "", my.df[,8],fixed=TRUE)
write.table(my.df, out_tab, row.names=FALSE,sep="\t",quote=FALSE)

## Now make a summary plot, and some plots of particular gene sets

topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=top_n), pathway]
topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=top_n), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))


pdf(summary_plot)
plotGseaTable(pathways[topPathways], ranks, fgseaRes, 
              gseaParam = 0.5)
dev.off()

pdf(individual_plot)
for(i in 1:top_n){
  
  p <- plotEnrichment(pathways[[topPathways[i]]],
                 ranks) + ggtitle(topPathways[i])
  print(p)
}
dev.off()

sessionInfo()

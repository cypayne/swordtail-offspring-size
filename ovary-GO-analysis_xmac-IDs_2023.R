## Project: Offspring size evolution
## Payne et al
##
## Gene Ontology enrichment analysis of X. malinche and X. birchmanni ovaries
## Takes in DESeq2 DGE output for all genes to build universe and for significantly
## differentially expressed genes between species
##
## cyp XII-2023

## load libraries
library("GOstats")
library("GSEABase")
library("biomaRt")

## set output directory
go_out_dir <- "Data/embryo_ovary_dge_combined_2023/ovary_xmac-gtf_dge/GO_analysis"

## load ovary DGE
ovary_dge <- read.csv("Data/embryo_ovary_dge_combined_2023/ovary_xmac-gtf_dge/ovary-xmacID-combined2023_dge_lfc-shr_all.csv_with-annots.csv",header=T)

## create Ensembl ID to GO term map
mart <- useEnsembl(biomart = "ensembl",
                   dataset = "xmaculatus_gene_ensembl",
                   mirror = "useast")

# match go ids to ensembl gene ids
results <- getBM(attributes = c("go_id","external_gene_name","ensembl_gene_id"), filters=c("ensembl_gene_id"),values=ovary_dge$Gene, mart = mart)

## create gene universe
# subset gene universe to only include genes with valid go ids and external gene names
gene_universe<-subset(results,nchar(results$go_id)>0 & nchar(results$external_gene_name) >0)
gene_universe$ensembl_id<-gene_universe[,3]
gene_universe[,3]<-as.numeric(as.factor(gene_universe[,3]))
gene_universe$Evidence<-rep("ISA",length(gene_universe[,3]))
colnames(gene_universe)<-c("frame.go_id","frame.gene_name","frame.gene_id","frame.ensembl","frame.Evidence")

goframeData <- data.frame(gene_universe$frame.go_id,gene_universe$frame.Evidence,gene_universe$frame.gene_id)

goFrame <- GOFrame(goframeData,organism="Xiphophorus")
goAllFrame <- GOAllFrame(goFrame)
gsc <- GeneSetCollection(goAllFrame, setType = GOCollection())

universe <- goframeData$gene_universe.frame.gene_id

## enrich for GO categories
# set p-value threshold
padj_cutoff = 0.05
# choose one of the following subsets to test for GO enrichment
# ex. by choosing padj0.05_malvbir_late, you will look for GO enrichment in the
#     set of genes that were differentially expressed between late-pregnancy
#     X. malinche and late-pregnancy X. birchmanni ovaries
padj0.05_malvbir_late <- subset(ovary_dge,padj_res.malvbirch_late<padj_cutoff)
padj0.05_malvbir_early <- subset(ovary_dge,padj_res.malvbirch_early<padj_cutoff)
padj0.05_malvbir_late_not_early <- padj0.05_malvbir_late[!(row.names(padj0.05_malvbir_late) %in% row.names(padj0.05_malvbir_early)),]
padj0.05_malvbir_late_and_early <- padj0.05_malvbir_late[(row.names(padj0.05_malvbir_late) %in% row.names(padj0.05_malvbir_early)),]

# XXX set this value to the subset you want to test for GO enrichment
ovary_dgesig <- padj0.05_malvbir_late_and_early
genes_sig <- ovary_dgesig$Gene
genes_match<-gene_universe[gene_universe$frame.ensembl %in% genes_sig,]
genes_match_sig <- genes_match$frame.gene_id

# optional: write the gene subset's GO information to an outfile
#write.csv(genes_match,"Data/embryo_ovary_dge_combined_2023/ovary_xmac-gtf_dge/ovary_dge_combined-FebAug2023_gene2GO-matchup.csv")

# set outfile header
of_header <- "GO_ovary-xmacID-combined2023_dge_padj0.05malvbir-late-and-early"

## loop through the three GO categories to output GO enrichment results for each
for(ontology_category in c("BP","MF","CC")) {
  # set params for hyperGTest
  # testDirection = "over" for overrepresented genes, = "under" for underrepresented
  # interested in overrep genes
  # three categories:
  #   cellular component (CC; where gene products are active)
  #   molecular function (MF; the biological function of gene or gene product)
  #   biological process (BP; pathways or larger processes that multiple gene products involved in).
  # Output:
  # ExpCount is the expected count and the Count is how many instances of that term were actually oberved
  # in your gene list while the Size is the number that could have been found in your gene list if every
  # instance had turned up. Values like the ExpCount and the Size are going to be affected by what is included
  # in the gene universe as well as by whether or not it was a conditional test.
  params <- GSEAGOHyperGParams(name="Xiphophorus maculatus genes",
                               geneSetCollection=gsc,
                               geneIds = genes_match_sig,
                               universeGeneIds = universe,
                               ontology = ontology_category,
                               pvalueCutoff = 0.05,
                               conditional = FALSE,
                               testDirection = "over")

  overrep <- hyperGTest(params)
  results_overrep <- summary(overrep)

  ## output an additional column with the genes falling under each GO category
  if(ontology_category=="BP") {
    res_GOID <- results_overrep$GOBPID
  } else if(ontology_category=="MF") {
    res_GOID <- results_overrep$GOMFID
  } else if(ontology_category=="CC") {
    res_GOID <- results_overrep$GOCCID
  }

  cats <- geneIdsByCategory(overrep)
  # get list of lists: all geneids and names per GOBPID
  extract_geneids <- sapply(res_GOID, function(go) cats[[go]])
  list_of_GOID_geneID_lists <- list()
  list_of_GOID_genename_lists <- list()
  for(gene_list in extract_geneids) {
    geneID_list <- unique(subset(genes_match, frame.gene_id %in% gene_list)$frame.ensembl)
    genename_list <- unique(subset(genes_match, frame.gene_id %in% gene_list)$frame.gene_name)
    geneID_list_str <- paste(unlist(geneID_list), sep=",", collapse=",")
    genename_list_str <- paste(unlist(genename_list), sep=",", collapse=",")
    list_of_GOID_geneID_lists <- append(list_of_GOID_geneID_lists,geneID_list_str)
    list_of_GOID_genename_lists <- append(list_of_GOID_genename_lists,genename_list_str)
  }
  # this gives dataframe of three columns: GOID, the list of gene ids, and the list of gene names falling under it
  gobpid2geneid <- tryCatch(
    {
      data.frame(GOID=names(extract_geneids), GENEID_list=unlist(list_of_GOID_geneID_lists), GENENAME_list=unlist(list_of_GOID_genename_lists))
    },
    error=function(e) {
      print(paste("no output for ", f, " category ",ontology_category))
      next
    } )
  # merge GO_GENEID_list with the rest of GO results
  if(ontology_category=="BP") {
    go_output <- merge(results_overrep,gobpid2geneid, by.x="GOBPID", by.y="GOID", all=TRUE)
  } else if(ontology_category=="MF") {
    go_output <- merge(results_overrep,gobpid2geneid, by.x="GOMFID", by.y="GOID", all=TRUE)
  } else if(ontology_category=="CC") {
    go_output <- merge(results_overrep,gobpid2geneid, by.x="GOCCID", by.y="GOID", all=TRUE)
  }
  # write output
  write.table(go_output, file.path(go_out_dir,paste0(of_header,"_pval0.05_",ontology_category,".tsv")), row.names = FALSE, quote=FALSE, sep="\t")
}

## ----setup, echo = FALSE--------------------------------------------------------------------------
knitr::opts_chunk$set(error = TRUE, cache = FALSE, eval = TRUE)

## ----annotate,echo=FALSE--------------------------------------------------------------------------
options(width=100)

## ----biomaRt--------------------------------------------------------------------------------------
library("biomaRt")
listMarts()

## ----ensembl1-------------------------------------------------------------------------------------
ensembl <- useMart("ensembl")

## ----listDatasets---------------------------------------------------------------------------------
datasets <- listDatasets(ensembl)
head(datasets)

## ----ensembl2, eval=TRUE--------------------------------------------------------------------------
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)

## ----ensembl3, eval = FALSE-----------------------------------------------------------------------
#  ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")

## ----filters--------------------------------------------------------------------------------------
filters = listFilters(ensembl)
filters[1:5,]

## ----attributes-----------------------------------------------------------------------------------
attributes = listAttributes(ensembl)
attributes[1:5,]

## ----getBM1, echo=TRUE, eval=TRUE-----------------------------------------------------------------
affyids=c("202763_at","209310_s_at","207500_at")
getBM(attributes=c('affy_hg_u133_plus_2', 'entrezgene_id'), 
      filters = 'affy_hg_u133_plus_2', 
      values = affyids, 
      mart = ensembl)

## ----cacheInfo------------------------------------------------------------------------------------
biomartCacheInfo()

## ----wormbase, echo=TRUE, eval=TRUE---------------------------------------------------------------
listMarts(host = "parasite.wormbase.org")
wormbase <- useMart(biomart = "parasite_mart", 
                   host = "https://parasite.wormbase.org", 
                   port = 443)

## ----wormbase-2, echo=TRUE, eval=TRUE-------------------------------------------------------------
listDatasets(wormbase)
wormbase <- useDataset(mart = wormbase, dataset = "wbps_gene")
head(listFilters(wormbase))
head(listAttributes(wormbase))
getBM(attributes = c("external_gene_id", "wbps_transcript_id", "transcript_biotype"), 
      filters = "gene_name", 
      values = c("unc-26","his-33"), 
      mart = wormbase)
     

## ---- phytozome-1, echo = TRUE, eval = TRUE-------------------------------------------------------
listMarts(host = "https://phytozome.jgi.doe.gov")
phytozome <- useMart(biomart = "phytozome_mart", 
                host = "https://phytozome.jgi.doe.gov")
listDatasets(phytozome)
wormbase <- useDataset(mart = phytozome, dataset = "phytozome")

## ---- pytozome-2----------------------------------------------------------------------------------
getBM(attributes = c("organism_name", "gene_name1"), 
      filters = "gene_name_filter", 
      values = "82092", 
      mart = phytozome)

## ---- phytozome-13, echo = TRUE, eval = TRUE------------------------------------------------------
phytozome_v13 <- useMart(biomart = "phytozome_mart", 
                dataset = "phytozome", 
                host = "https://phytozome-next.jgi.doe.gov")

## ----localCopy, eval = FALSE----------------------------------------------------------------------
#  listMarts(host="www.myLocalHost.org", path="/myPathToWebservice/martservice")
#  mart=useMart("nameOfMyMart",dataset="nameOfMyDataset",host="www.myLocalHost.org", path="/myPathToWebservice/martservice")

## ----columnsAndKeyTypes---------------------------------------------------------------------------
mart <- useMart(dataset="hsapiens_gene_ensembl",biomart='ensembl')
head(keytypes(mart), n=3)
head(columns(mart), n=3)

## ----keys1----------------------------------------------------------------------------------------
k = keys(mart, keytype="chromosome_name")
head(k, n=3)

## ----keys2----------------------------------------------------------------------------------------
k = keys(mart, keytype="chromosome_name", pattern="LRG")
head(k, n=3)

## ----select---------------------------------------------------------------------------------------
affy=c("202763_at","209310_s_at","207500_at")
select(mart, keys=affy, columns=c('affy_hg_u133_plus_2','entrezgene_id'),
  keytype='affy_hg_u133_plus_2')

## ----putenv, eval = FALSE-------------------------------------------------------------------------
#  Sys.setenv("http_proxy" = "http://my.proxy.org:9999")

## ----rCurlOptions, eval = FALSE-------------------------------------------------------------------
#  options(RCurlOptions = list(proxy="uscache.kcc.com:80",proxyuserpwd="------:-------"))

## ----sessionInfo----------------------------------------------------------------------------------
sessionInfo()
warnings()


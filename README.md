R-markdown_pajekfile
====================
---
title: "Markdown_from_downloadingdata_to_net_pajekfile"
author: "Sini Nagpal"
date: "Saturday, August 09, 2014"
output:
  html_document:
    toc: yes
---
###First install the GEOmetadb package:-
source("http://bioconductor.org/biocLite.R")
biocLite("GEOmetadb")

####Load the library
```{r}
library(GEOmetadb)
```
####Download complete metadb
getSQLiteFile(destdir = getwd(), destfile = "GEOmetadb.sqlite.gz")

####Set working directory
```{r}
setwd("E:/Work/new dataset for gpl80 platform")
```
####Connect with the sqlite file
```{r}
con <- dbConnect(SQLite(),'GEOmetadb.sqlite')
```
####Query to sqlite using specific keywords
```{r}
gds_subset=dbGetQuery(con,"SELECT * FROM gds WHERE sample_organism LIKE 
        '%Homo sapiens%' AND description LIKE '%lung%' AND sample_type LIKE '%RNA%'")
```
####Select gds id's
```{r}
gds_id<-gds_subset[,2]
```

####Download data corresponding to a specific gds_id
gds <- getGEO(gds_id[32])

####Extract metadata and expression data
metadata <- Columns(gds)
expressiondata <- Table(gds)
write.csv(metadata,"metadata32.csv")
write.csv(expressiondata,"exprsdata12.csv")


setwd("E:/GDS_IDS")
path = getwd()
for(i in 2:3){
  dir.create(gds_id[i])
  gds <- getGEO(gds_id[i])
  setwd(gds_id[i])
  metadata <- Columns(gds)
  expressiondata <- Table(gds)
  write.csv(metadata,"metadata.csv")
  write.csv(expressiondata, "exprsdata.csv")
 


 


###Normalize the expression data
```{r}
setwd("E:/Work/new dataset for gpl80 platform")
 mydata <- read.csv("metadata5.csv")
 ndt <- read.csv("exprsdatagpl80.csv")

boxplot(ndt)
log_transform<- log(ndt[,4:ncol(ndt)],base=2)

library(limma)
test <- normalizeQuantiles(log_transform)
boxplot(test)
test[is.na(test)] <- as.numeric(0)
```
####Split the groups from metadata to perform t test or anova
```{r}
rownames(test) <-ndt[,2]
col2 <- mydata[,3]
splitted_list <- split(mydata,col2)
```
###For t test get expression data corresponding to the 2 splitted groups
```{r}


if(nlevels(col2)=='2'){
  controlid<-splitted_list[[1]][,2]
               
               Control_ID <-colnames(test[,-1])[controlid]
              
               treatedid <-splitted_list[[2]][,2]
               
               treated_ID <-colnames(test[,-1])[treatedid]
               
Control_IDs <- colnames(test) %in% Control_ID
treated_IDs <- !colnames(test) %in% Control_ID 


results<-apply(test[,-1], 1, function(x){t.test(as.numeric(x[Control_IDs]), as.numeric(x[treated_IDs]))$p.value})

}
```
###Anova for more than 2 splitted groups
```{r}
if(nlevels(col2) > '2'){
new<- t(test)
splitted_groups <- mydata[,3]
required_format<- cbind(splitted_groups, new)
tempdf <- as.data.frame(required_format)
tempdf$splitted_groups <- as.factor(tempdf$splitted_groups)
results <- apply(tempdf[,-1], 2, function(x){
  oneway.test(x ~ tempdf$splitted_groups)$p.value
})
}
```

####Add a new column to your expressiondata containing p values
```{r}
testdata <- cbind(results,test)
```
####For diffrentially expressed genes(deg) p value is less than 0.05. Filter the results.               
```{r}
deg<-subset(testdata,testdata$results < 0.05)
nrow(deg)              
nrow(testdata)              
newdeg<-deg[,-(1:2)]
```

####Transpose of matrix containing expression values of diffrentially expressed genes               
```{r}             
tdeg<-t(newdeg)
```
###Obtain correlation matrix 
```{r}
rownames(tdeg) <- NULL       #deletes rownames
newtdeg <- apply(tdeg, 2, as.numeric) 
cormatrix <- cor(newtdeg)
absolute_value_cormatrix<- abs(cormatrix) # all values positive
```

###Obtain p values matrix from correlation matrix
```{r}
require(Hmisc)
pval<-rcorr(newtdeg)
class(pval)
sapply(pval,dim)

####Note: pval[[1]] contains correlation values, pval[[2]] containes dimensions of matrix and pval[[3]] contains p-values. Perform negative log transformation to obtain the final adjacency matrix.
```{r}
pvalue_matrix <- pval[[3]]
diag(pvalue_matrix) <- 1

```
####Take only significant p values and replace values>0.001 by zero 
```{r}
absolute_value_cormatrix[pvalue_matrix>0.005] <- 0 
```
####Load the necessary libraries
```{r}
require(reshape2)
```
###Create edge list
```{r}
Edgelist <- melt(absolute_value_cormatrix)
```
####Create edgelist in pajek format and save as text file
```{r}
Edgelist <-  subset(Edgelist, Edgelist$value != 0)

```
####Replace all punctuations with "_"
```{r}
Edgelist$Var1 <-gsub("[[:punct:]]", "_", Edgelist$Var1)
Edgelist$Var2 <-gsub("[[:punct:]]", "_", Edgelist$Var2)
```

###Save as .csv as well as .net file
```{r}
write.csv(Edgelist, "Edgelist.csv", row.names=FALSE)
```


# ---
#   title: "First Set of Exercises to Deliver"
# author: "Álvaro Ponce Cabrera"
# date: "April 04, 2016"
# output: pdf_document
# ---

# setwd("C:/Users/alvaro/Desktop")





#Exercise 1. 

generateExpressions<-function(nrow=5,ncol.1=1,mean.1=0,sd.1=1,samples.2=F,
                              ncol.2=NULL,mean.2=NULL,sd.2=NULL){
  out<-matrix(rnorm(nrow*ncol.1,mean.1,sd.1),nrow,ncol=ncol.1)
  if(samples.2){
    if(is.null(ncol.2)) ncol.2<-ncol.1
    if(is.null(mean.2)) mean.2<-mean.1
    if(is.null(sd.2)) sd.2<-sd.1
    matrix.2<-matrix(rnorm(nrow*ncol.2,mean.2,sd.2), nrow, ncol=ncol.2)
    out<-cbind(out,matrix.2)
  }else{
    ncol.2<- 0
  }
  colnames(out)<-paste("sample",1:(ncol.1+ncol.2),sep=".")
  rownames(out)<-paste("gene",1:nrow,sep=".")
  return(out)
}


##1.1 
generateExpressions(3, 2, 12, 1)
generateExpressions(2,5,1,1,samples.2=TRUE, ncol.2=3)





##1.2 

generateExpressions<-function(nrow=5, ncol.1=1, mean.1=0, sd.1=1, samples.2=F,
                              ncol.2=NULL, mean.2=NULL, sd.2=NULL)
{
  # Generate a matrix with nrow and ncol indicate in the firt and second arguments of the function
  # Elements = rnorm(nrow * ncol) (first argument of matrix, generate a normal distribution with
  #mean = mean.1 and sd = sd.1, arguments of the function)
  out<-matrix (rnorm (nrow*ncol.1, mean.1, sd.1), nrow, ncol=ncol.1)
  
  # Data for a second matrix if samples.2 = TRUE. 
  
  #If data of second matrix is null, the function use this data from matrix 1.
  if(samples.2){
    if(is.null(ncol.2)) ncol.2<-ncol.1
    if(is.null(mean.2)) mean.2<-mean.1
    if(is.null(sd.2)) sd.2<-sd.1
    
    #Creating a second matrix binded to the first one. 
    #This matrix is generated exactly like the first one but with arguments ended in ".2".
    matrix.2<-matrix (rnorm (nrow*ncol.2, mean.2, sd.2), nrow, ncol=ncol.2)
    out<-cbind (out, matrix.2)
  }
  
  else{
    ncol.2<- 0
  }
  
  #The function automaticaly name the columns and the rows like-> sample.numberofcolum and gene.numberofrow
  colnames(out)<-paste("sample",1:(ncol.1+ncol.2),sep=".")
  rownames(out)<-paste("gene",1:nrow,sep=".")
  return(out)
}

generateExpressions(6,3)

generateExpressions( )


##1.3 



#Exercise 2. 


##2.1 


gene <- read.table("gene_expr.data", header = TRUE)
head(gene)

##2.2 

length(gene$Probe)

##2.3

length(unique(gene$Symbol))



##2.4


(length(gene$Probe)) / (length(unique(gene$Symbol)))



##2.5 
numbers.f<- function (file) {
  
  gene.f <- read.table(file, header = TRUE)
  ProbeNumber <- length((gene.f$Probe))
  GeneNumber <- length(unique(gene.f$Symbol))
  ratio <- (length(gene.f$Probe)) / (length(unique(gene.f$Symbol)))
  return(c("ProbeNumber" = ProbeNumber,"GeneNumber" = GeneNumber, "Probe per Gene" = ratio))
  
}
numbers.f("gene_expr.data")


#Exercise 3. 


descriptive <- function (vector) {
  
  n<- length(vector)
  missing <- sum(is.na(vector))
  Mean <- mean(vector)
  Std.Dev<- sd(vector)
  Min <- min(vector)
  Q1 <- quantile(vector, 0.25, na.rm=TRUE)
  Median <- median(vector)
  Q3<- quantile(vector, 0.75,na.rm=TRUE)
  Max<- max(vector)
  
  return(c("number of measurements"= n, "NA's Number" = missing,"Mean" = Mean,
           "Standar deviation" = Std.Dev, "Minumun value" = Min, "Quantile 1" = Q1, "Median" = Median,
           "Quantile 3"= Q3, "Maximun value" = Max))
}

v<- c(rnorm(15,1,1))
descriptive(v)

#Exercise 4. 



exprDescriptive <- function (data = NULL, colpos = -(1:2), row.col) {
  data <- gene [,colpos]
  if (row.col== "row") {
    d.row <- c()
    for (i in 1:nrow(data)) {
      d.row<-c(d.row, descriptive(as.numeric(data[i,])))
    }
    out <- matrix (d.row, nrow = nrow (data), byrow = TRUE)
  }
  if (row.col == "col") {
    d.col <- c()
    for (i in 1:ncol(data)) {
      d.col<-c(d.col, descriptive(as.numeric(data[,i])))
    }
    
    out <- matrix (d.col, nrow= ncol(data), byrow = TRUE)
  }
  return(out)
}

head(exprDescriptive(gene, row.col  ="col"))


#Exercise 5. 

s<- "gccaattcgtatcgatctgacgt"
gc <- function (sequence) {
  c<-length(gregexpr("c",sequence)[[1]])
  g<-length(gregexpr("g",sequence)[[1]])
  out <- (g+c) / length(unlist(strsplit(sequence,"")))*100
  return(out)
}
gc(s)


#Exercise 6. 

wc <- function (file) {
  out <- length (readLines (file))
  return (out)
}
wc ("dna_seqs.fasta")

#Exercise 7. 

summarySingleFasta <- function (file) {
  my.sequence <- readLines (file)
  description <- my.sequence[1] 
  sequence <- my.sequence[2] 
  width <- length (unlist (strsplit (my.sequence[2], ""))) 
  gc.content <- gc (my.sequence[2])
  return (list ("Description" = description, "Sequence" = sequence, "Width" = width, "GC.content" = gc.content))
}
summarySingleFasta ("single_sequence.fasta")



#Exercise 8. 


summaryMultipleFasta <- function (file) {
  multi.fasta <- readLines (file)
  description <- strsplit(multi.fasta[seq(1,length(multi.fasta), by = 2)]," ")
  gene <- c(NULL)
  for (i in 1:length (description)) 
  { 
    gene <- c(gene, description[[i]][1]) 
  }
  specie <- c(NULL)
  
  for (i in 1:length (description)) { 
    specie <- c(specie, description[[i]][2]) 
  }
  
  sequence <- multi.fasta[seq(2,length(multi.fasta), by = 2)]
  width <- c(NULL)
  
  for (i in 1:length (sequence)) { 
    width <- c(width, length (unlist(strsplit(sequence[i], ""))))
  }
  gc.content <- c(NULL)
  
  for (i in 1:length (sequence)) { 
    gc.content <- c(gc.content, gc(sequence[i]))
  }
  
  out <- data.frame ("Gene" = gene, "Specie" = specie, "Sequence" = sequence, "Width" = width, "GC.content" = gc.content) 
  return(out)
}
head (summaryMultipleFasta ("dna_seqs.fasta"))


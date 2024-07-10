#!/usr/bin/env Rscript
#source activate R-Env
library(tidyverse)
library(stringr)
library(tidyr)
menangle_coverage2<- readr::read_tsv((commandArgs(TRUE)[1]),col_names = FALSE)
colnames(menangle_coverage2) <- c("Virus", "Nucleotide","Coverage")

genomelength <- read_tsv((commandArgs(TRUE)[2]),col_names = FALSE)

colnames(genomelength) <- c("genomelength")

menangle_coverage2 <- menangle_coverage2%>%
  slice(rep(1:nrow(menangle_coverage2), times = menangle_coverage2$Coverage))

annotated <- (commandArgs(TRUE)[5])
if(annotated==TRUE){
  geneannotations <- read_tsv((commandArgs(TRUE)[6]))
  colnames(geneannotations) <- c("Virus")
  geneannotations[ , 'gene']=NA
  geneannotations
  for (x in 1:nrow(geneannotations)){
    if(str_sub(geneannotations[x,1],-4,-1)=="gene"){
      geneannotations[x,2] <- str_sub(str_replace(geneannotations[x+1,1], '(.*?)\t(.*?)', ''))
    }else{
    }
  }
  geneannotations <- geneannotations[!is.na(geneannotations$gene),]
  geneannotations <- geneannotations %>% separate_wider_delim(Virus, delim = "\t", names = c("start", "end","annotation"))
  geneannotations$start <- as.numeric(geneannotations$start)
  geneannotations$end <- as.numeric(geneannotations$end)
  annotations <- geneannotations
}


#annotations<- data.frame(list(gene=c("N","V/P","F","L","HN","M"),
 #                                 start=c(56,1861,4818,8634,6691,3289),
  #                                end=c(1842,3283,6556,15485,8465,4714)))
  if(exists("annotations")) {
    h <- ggplot()+
      geom_histogram(data=menangle_coverage2,mapping=ggplot2::aes(x=Nucleotide),fill="black",binwidth = 1)+
      scale_x_continuous(limits = c(0, genomelength$genomelength),expand=c(0,0))+
      geom_rect(data=annotations,mapping=aes(xmin=start,xmax=end,ymin=-.5,ymax=-8,fill=gene),alpha=.5)+
        geom_text(data=geneannotations,mapping=aes(x=(start+end)/2,y=-4,label=gene))+
     theme(text=element_text(size=14,family="sans",color="Black",face="bold"))+
    theme_classic()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    scale_y_continuous(expand = c(0,0),name="Coverage")+
    theme(axis.text = element_text(colour = "black"))
  }else{
    h <- ggplot()+
      geom_histogram(data=menangle_coverage2,mapping=ggplot2::aes(x=Nucleotide),fill="black",binwidth = 1)+
      scale_x_continuous(limits = c(0, genomelength$genomelength),expand=c(0,0))+
      scale_y_continuous(expand = c(0,0),name="Coverage")+
      theme(text=element_text(size=14,family="sans",color="Black",face="bold"))+
      theme_classic()+
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
      scale_y_continuous(expand = c(0,0),name="Coverage")+
      theme(axis.text = element_text(colour = "black"))
  }
  
#h <- ggplot()+
#  geom_histogram(data=menangle_coverage2,mapping=aes(x=Nucleotide),fill="black",binwidth = 1)#+
 # scale_x_continuous(limits = c(0, genomelength$genomelength),expand=c(0,0))+
 # scale_y_continuous(expand = c(0,0),name="Coverage")

menangle_readlength<- readr::read_tsv((commandArgs(TRUE)[3]),col_names = FALSE)
colnames(menangle_readlength) <- c("ReadLength")
l <- ggplot(menangle_readlength, aes(x = "", y = ReadLength)) + 
  geom_violin(trim=FALSE,fill="gray")+
  theme_classic()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

pdf(file = (commandArgs(TRUE)[4])) # The height of the plot in inches


plot(h)
plot(l)

dev.off()

library(plotly)

p <- ggplot() +geom_histogram(data=menangle_coverage2,aes(x=Nucleotide, y=after_stat(count),weight=Coverage),binwidth=1)
p <- ggplotly(p)

htmlwidgets::saveWidget(as_widget(p), "/Users/greenebp/Desktop/Sequencing2/Menangle/Targeted/fastq_pass/barcode05/sequenceresults/coverage.html")
getwd()

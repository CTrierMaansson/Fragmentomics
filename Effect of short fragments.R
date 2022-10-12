setwd("C:/Users/Christoffer/OneDrive/1PhD/R files/Important excel and text documents")
gr <- function(x){
    library(GenomicRanges)
    y <- read.table(x, header = T)
    y$start <- as.numeric(y$start) 
    y$end <- as.numeric(y$end) 
    if ("Regions" %in% colnames(y)){
      y$Regions <- as.numeric(y$Regions) 
    }
    granges <- GRanges(y$seqnames, 
                       IRanges(start = y$start, end = y$end),
                       strand = y$strand)
    if ("gene_start" %in% colnames(y)){
      y$gene_start <- as.numeric(y$gene_start)
      y$gene_end <- as.numeric(y$gene_end)
    }
    mcols(granges)$SYMBOL <- y$SYMBOL 
    return(granges)
}
grs <- gr("AVENIO_genes.txt")


library(ggplot2)
th <- theme(
  legend.position = 'right',
  legend.background = element_rect(),
  plot.title = element_text(angle = 0, size = 16, face = 'bold', vjust = 1),
  plot.caption = element_text(angle = 0, size = 10, vjust = 1, hjust = 0.37),
  axis.text.x = element_text(angle = 0, size = 12),
  axis.text.y = element_text(angle = 0, size = 12),
  axis.title = element_text(size = 14, face = 'bold'),
  axis.title.x = element_text(size = 14, face = 'bold'),
  axis.title.y = element_text(size = 14, face = 'bold'),
  axis.line = element_line(colour = 'black'),
  legend.text = element_text(size = 12), 
  title = element_text(size = 14, face = "bold"))

####Fragment lengths####
x #Vector of names for BAM files for group 1
y #Vector of names for BAM files for group 2
z #vector of names for group 1 and 2
p #Title of plot
fragment_length <- function(x,y,z,p){
  library(Rsamtools)
  library(IRanges)
  library(dplyr)
  library(ggplot2)
  library(plyr)
  what <- c("isize")
  param <- ScanBamParam(which = grs, what = what)
  bam1 <- list()
  bam2 <- list()
  for (i in 1:length(x)){
    bam <- scanBam(x[i], param = param)
    bam_lengths <- unname(unlist(bam))
    bam_lengths <- bam_lengths[bam_lengths>0]
    bam1[[i]] <- bam_lengths
    print(paste("progress =",round(i/(length(x)+length(y))*100,2),"%"))
  }
  for (i in 1:length(y)){
    bam <- scanBam(y[i], param = param)
    bam_lengths <- unname(unlist(bam))
    bam_lengths <- bam_lengths[bam_lengths>0]
    bam2[[i]] <- bam_lengths
    print(paste("progress =",round((length(x)+i)/(length(x)+length(y))*100,2),"%"))
  }
  bam1 <- unlist(bam1)
  bam2 <- unlist(bam2)
  df <- data.frame(len = c(bam1,bam2),
                   Sample = c(rep(z[1],length(bam1)),rep(z[2],length(bam2))))
  mu <- ddply(df, "Sample", summarise, grp.median=median(len))
  print("Analysis done, printing plot")
  gg <- ggplot(data = df, aes(x = len, color = Sample))+
    geom_density(size = 1)+
    geom_vline(data=mu, aes(xintercept=grp.median, color=Sample),
               linetype="dashed", size = 1)+
    theme_bw(base_size = 15)+
    scale_color_manual(values = c("#6a00fc","#ffa10c"))+
    labs(title = p, y = "Density", x = "Fragment length (bp)")+
    scale_x_continuous(breaks = seq(0,400,50), limits = c(0,400))+
    th
  return(gg)
}
setwd("C:/Users/Christoffer/OneDrive/1PhD/Adeno, plano and SCLC article/cfChIP-seq/PosDeduped")
fragment_length(c("PosDeduped-A-1279-cfChIP.bam",
                  "PosDeduped-B-1288-cfChIP.bam",
                  "PosDeduped-C-1475-cfChIP.bam",
                  "PosDeduped-D-1578-cfChIP.bam",
                  "PosDeduped-E-439-cfChIP.bam",
                  "PosDeduped-F-1449-cfChIP.bam",
                  "PosDeduped-I-645-cfChIP.bam",
                  "PosDeduped-J-1663-cfChIP.bam",
                  "PosDeduped-G-514-cfChIP.bam",
                  "PosDeduped-H-1169-cfChIP.bam",
                  "PosDeduped-K-440-cfChIP.bam",
                  "PosDeduped-L-1100-cfChIP.bam"),
                c("PosDeduped-A-1279-input.bam",
                  "PosDeduped-B-1288-input.bam",
                  "PosDeduped-C-1475-input.bam",
                  "PosDeduped-D-1578-input.bam",
                  "PosDeduped-E-439-input.bam",
                  "PosDeduped-F-1449-input.bam",
                  "PosDeduped-I-645-input.bam",
                  "PosDeduped-J-1663-input.bam",
                  "PosDeduped-G-514-input.bam",
                  "PosDeduped-H-1169-input.bam",
                  "PosDeduped-K-440-input.bam",
                  "PosDeduped-L-1100-input.bam"),
                c("cfChIP", "Input"), 
                "Fragment lengths input and cfChIP cancer samples")


fragment_length(c("PosDeduped-B-1288-input.bam",
                  "PosDeduped-I-645-input.bam",
                  "PosDeduped-J-1663-input.bam"),
                c("PosDeduped-A-1279-input.bam",
                  "PosDeduped-C-1475-input.bam",
                  "PosDeduped-D-1578-input.bam",
                  "PosDeduped-E-439-input.bam",
                  "PosDeduped-F-1449-input.bam",
                  "PosDeduped-G-514-input.bam",
                  "PosDeduped-H-1169-input.bam",
                  "PosDeduped-K-440-input.bam",
                  "PosDeduped-L-1100-input.bam"),
                c("ctDNA negative", "ctDNA positive"), 
                "Fragment lengths ctDNA positive/negative cancer input samples")

gc()

x #Table consisting of BAM file names, genes mutated and the position of the mutation
y #Number of the tumor sample to plot
fragment_length_mut <- function(x,y=NULL){
  library(Rsamtools)
  library(IRanges)
  library(dplyr)
  library(ggplot2)
  library(plyr)
  `%ni%` <- Negate(`%in%`)
  what <- c("isize","pos","cigar","seq")
  df <- read.table(x,header = T)
  df <- df[!is.na(df$genes),]
  df6 <- data.frame(len = c(),
                    Sample = c())
  n_samples <- length(df$file)
  tot_genes <- 0
  for (i in 1:n_samples){
    genes <- strsplit(df$genes[i],split = ",")[[1]]
    which <- grs[elementMetadata(grs)[,1] %in% genes]
    index1 <- match(genes,elementMetadata(which)[,1])
    which <- which[index1]
    n_genes <- length(genes)
    tot_genes <- tot_genes + n_genes
    positions <- strsplit(df$positions[i],split = ",")[[1]]
    string <- strsplit(positions,split = ":")[1:n_genes]
    string <- string
    types <- strsplit(df$Type[i],split = ",")[[1]]
    types <- types
    chrs <- c()
    poss <- c()
    param <- ScanBamParam(which = which, what = what, tag = "MD")
    bam <- scanBam(df$file[i], param = param)
    name <- strsplit(names(bam),split = ":")
    spans <- unlist(lapply(name,"[[",2))
    sort_df <- data.frame(name = unlist(lapply(name,"[[",1)))
    starts <- c()
    ends <- c()
    for(j in 1:n_genes){
      span <- strsplit(spans[[j]],split = "-")
      span[[1]][1] <- as.numeric(span[[1]][1])
      span[[1]][2] <- as.numeric(span[[1]][2])
      starts[j] <- span[[1]][1]
      ends[j] <- span[[1]][2]
    }
    sort_df$start <- starts
    sort_df$end <- ends
    input <- as.data.frame(ranges(which))
    index <- match(input$start,sort_df$start)
    bam <- bam[index]
    df4 <- data.frame(len = c(),
                      Sample = c())
    for(j in 1:n_genes){
      chrs[j] <- string[[j]][1]
      poss[j] <- as.numeric(string[[j]][2])
      bam_gene <- bam[[j]]
      type <- types[j]
      if (type == "indel"){
        df1 <- data.frame(position = bam_gene$pos[(poss[j]-96) < bam_gene$pos & bam_gene$pos < poss[j]],
                          size = bam_gene$isize[(poss[j]-96) < bam_gene$pos & bam_gene$pos < poss[j]],
                          cigar = bam_gene$cigar[(poss[j]-96) < bam_gene$pos & bam_gene$pos < poss[j]])
        df1 <- df1[df1$size >0,]
        df1 <- df1[!is.na(df1$position),]
        df2 <- df1 %>% filter(grepl("D",cigar) | grepl("I",cigar))
        df3 <- df1 %>% filter(cigar %in% "96M")
        df4 <- rbind(df4,data.frame(len = c(df2$size,df3$size),
                         Sample = c(rep("ctDNA",length(df2$size)),
                                    rep("WT",length(df3$size)))))
      }
      else{
        df1 <- data.frame(position = bam_gene$pos[(poss[j]-96) < bam_gene$pos & bam_gene$pos < poss[j]],
                          size = bam_gene$isize[(poss[j]-96) < bam_gene$pos & bam_gene$pos < poss[j]],
                          miss = bam_gene$tag[[1]][(poss[j]-96) < bam_gene$pos & bam_gene$pos < poss[j]])
        df1 <- df1 %>% filter(position < poss[j])
        df1 <- df1[!is.na(df1$position),]
        df1 <- df1[df1$size >0,]
        df2 <- df1 %>% filter(miss %ni% "96")
        df3 <- df1 %>% filter(miss %in% "96")
        df4 <- rbind(df4,data.frame(len = c(df2$size,df3$size),
                                    Sample = c(rep("ctDNA",length(df2$size)),
                                               rep("WT",length(df3$size)))))
      }
    }  
    if (!is.null(y)){
      if (i == y){
        df4 <- na.omit(df4)
        mu <- ddply(df4, "Sample", summarise, grp.median=median(len))
        df_ctDNA <- df4 %>% filter(len < 150) %>% filter(Sample == "ctDNA")
        df_WT <- df4 %>% filter(len < 150) %>% filter(Sample == "WT")
        ct <- df4 %>% filter(Sample == "ctDNA")
        wt <- df4 %>% filter(Sample == "WT")
        under100 <- round(length(df_ctDNA$len)/(length(df_ctDNA$len)+length(df_WT$len))*100,2)
        gg <- ggplot(data = df4, aes(x = len, color = Sample))+
          geom_density(size = 1)+
          geom_vline(data=mu, aes(xintercept=grp.median, color=Sample),
                     linetype="dashed", size = 1)+
          theme_bw(base_size = 15)+
          scale_color_manual(values = c("#6a00fc","#ffa10c"))+
          labs(title = paste("ctDNA vs. WT for",df$file[i]), y = "Density", x = "Fragment length (bp)",
               subtitle = paste("ctDNA/WT fraction below 150 bp =",under100,"%"),
               caption = paste(length(ct$len),"ctDNA and",length(wt$len),
                               "WT reads in",df$genes[i]))+
          scale_x_continuous(breaks = seq(0,400,50), limits = c(0,400))+
          th+
          theme(plot.subtitle = element_text(size = 12, face = "bold"))
        return(gg)
      }
    }
    else{
      df6 <- rbind(df6,df4)
    }
  }
  df6 <- na.omit(df6)
  mu <- ddply(df6, "Sample", summarise, grp.median=median(len))
  df_ctDNA <- df6 %>% filter(len < 150) %>% filter(Sample == "ctDNA")
  df_WT <- df6 %>% filter(len < 150) %>% filter(Sample == "WT")
  ct <- df6 %>% filter(Sample == "ctDNA")
  wt <- df6 %>% filter(Sample == "WT")
  under100 <- round(length(df_ctDNA$len)/(length(df_ctDNA$len)+length(df_WT$len))*100,2)
  gg <- ggplot(data = df6, aes(x = len, color = Sample))+
    geom_density(size = 1)+
    geom_vline(data=mu, aes(xintercept=grp.median, color=Sample),
               linetype="dashed", size = 1)+
    theme_bw(base_size = 15)+
    scale_color_manual(values = c("#6a00fc","#ffa10c"))+
    labs(title = "ctDNA vs. WT", y = "Density", x = "Fragment length (bp)",
         subtitle = paste("ctDNA/WT fraction below 150 bp =",under100,"%"),
         caption = paste("samples =",n_samples,"      genes =",tot_genes,
                         "      ctDNA reads =",length(ct$len),
                         "      WT reads =",length(wt$len)))+
    scale_x_continuous(breaks = seq(0,400,50), limits = c(0,400))+
    th+
    theme(plot.subtitle = element_text(size = 12, face = "bold"))
  return(gg)
}
fragment_length_mut("cfChIP mutations.txt")
fragment_length_mut("BL mutations input files.txt")
setwd("C:/Users/Christoffer/OneDrive/1PhD/Adeno, plano and SCLC article/cfChIP-seq/PosDeduped")
fragment_length_mut("BL mutations.txt")
setwd("C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/PosDeduped files")

Single_fragment <- function(x){
  library(Rsamtools)
  library(IRanges)
  library(dplyr)
  library(ggplot2)
  library(plyr)
  `%ni%` <- Negate(`%in%`)
  what <- c("isize","pos","cigar","seq")
  df <- read.table(x,header = T)
  df <- df[!is.na(df$genes),]
  df6 <- data.frame(len = c(),
                    Sample = c())
  n_samples <- length(df$file)
  tot_genes <- 0
  myplots <- list()
  for (i in 1:n_samples){
    genes <- strsplit(df$genes[i],split = ",")[[1]]
    which <- grs[elementMetadata(grs)[,1] %in% genes]
    index1 <- match(genes,elementMetadata(which)[,1])
    which <- which[index1]
    n_genes <- length(genes)
    tot_genes <- tot_genes + n_genes
    positions <- strsplit(df$positions[i],split = ",")[[1]]
    string <- strsplit(positions,split = ":")[1:n_genes]
    string <- string
    types <- strsplit(df$Type[i],split = ",")[[1]]
    types <- types
    chrs <- c()
    poss <- c()
    param <- ScanBamParam(which = which, what = what, tag = "MD")
    bam <- scanBam(df$file[i], param = param)
    name <- strsplit(names(bam),split = ":")
    spans <- unlist(lapply(name,"[[",2))
    sort_df <- data.frame(name = unlist(lapply(name,"[[",1)))
    starts <- c()
    ends <- c()
    for(j in 1:n_genes){
      span <- strsplit(spans[[j]],split = "-")
      span[[1]][1] <- as.numeric(span[[1]][1])
      span[[1]][2] <- as.numeric(span[[1]][2])
      starts[j] <- span[[1]][1]
      ends[j] <- span[[1]][2]
    }
    sort_df$start <- starts
    sort_df$end <- ends
    input <- as.data.frame(ranges(which))
    index <- match(input$start,sort_df$start)
    bam <- bam[index]
    df4 <- data.frame(len = c(),
                      Sample = c())
    for(j in 1:n_genes){
      chrs[j] <- string[[j]][1]
      poss[j] <- as.numeric(string[[j]][2])
      bam_gene <- bam[[j]]
      type <- types[j]
      if (type == "indel"){
        df1 <- data.frame(position = bam_gene$pos[(poss[j]-96) < bam_gene$pos & bam_gene$pos < poss[j]],
                          size = bam_gene$isize[(poss[j]-96) < bam_gene$pos & bam_gene$pos < poss[j]],
                          cigar = bam_gene$cigar[(poss[j]-96) < bam_gene$pos & bam_gene$pos < poss[j]])
        df1 <- df1[df1$size >0,]
        df1 <- df1[!is.na(df1$position),]
        df2 <- df1 %>% filter(grepl("D",cigar) | grepl("I",cigar))
        df3 <- df1 %>% filter(cigar %in% "96M")
        df4 <- rbind(df4,data.frame(len = c(df2$size,df3$size),
                                    Sample = c(rep("ctDNA",length(df2$size)),
                                               rep("WT",length(df3$size)))))
      }
      else{
        df1 <- data.frame(position = bam_gene$pos[(poss[j]-96) < bam_gene$pos & bam_gene$pos < poss[j]],
                          size = bam_gene$isize[(poss[j]-96) < bam_gene$pos & bam_gene$pos < poss[j]],
                          miss = bam_gene$tag[[1]][(poss[j]-96) < bam_gene$pos & bam_gene$pos < poss[j]])
        df1 <- df1 %>% filter(position < poss[j])
        df1 <- df1[!is.na(df1$position),]
        df1 <- df1[df1$size >0,]
        df2 <- df1 %>% filter(miss %ni% "96")
        df3 <- df1 %>% filter(miss %in% "96")
        df4 <- rbind(df4,data.frame(len = c(df2$size,df3$size),
                                    Sample = c(rep("ctDNA",length(df2$size)),
                                               rep("WT",length(df3$size)))))
      }
    }  
        df4 <- na.omit(df4)
        mu <- ddply(df4, "Sample", summarise, grp.median=median(len))
        df_ctDNA <- df4 %>% filter(len < 150) %>% filter(Sample == "ctDNA")
        df_WT <- df4 %>% filter(len < 150) %>% filter(Sample == "WT")
        ct <- df4 %>% filter(Sample == "ctDNA")
        wt <- df4 %>% filter(Sample == "WT")
        under100 <- round(length(df_ctDNA$len)/(length(df_ctDNA$len)+length(df_WT$len))*100,2)
        p1 <- eval(substitute(
          ggplot(data = df4, aes(x = len, color = Sample))+
          geom_density(size = 1)+
          geom_vline(data=mu, aes(xintercept=grp.median, color=Sample),
                     linetype="dashed", size = 1)+
          theme_bw(base_size = 15)+
          scale_color_manual(values = c("#6a00fc","#ffa10c"))+
          labs(title = paste("ctDNA vs. WT for",df$file[i]), y = "Density", x = "Fragment length (bp)",
               subtitle = paste("ctDNA/WT fraction below 150 bp =",under100,"%"),
               caption = paste(length(ct$len),"ctDNA and",length(wt$len),
                               "WT reads in",df$genes[i]))+
          scale_x_continuous(breaks = seq(0,400,50), limits = c(0,400))+
          th+
          theme(plot.subtitle = element_text(size = 12, face = "bold")),
          list(i=i)))
        print(i)
        print(p1)
        myplots[[i]] <- p1
  }
  return(multiplot(plotlist = myplots, cols=3))
}
Single_fragment("cfChIP mutations.txt")
Single_fragment("BL mutations input files.txt")
x #Table consisting of BAM file names, genes mutated and the position of the mutation

MAF_under_length <- function(x,y=NULL){
  library(Rsamtools)
  library(IRanges)
  library(dplyr)
  library(ggplot2)
  library(plyr)
  `%ni%` <- Negate(`%in%`)
  what <- c("isize","pos","cigar","seq")
  df <- read.table(x,header = T)
  df <- df[!is.na(df$genes),]
  df6 <- data.frame(len = c(),
                    Sample = c())
  n_samples <- length(df$file)
  tot_genes <- 0
  intervals <- seq(50,400,10)
  df5 <- data.frame(percent = c(),n_fragments = c())
  for (i in 1:n_samples){
    genes <- strsplit(df$genes[i],split = ",")[[1]]
    which <- grs[elementMetadata(grs)[,1] %in% genes]
    index1 <- match(genes,elementMetadata(which)[,1])
    which <- which[index1]
    n_genes <- length(genes)
    tot_genes <- tot_genes + n_genes
    positions <- strsplit(df$positions[i],split = ",")[[1]]
    string <- strsplit(positions,split = ":")[1:n_genes]
    string <- string
    types <- strsplit(df$Type[i],split = ",")[[1]]
    types <- types
    chrs <- c()
    poss <- c()
    param <- ScanBamParam(which = which, what = what, tag = "MD")
    bam <- scanBam(df$file[i], param = param)
    name <- strsplit(names(bam),split = ":")
    spans <- unlist(lapply(name,"[[",2))
    sort_df <- data.frame(name = unlist(lapply(name,"[[",1)))
    starts <- c()
    ends <- c()
    for(j in 1:n_genes){
      span <- strsplit(spans[[j]],split = "-")
      span[[1]][1] <- as.numeric(span[[1]][1])
      span[[1]][2] <- as.numeric(span[[1]][2])
      starts[j] <- span[[1]][1]
      ends[j] <- span[[1]][2]
    }
    sort_df$start <- starts
    sort_df$end <- ends
    input <- as.data.frame(ranges(which))
    index <- match(input$start,sort_df$start)
    bam <- bam[index]
    df4 <- data.frame(len = c(),
                      Sample = c())
    for(j in 1:n_genes){
      chrs[j] <- string[[j]][1]
      poss[j] <- as.numeric(string[[j]][2])
      bam_gene <- bam[[j]]
      type <- types[j]
      if (type == "indel"){
        df1 <- data.frame(position = bam_gene$pos[(poss[j]-96) < bam_gene$pos & bam_gene$pos < poss[j]],
                          size = bam_gene$isize[(poss[j]-96) < bam_gene$pos & bam_gene$pos < poss[j]],
                          cigar = bam_gene$cigar[(poss[j]-96) < bam_gene$pos & bam_gene$pos < poss[j]])
        df1 <- df1[df1$size >0,]
        df1 <- df1[!is.na(df1$position),]
        df2 <- df1 %>% filter(grepl("D",cigar))
        df3 <- df1 %>% filter(cigar %in% "96M")
        df4 <- rbind(df4,data.frame(len = c(df2$size,df3$size),
                                    Sample = c(rep("ctDNA",length(df2$size)),
                                               rep("WT",length(df3$size)))))
      }
      else{
        df1 <- data.frame(position = bam_gene$pos[(poss[j]-96) < bam_gene$pos & bam_gene$pos < poss[j]],
                          size = bam_gene$isize[(poss[j]-96) < bam_gene$pos & bam_gene$pos < poss[j]],
                          miss = bam_gene$tag[[1]][(poss[j]-96) < bam_gene$pos & bam_gene$pos < poss[j]])
        df1 <- df1 %>% filter(position < poss[j])
        df1 <- df1[!is.na(df1$position),]
        df1 <- df1[df1$size >0,]
        df2 <- df1 %>% filter(miss %ni% "96")
        df3 <- df1 %>% filter(miss %in% "96")
        df4 <- rbind(df4,data.frame(len = c(df2$size,df3$size),
                                    Sample = c(rep("ctDNA",length(df2$size)),
                                               rep("WT",length(df3$size)))))
      }
    }  
    if (!is.null(y)){
      if (i == y){
        for (i in 1:length(intervals)){
          frag_len <- intervals[i]
          df_ctDNA <- df4 %>% filter(len < frag_len) %>% filter(Sample == "ctDNA")
          df_WT <- df4 %>% filter(len < frag_len) %>% filter(Sample == "WT")
          under <- length(df_ctDNA$len)/(length(df_ctDNA$len)+length(df_WT$len))*100
          n_fragments <-length(df_ctDNA$len)+length(df_WT$len)
          df5 <- rbind(df5,data.frame(percent = under,n_fragments = n_fragments))
        }
        df5$intervals <- intervals
        gg <- ggplot(data = df5, aes(x = intervals, y = percent, size = n_fragments))+
          geom_point(color = "#6a00fc")+
          scale_size_continuous(name = "Number of fragments")+
          theme_bw(base_size = 15)+
          labs(title = "MAF at fragment lengths below cutoffs", 
               y = "MAF %", x = "Fragment length < cutoff")+
          scale_y_continuous(breaks = seq(0,100,20), limits = c(0,100))+
          th
        return(gg)
      }
    }
    else{
      df6 <- rbind(df6,df4)
    }
  }  
  for (i in 1:length(intervals)){
    frag_len <- intervals[i]
    df_ctDNA <- df6 %>% filter(len < frag_len) %>% filter(Sample == "ctDNA")
    df_WT <- df6 %>% filter(len < frag_len) %>% filter(Sample == "WT")
    under <- length(df_ctDNA$len)/(length(df_ctDNA$len)+length(df_WT$len))*100
    n_fragments <-length(df_ctDNA$len)+length(df_WT$len)
    df5 <- rbind(df5,data.frame(percent = under,n_fragments = n_fragments))
  }
  df5$intervals <- intervals
  gg <- ggplot(data = df5, aes(x = intervals, y = percent, size = n_fragments))+
    geom_point(color = "#6a00fc")+
    scale_size_continuous(name = "Number of fragments")+
    theme_bw(base_size = 15)+
    labs(title = "MAF at fragment lengths below cutoffs", 
         y = "MAF %", x = "Fragment length < cutoff")+
    scale_y_continuous(breaks = seq(0,100,20), limits = c(0,100))+
    th
  return(gg)
}

MAF_under_length("BL mutations.txt")

MAF_over_length <- function(x){
  library(Rsamtools)
  library(IRanges)
  library(dplyr)
  library(ggplot2)
  library(plyr)
  `%ni%` <- Negate(`%in%`)
  what <- c("isize","pos","cigar")
  df <- read.table(x,header = T)
  df <- df[!is.na(df$genes),]
  df4 <- data.frame(len = c(),
                    Sample = c())
  n_samples <- length(df$file)
  tot_genes <- 0
  for (i in 1:n_samples){
    genes <- strsplit(df$genes[i],split = ",")[[1]]
    n_genes <- length(genes)
    tot_genes <- tot_genes + n_genes
    positions <- strsplit(df$positions[i],split = ",")[[1]]
    string <- strsplit(positions,split = ":")[1:n_genes]
    types <- strsplit(df$Type[i],split = ",")[[1]]
    chrs <- c()
    poss <- c()
    which <- grs[elementMetadata(grs)[,1] %in% genes]
    param <- ScanBamParam(which = which, what = what, tag = "MD")
    bam <- scanBam(df$file[i], param = param)
    print(df$file[i])
    for(j in 1:n_genes){
      chrs[j] <- string[[j]][1]
      poss[j] <- as.numeric(string[[j]][2])
      bam_gene <- bam[[j]]
      type <- types[j]
      if (type == "indel"){
        df1 <- data.frame(position = bam_gene$pos[(poss[j]-96) < bam_gene$pos & bam_gene$pos < poss[j]],
                          size = bam_gene$isize[(poss[j]-96) < bam_gene$pos & bam_gene$pos < poss[j]],
                          cigar = bam_gene$cigar[(poss[j]-96) < bam_gene$pos & bam_gene$pos < poss[j]])
        df1 <- df1[df1$size >0,]
        df2 <- df1 %>% filter(cigar %ni% "96M")
        df2 <- df2[!is.na(df2$cigar),]
        df3 <- df1 %>% filter(cigar %in% "96M")
        df4 <- rbind(df4,data.frame(len = c(df2$size,df3$size),
                                    Sample = c(rep("ctDNA",length(df2$size)),
                                               rep("WT",length(df3$size)))))
      }
      
      else{
        df1 <- data.frame(position = bam_gene$pos[(poss[j]-96) < bam_gene$pos & bam_gene$pos < poss[j]],
                          size = bam_gene$isize[(poss[j]-96) < bam_gene$pos & bam_gene$pos < poss[j]],
                          miss = bam_gene$tag[[1]][(poss[j]-96) < bam_gene$pos & bam_gene$pos < poss[j]])
        df1 <- df1[!is.na(df1$position),]
        df1 <- df1[df1$size >0,]
        df2 <- df1 %>% filter(miss %ni% "96")
        df3 <- df1 %>% filter(miss %in% "96")
        df4 <- rbind(df4,data.frame(len = c(df2$size,df3$size),
                                    Sample = c(rep("ctDNA",length(df2$size)),
                                               rep("WT",length(df3$size)))))
      }
    }
  }
  intervals <- seq(50,400,10)
  df5 <- data.frame(percent = c(),n_fragments = c())
  for (i in 1:length(intervals)){
    frag_len <- intervals[i]
    df_ctDNA <- df4 %>% filter(len > frag_len) %>% filter(Sample == "ctDNA")
    df_WT <- df4 %>% filter(len > frag_len) %>% filter(Sample == "WT")
    under <- length(df_ctDNA$len)/(length(df_ctDNA$len)+length(df_WT$len))*100
    n_fragments <-length(df_ctDNA$len)+length(df_WT$len)
    df5 <- rbind(df5,data.frame(percent = under,n_fragments = n_fragments))
  }
  df5$intervals <- intervals
  gg <- ggplot(data = df5, aes(x = intervals, y = percent, size = n_fragments))+
    geom_point(color = "#6a00fc")+
    scale_size_continuous(name = "Number of fragments")+
    theme_bw(base_size = 15)+
    labs(title = "MAF at fragment lengths above cutoffs", 
         y = "MAF %", x = "Fragment length > cutoff")+
    scale_y_continuous(breaks = seq(0,100,20), limits = c(0,100))+
    th
  return(gg)
}
MAF_over_length("BL mutations.txt")

####Effect of fragment size isolation on EGFR enrichment in Adeno1####
x #Name of BAM file
y #Granges object returned by gr() for the 197 AVENIO target genes
z #cutoff for maximum fragment length 
gene_count_low_length <- function(x,y,z){
  library(Rsamtools)
  library(dplyr)
  bf <- BamFile(x)
  which = y
  counts <- scanBam(bf, param = ScanBamParam(what = c("pos","isize"), which = which))
  reads <- c()
  name <- c()
  lengths <- c()
  for (i in 1:length(counts)){
    df <- data.frame(len = counts[[i]]$isize,
                     pos = counts[[i]]$pos)
    df <- df %>% filter(len>0) %>% filter(len<z)
    reads[i] <- length(df$pos)
    string <- names(counts)[i]
    splits <- strsplit(string, ":")
    chr <- splits[[1]][1]
    position <- as.numeric(strsplit(splits[[1]][2],"-")[[1]][1])
    g <- GRanges(chr, IRanges(start = position+10, end = position+110))
    gene <- which[which %over% g]
    name[i] <- gene$SYMBOL
  }
  df <- data.frame(genes = name, readcounts = reads)
  df <- df[order(df$genes),]
  return(df)
}
Adeno1_low <- gene_count_low_length("A-1279-cfChIP.bam",grs,150)
setwd("C:/Users/Christoffer/OneDrive/1PhD/Adeno, plano and SCLC article/cfChIP-seq")
Adeno1_low
x #Gene count table returned by gene_count()
y #.txt file of the AVENIO genes with the number of bases sequenced in each gene. based on data.frame returned by coverages()
z #Granges object returned by gr() for the 197 AVENIO target genes
h #data.frame object returned by healthy() for normilization to healthy. Otherwise NULL
e.score <- function(x,y,z,h = NULL){
  library(GenomicAlignments)
  library(GenomicRanges)
  library(BiocParallel)
  library(dplyr)
  y <- read.table(y, header=T)
  y$coverage <- as.numeric(y$coverage)
  ChIP_reads <- x
  ChIP_reads <- ChIP_reads %>% filter(readcounts > 10)
  if(is.null(h)){
  }
  else {
    h <- h %>% filter(genes %in% ChIP_reads$genes)
    h <- h[match(ChIP_reads$genes, h$genes),]
    ChIP_reads$readcounts <- ChIP_reads$readcounts - h$readcounts 
  }
  y <- y[match(ChIP_reads$genes, y$SYMBOL),]
  len <- sum(ChIP_reads$readcounts)
  RPKM <- c()
  for (i in 1:length(ChIP_reads$genes)){
    RPKM[i] <- (ChIP_reads$readcounts[i]*1000*1000000)/(len*y$coverage[i])
  }
  ChIP_reads$e <- RPKM
  res <- data.frame(genes = ChIP_reads$genes, enrichment=ChIP_reads$e)
  res <- res[order(res$enrichment),]
  return(res)
}

setwd("C:/Users/Christoffer/OneDrive/1PhD/R files/Important excel and text documents")
Adeno1_low_enrichment <- e.score(Adeno1_low, "Coverage of AVENIO genes.txt",gr("AVENIO_genes.txt"))
Adeno1_low_enrichment

x #name of BAM file
l #lower limit for fragment length cutoff
u #upper limit for fragment length cutoff
EGFR_enrichment_over_cutoff <- function(x,l,u){
  intervals <- seq(l,u,10)
  enrichment <- c()
  for (i in 1:length(intervals)){
    setwd("C:/Users/Christoffer/OneDrive/1PhD/Adeno, plano and SCLC article/cfChIP-seq")
    gene_counts <- gene_count_low_length(x,grs,intervals[i])
    setwd("C:/Users/Christoffer/OneDrive/1PhD/R files/Important excel and text documents")
    enrich <- e.score(gene_counts, "Coverage of AVENIO genes.txt",grs)
    enrich_EGFR <- enrich %>% filter(genes == "EGFR")
    enrichment[i] <- enrich_EGFR$enrichment[1]
    print(intervals[i])
    print(enrichment)
  }
  setwd("C:/Users/Christoffer/OneDrive/1PhD/Adeno, plano and SCLC article/cfChIP-seq")
  dfx <- data.frame(intervals = intervals,
                   enrichment = enrichment)
  gg <- ggplot(dfx, aes(x = intervals, y = enrichment))+
    geom_point()
  return(gg)
}
EGFR_enrichment_over_cutoff("D-1578-cfChIP.bam",100,150)

Adeno1_low <- gene_count_low_length("A-1279-cfChIP.bam",grs,150)
Adeno2_low <- gene_count_low_length("B-1288-cfChIP.bam",grs,150)
Adeno3_low <- gene_count_low_length("C-1475-cfChIP.bam", grs,150)
Adeno4_low <- gene_count_low_length("D-1578-cfChIP.bam", grs,150)
Plano1_low <- gene_count_low_length("E-439-cfChIP.bam", grs,150)
Plano2_low <- gene_count_low_length("F-1449-cfChIP.bam", grs,150)
Plano3_low <- gene_count_low_length("I-645-cfChIP.bam",grs,150)
Plano4_low <- gene_count_low_length("J-1663-cfChIP.bam",grs,150)
SCLC1_low <- gene_count_low_length("G-514-cfChIP.bam", grs,150)
SCLC2_low <- gene_count_low_length("H-1169-cfChIP.bam", grs,150)
SCLC3_low <- gene_count_low_length("K-440-cfChIP.bam", grs,150)
SCLC4_low <- gene_count_low_length("L-1100-cfChIP.bam", grs,150)



setwd("C:/Users/Christoffer/OneDrive/1PhD/R files/Important excel and text documents")
Adeno1_low_enrichment <- e.score(Adeno1_low, "Coverage of AVENIO genes.txt",gr("AVENIO_genes.txt"))
Adeno2_low_enrichment <- e.score(Adeno2_low, "Coverage of AVENIO genes.txt",gr("AVENIO_genes.txt"))
Adeno3_low_enrichment <- e.score(Adeno3, "Coverage of AVENIO genes.txt",gr("AVENIO_genes.txt"))
Adeno4_low_enrichment <- e.score(Adeno4, "Coverage of AVENIO genes.txt",gr("AVENIO_genes.txt"))
Plano1_low_enrichment <- e.score(Plano1, "Coverage of AVENIO genes.txt",gr("AVENIO_genes.txt"))
Plano2_low_enrichment <- e.score(Plano2, "Coverage of AVENIO genes.txt",gr("AVENIO_genes.txt"))
Plano3_low_enrichment <- e.score(Plano3, "Coverage of AVENIO genes.txt",gr("AVENIO_genes.txt"))
Plano4_low_enrichment <- e.score(Plano4, "Coverage of AVENIO genes.txt",gr("AVENIO_genes.txt"))
SCLC1_low_enrichment <- e.score(SCLC1, "Coverage of AVENIO genes.txt",gr("AVENIO_genes.txt"))
SCLC2_low_enrichment <- e.score(SCLC2, "Coverage of AVENIO genes.txt",gr("AVENIO_genes.txt"))
SCLC3_low_enrichment <- e.score(SCLC3, "Coverage of AVENIO genes.txt",gr("AVENIO_genes.txt"))
SCLC4_low_enrichment <- e.score(SCLC4, "Coverage of AVENIO genes.txt",gr("AVENIO_genes.txt"))

Adeno_tot_low <- list(Adeno1_low_enrichment,Adeno2_low_enrichment,Adeno3_low_enrichment,Adeno4_low_enrichment)
Plano_tot_low <- list(Plano1_low_enrichment,Plano2_low_enrichment,Plano3_low_enrichment,Plano4_low_enrichment)
SCLC_tot_low <- list(SCLC1_low_enrichment,SCLC2_low_enrichment,SCLC3_low_enrichment,SCLC4_low_enrichment)
NSCLC_tot_low <- list(Adeno1_low_enrichment,Adeno2_low_enrichment,Adeno3_low_enrichment,Adeno4_low_enrichment,
                  Plano1_low_enrichment,Plano2_low_enrichment,Plano3_low_enrichment,Plano4_low_enrichment)
Adeno_tot_enrichment_low <- average_enrichment(Adeno_tot_low)
Plano_tot_enrichment_low <- average_enrichment(Plano_tot_low)
SCLC_tot_enrichment_low <- average_enrichment(SCLC_tot_low)
NSCLC_tot_enrichment_low <- average_enrichment(NSCLC_tot_low)

ChIPcorr(Adeno_tot_enrichment_low, Plano_tot_enrichment_low, "Adeno avg low", "Plano avg low", bad)
ChIPcorr(NSCLC_tot_enrichment_low, SCLC_tot_enrichment_low, "NSCLC avg low", "SCLC avg low", bad)

ChIPcorr(Adeno1_low_enrichment, Adeno2_low_enrichment, "Adeno 1 low", "Adeno 2 low",bad)
ChIPcorr(Adeno1_low_enrichment, Adeno3_low_enrichment, "Adeno 1 low", "Adeno 3 low",bad)
ChIPcorr(Adeno1_low_enrichment, Adeno4_low_enrichment, "Adeno 1 low", "Adeno 4 low",bad)
ChIPcorr(Adeno2_low_enrichment, Adeno3_low_enrichment, "Adeno 2 low", "Adeno 3 low",bad)
ChIPcorr(Adeno2_low_enrichment, Adeno4_low_enrichment, "Adeno 2 low", "Adeno 4 low",bad)
ChIPcorr(Adeno3_low_enrichment, Adeno4_low_enrichment, "Adeno 3 low", "Adeno 4 low",bad)

NSCLC_tot_enrichment_low %>% filter(genes == "CACNA1E")
SCLC_tot_enrichment_low %>% filter(genes == "CACNA1E")

####Fragment lengths of top 15 and bottom 15 genes for an individual####

x #name of cfChIP BAM file
y #Granges object returned by gr() for the 197 AVENIO target genes
z #.txt file of the AVENIO genes with the number of bases sequenced in each gene. based on data.frame returned by coverages()
i #name of input BMA file
h #data.frame object returned by healthy() for normilization to healthy. Otherwise NULL
t #name of sample

fragment_active_vs_inactive <- function(x,y,z,i,h,t){
  setwd("C:/Users/Christoffer/OneDrive/1PhD/Adeno, plano and SCLC article/cfChIP-seq")
  c <- gene_count(x,y)
  setwd("C:/Users/Christoffer/OneDrive/1PhD/R files/Important excel and text documents")
  getwd()
  e <- e.score(c,z,y,h)
  top15 <- e$genes[(length(e$genes)-14):length(e$genes)]
  bottom15 <- e$genes[1:15]
  which <- grs[elementMetadata(y)[,1] %in% c(top15)]
  index1 <- match(top15,elementMetadata(which)[,1])
  which <- which[index1]
  what <- c("isize")
  param <- ScanBamParam(which = which, what = what)
  setwd("C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/PosDeduped files")
  bam <- scanBam(i, param = param)
  bam <- unlist(bam)
  bam <- bam[bam>0]
  bam <- bam[!is.na(bam)]
  df <- data.frame(len = bam, Sample = rep("Top15",length(bam)))
  which <- grs[elementMetadata(y)[,1] %in% c(bottom15)]
  index1 <- match(bottom15,elementMetadata(which)[,1])
  which <- which[index1]
  what <- c("isize")
  param <- ScanBamParam(which = which, what = what)
  bam <- scanBam(i, param = param)
  bam <- unlist(bam)
  bam <- bam[bam>0]
  bam <- bam[!is.na(bam)]
  df <- rbind(df,data.frame(len = bam, Sample = rep("Bottom15",length(bam))))
  p <- wilcox.test(df$len[df$Sample=="Top15"], df$len[df$Sample=="Bottom15"],
              var.equal = T, alternative = "two.sided")$p.value
  mu <- ddply(df, "Sample", summarise, grp.median=median(len))
  gg <- ggplot(data = df, aes(x = len, color = Sample))+
    geom_density(size = 1)+
    geom_vline(data=mu, aes(xintercept=grp.median, color=Sample),
               linetype="dashed", size = 1)+
    theme_bw(base_size = 15)+
    scale_color_manual(values = c("#6a00fc","#ffa10c"))+
    labs(title = paste("Top 15 vs. Bottom 15 for",t), y = "Density", x = "Fragment length (bp)",
         subtitle = paste("P-value =",round(p,5)))+
    scale_x_continuous(breaks = seq(0,400,50), limits = c(0,400))+
    th+
    theme(plot.subtitle = element_text(size = 12, face = "bold"))
  return(gg)
}
fragment_active_vs_inactive("A-1279-cfChIP.bam",
                            grs,
                            "Coverage of AVENIO genes.txt",
                            "PosDeduped-A-1279-input.bam",
                            healthy_reads,
                            "Adeno 1")
fragment_active_vs_inactive("B-1288-cfChIP.bam",
                            grs,
                            "Coverage of AVENIO genes.txt",
                            "PosDeduped-B-1288-input.bam",
                            healthy_reads,
                            "Adeno 2")
fragment_active_vs_inactive("C-1475-cfChIP.bam",
                            grs,
                            "Coverage of AVENIO genes.txt",
                            "PosDeduped-C-1475-input.bam",
                            healthy_reads,
                            "Adeno 3")
fragment_active_vs_inactive("D-1578-cfChIP.bam",
                            grs,
                            "Coverage of AVENIO genes.txt",
                            "PosDeduped-D-1578-input.bam",
                            healthy_reads,
                            "Adeno 4")
fragment_active_vs_inactive("E-439-cfChIP.bam",
                            grs,
                            "Coverage of AVENIO genes.txt",
                            "PosDeduped-E-439-input.bam",
                            healthy_reads,
                            "Plano 1")
fragment_active_vs_inactive("F-1449-cfChIP.bam",
                            grs,
                            "Coverage of AVENIO genes.txt",
                            "PosDeduped-F-1449-input.bam",
                            healthy_reads,
                            "Plano 2")
fragment_active_vs_inactive("I-645-cfChIP.bam",
                            grs,
                            "Coverage of AVENIO genes.txt",
                            "PosDeduped-I-645-input.bam",
                            healthy_reads,
                            "Plano 3")
fragment_active_vs_inactive("J-1663-cfChIP.bam",
                            grs,
                            "Coverage of AVENIO genes.txt",
                            "PosDeduped-J-1663-input.bam",
                            healthy_reads,
                            "Plano 4")
fragment_active_vs_inactive("G-514-cfChIP.bam",
                            grs,
                            "Coverage of AVENIO genes.txt",
                            "PosDeduped-G-514-input.bam",
                            healthy_reads,
                            "SCLC 1")
fragment_active_vs_inactive("H-1169-cfChIP.bam",
                            grs,
                            "Coverage of AVENIO genes.txt",
                            "PosDeduped-H-1169-input.bam",
                            healthy_reads,
                            "SCLC 2")
fragment_active_vs_inactive("K-440-cfChIP.bam",
                            grs,
                            "Coverage of AVENIO genes.txt",
                            "PosDeduped-K-440-input.bam",
                            healthy_reads,
                            "SCLC 3")
fragment_active_vs_inactive("L-1100-cfChIP.bam",
                            grs,
                            "Coverage of AVENIO genes.txt",
                            "PosDeduped-K-440-input.bam",
                            healthy_reads,
                            "SCLC 4")
gc()


x # name of BAM file
proportion_sub150 <- function(x){
  library(Rsamtools)
  what <- c("isize")
  param <- ScanBamParam(which = grs, what = what)
  bam <- scanBam(x, param = param)
  bam_lengths <- unname(unlist(bam))
  bam_lengths <- bam_lengths[bam_lengths>0]
  bam_lengths <- bam_lengths[!is.na(bam_lengths)]
  bam_lengths_sub150 <- bam_lengths[bam_lengths<150]
  res <- length(bam_lengths_sub150)/length(bam_lengths)
  return(res)
}

setwd("C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/PosDeduped files")
proportion_sub150("PosDeduped-B12-311.bam")
proportion_sub150("PosDeduped-A11-113.bam")
proportion_sub150("PosDeduped-A41-103.bam")
proportion_sub150("PosDeduped-A11-115.bam")
proportion_sub150("PosDeduped-A12-113.bam")
proportion_sub150("PosDeduped-A12-111.bam")
proportion_sub150("PosDeduped-B12-315.bam")
proportion_sub150("PosDeduped-A11-116.bam")
proportion_sub150("PosDeduped-B12-322.bam")
proportion_sub150("PosDeduped-B21-303.bam")

proportion_sub150("PosDeduped-B11-304.bam")
proportion_sub150("PosDeduped-A11-104.bam")
proportion_sub150("PosDeduped-A12-101.bam")
proportion_sub150("PosDeduped-A11-106.bam")
proportion_sub150("PosDeduped-B41-308.bam")
proportion_sub150("PosDeduped-A12-115.bam")
proportion_sub150("PosDeduped-A21-102.bam")

proportion_sub150("PosDeduped-Rask_kontrol_1_input.bam")
proportion_sub150("PosDeduped-Rask_kontrol_2_input.bam")
proportion_sub150("PosDeduped-Rask_kontrol_3_input.bam")
proportion_sub150("PosDeduped-Rask_kontrol_4_input.bam")
proportion_sub150("PosDeduped-Rask_kontrol_5_input.bam")
proportion_sub150("PosDeduped-Rask_kontrol_6_input.bam")
proportion_sub150("PosDeduped-Rask_kontrol_7_input.bam")

proportion_sub150("PosDeduped-Rask_kontrol_1_cfChIP.bam")
proportion_sub150("PosDeduped-Rask_kontrol_2_cfChIP.bam")
proportion_sub150("PosDeduped-Rask_kontrol_3_cfChIP.bam")
proportion_sub150("PosDeduped-Rask_kontrol_4_cfChIP.bam")

setwd("C:/Users/Christoffer/OneDrive/1PhD/Adeno, plano and SCLC article/cfChIP-seq/PosDeduped")

proportion_sub150("PosDeduped-A-1279-cfChIP.bam")
proportion_sub150("PosDeduped-B-1288-cfChIP.bam")
proportion_sub150("PosDeduped-C-1475-cfChIP.bam")
proportion_sub150("PosDeduped-D-1578-cfChIP.bam")
proportion_sub150("PosDeduped-E-439-cfChIP.bam")
proportion_sub150("PosDeduped-F-1449-cfChIP.bam")
proportion_sub150("PosDeduped-I-645-cfChIP.bam")
proportion_sub150("PosDeduped-J-1663-cfChIP.bam")
proportion_sub150("PosDeduped-G-514-cfChIP.bam")
proportion_sub150("PosDeduped-H-1169-cfChIP.bam")
proportion_sub150("PosDeduped-K-440-cfChIP.bam")
proportion_sub150("PosDeduped-L-1100-cfChIP.bam")

proportion_sub150("PosDeduped-A-1279-input.bam")
proportion_sub150("PosDeduped-B-1288-input.bam")
proportion_sub150("PosDeduped-C-1475-input.bam")
proportion_sub150("PosDeduped-D-1578-input.bam")
proportion_sub150("PosDeduped-E-439-input.bam")
proportion_sub150("PosDeduped-F-1449-input.bam")
proportion_sub150("PosDeduped-I-645-input.bam")
proportion_sub150("PosDeduped-J-1663-input.bam")
proportion_sub150("PosDeduped-G-514-input.bam")
proportion_sub150("PosDeduped-H-1169-input.bam")
proportion_sub150("PosDeduped-K-440-input.bam")
proportion_sub150("PosDeduped-L-1100-input.bam")

proportion_sub150_mut <- function(x,y=NULL){
  library(Rsamtools)
  library(IRanges)
  library(dplyr)
  library(ggplot2)
  library(plyr)
  `%ni%` <- Negate(`%in%`)
  what <- c("isize","pos","cigar","seq")
  df <- read.table(x,header = T)
  df <- df[!is.na(df$genes),]
  df6 <- data.frame(len = c(),
                    Sample = c())
  n_samples <- length(df$file)
  tot_genes <- 0
  for (i in 1:n_samples){
    genes <- strsplit(df$genes[i],split = ",")[[1]]
    which <- grs[elementMetadata(grs)[,1] %in% genes]
    index1 <- match(genes,elementMetadata(which)[,1])
    which <- which[index1]
    n_genes <- length(genes)
    tot_genes <- tot_genes + n_genes
    positions <- strsplit(df$positions[i],split = ",")[[1]]
    string <- strsplit(positions,split = ":")[1:n_genes]
    string <- string
    types <- strsplit(df$Type[i],split = ",")[[1]]
    types <- types
    chrs <- c()
    poss <- c()
    param <- ScanBamParam(which = which, what = what, tag = "MD")
    bam <- scanBam(df$file[i], param = param)
    name <- strsplit(names(bam),split = ":")
    spans <- unlist(lapply(name,"[[",2))
    sort_df <- data.frame(name = unlist(lapply(name,"[[",1)))
    starts <- c()
    ends <- c()
    for(j in 1:n_genes){
      span <- strsplit(spans[[j]],split = "-")
      span[[1]][1] <- as.numeric(span[[1]][1])
      span[[1]][2] <- as.numeric(span[[1]][2])
      starts[j] <- span[[1]][1]
      ends[j] <- span[[1]][2]
    }
    sort_df$start <- starts
    sort_df$end <- ends
    input <- as.data.frame(ranges(which))
    index <- match(input$start,sort_df$start)
    bam <- bam[index]
    df4 <- data.frame(len = c(),
                      Sample = c())
    for(j in 1:n_genes){
      chrs[j] <- string[[j]][1]
      poss[j] <- as.numeric(string[[j]][2])
      bam_gene <- bam[[j]]
      type <- types[j]
      if (type == "indel"){
        df1 <- data.frame(position = bam_gene$pos[(poss[j]-96) < bam_gene$pos & bam_gene$pos < poss[j]],
                          size = bam_gene$isize[(poss[j]-96) < bam_gene$pos & bam_gene$pos < poss[j]],
                          cigar = bam_gene$cigar[(poss[j]-96) < bam_gene$pos & bam_gene$pos < poss[j]])
        df1 <- df1[df1$size >0,]
        df1 <- df1[!is.na(df1$position),]
        df2 <- df1 %>% filter(grepl("D",cigar) | grepl("I",cigar))
        df3 <- df1 %>% filter(cigar %in% "96M")
        df4 <- rbind(df4,data.frame(len = c(df2$size,df3$size),
                                    Sample = c(rep("ctDNA",length(df2$size)),
                                               rep("WT",length(df3$size)))))
      }
      else{
        df1 <- data.frame(position = bam_gene$pos[(poss[j]-96) < bam_gene$pos & bam_gene$pos < poss[j]],
                          size = bam_gene$isize[(poss[j]-96) < bam_gene$pos & bam_gene$pos < poss[j]],
                          miss = bam_gene$tag[[1]][(poss[j]-96) < bam_gene$pos & bam_gene$pos < poss[j]])
        df1 <- df1 %>% filter(position < poss[j])
        df1 <- df1[!is.na(df1$position),]
        df1 <- df1[df1$size >0,]
        df2 <- df1 %>% filter(miss %ni% "96")
        df3 <- df1 %>% filter(miss %in% "96")
        df4 <- rbind(df4,data.frame(len = c(df2$size,df3$size),
                                    Sample = c(rep("ctDNA",length(df2$size)),
                                               rep("WT",length(df3$size)))))
      }
    }  
    if (!is.null(y)){
      if (i == y){
        df4 <- na.omit(df4)
        df_ctDNA <- df4 %>% filter(len < 150) %>% filter(Sample == "ctDNA")
        df_WT <- df4 %>% filter(len < 150) %>% filter(Sample == "WT")
        ct <- df4 %>% filter(Sample == "ctDNA")
        wt <- df4 %>% filter(Sample == "WT")
        ct_150 <- length(df_ctDNA$len)/length(ct$len)
        WT_150 <- length(df_WT$len)/length(wt$len)
        return(list(WT = WT_150,
               ctDNA = ct_150,
               name = df$file[i]))
      }
    }
  }
}
proportion_sub150_mut("cfChIP mutations.txt")
setwd("C:/Users/Christoffer/OneDrive/1PhD/Adeno, plano and SCLC article/cfChIP-seq/PosDeduped")


proportion_sub150_mut("BL mutations input files.txt",1)
proportion_sub150_mut("BL mutations input files.txt",2)
proportion_sub150_mut("BL mutations input files.txt",3)
proportion_sub150_mut("BL mutations input files.txt",4)
proportion_sub150_mut("BL mutations input files.txt",5)
proportion_sub150_mut("BL mutations input files.txt",6)
proportion_sub150_mut("BL mutations input files.txt",7)
proportion_sub150_mut("BL mutations input files.txt",8)
proportion_sub150_mut("BL mutations input files.txt",9)



setwd("C:/Users/Christoffer/OneDrive/1PhD/Adeno, plano and SCLC article/cfChIP-seq/PosDeduped")
proportion_sub150_mut("cfChIP mutations.txt",1)
proportion_sub150_mut("cfChIP mutations.txt",2)
proportion_sub150_mut("cfChIP mutations.txt",3)
proportion_sub150_mut("cfChIP mutations.txt",4)
proportion_sub150_mut("cfChIP mutations.txt",5)
proportion_sub150_mut("cfChIP mutations.txt",6)
proportion_sub150_mut("cfChIP mutations.txt",7)
proportion_sub150_mut("cfChIP mutations.txt",8)
proportion_sub150_mut("cfChIP mutations.txt",9)
proportion_sub150_mut("cfChIP mutations.txt",10)
proportion_sub150_mut("cfChIP mutations.txt",11)
proportion_sub150_mut("cfChIP mutations.txt",12)



x #name of cfChIP BAM file
y #Granges object returned by gr() for the 197 AVENIO target genes
z #.txt file of the AVENIO genes with the number of bases sequenced in each gene. based on data.frame returned by coverages()
i #name of input BMA file
h #data.frame object returned by healthy() for normilization to healthy. Otherwise NULL
t #name of sample



fragment_active_vs_inactive_sub150 <- function(x,y,z,i,h,t){
  setwd("C:/Users/Christoffer/OneDrive/1PhD/Adeno, plano and SCLC article/cfChIP-seq")
  c <- gene_count(x,y)
  setwd("C:/Users/Christoffer/OneDrive/1PhD/R files/Important excel and text documents")
  getwd()
  e <- e.score(c,z,y,h)
  top15 <- e$genes[(length(e$genes)-14):length(e$genes)]
  bottom15 <- e$genes[1:15]
  which <- grs[elementMetadata(y)[,1] %in% c(top15)]
  index1 <- match(top15,elementMetadata(which)[,1])
  which <- which[index1]
  what <- c("isize")
  param <- ScanBamParam(which = which, what = what)
  setwd("C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/PosDeduped files")
  bam <- scanBam(i, param = param)
  bam <- unlist(bam)
  bam <- bam[bam>0]
  bam <- bam[!is.na(bam)]
  df <- data.frame(len = bam, Sample = rep("Top15",length(bam)))
  which <- grs[elementMetadata(y)[,1] %in% c(bottom15)]
  index1 <- match(bottom15,elementMetadata(which)[,1])
  which <- which[index1]
  what <- c("isize")
  param <- ScanBamParam(which = which, what = what)
  bam <- scanBam(i, param = param)
  bam <- unlist(bam)
  bam <- bam[bam>0]
  bam <- bam[!is.na(bam)]
  df <- rbind(df,data.frame(len = bam, Sample = rep("Bottom15",length(bam))))
  df_top <- df %>% filter(len < 150) %>% filter(Sample == "Top15")
  df_bottom <- df %>% filter(len < 150) %>% filter(Sample == "Bottom15")
  top <- df %>% filter(Sample == "Top15")
  bottom <- df %>% filter(Sample == "Bottom15")
  top_150 <- length(df_top$len)/length(top$len)
  bottom_150 <- length(df_bottom$len)/length(bottom$len)
  return(list(Bottom = bottom_150,
              Top = top_150))
}
fragment_active_vs_inactive_sub150("A-1279-cfChIP.bam",
                            grs,
                            "Coverage of AVENIO genes.txt",
                            "PosDeduped-A-1279-input.bam",
                            healthy_reads,
                            "Adeno 1")
fragment_active_vs_inactive_sub150("B-1288-cfChIP.bam",
                            grs,
                            "Coverage of AVENIO genes.txt",
                            "PosDeduped-B-1288-input.bam",
                            healthy_reads,
                            "Adeno 2")
fragment_active_vs_inactive_sub150("C-1475-cfChIP.bam",
                            grs,
                            "Coverage of AVENIO genes.txt",
                            "PosDeduped-C-1475-input.bam",
                            healthy_reads,
                            "Adeno 3")
fragment_active_vs_inactive_sub150("D-1578-cfChIP.bam",
                            grs,
                            "Coverage of AVENIO genes.txt",
                            "PosDeduped-D-1578-input.bam",
                            healthy_reads,
                            "Adeno 4")
fragment_active_vs_inactive_sub150("E-439-cfChIP.bam",
                            grs,
                            "Coverage of AVENIO genes.txt",
                            "PosDeduped-E-439-input.bam",
                            healthy_reads,
                            "Plano 1")
fragment_active_vs_inactive_sub150("F-1449-cfChIP.bam",
                            grs,
                            "Coverage of AVENIO genes.txt",
                            "PosDeduped-F-1449-input.bam",
                            healthy_reads,
                            "Plano 2")
fragment_active_vs_inactive_sub150("I-645-cfChIP.bam",
                            grs,
                            "Coverage of AVENIO genes.txt",
                            "PosDeduped-I-645-input.bam",
                            healthy_reads,
                            "Plano 3")
fragment_active_vs_inactive_sub150("J-1663-cfChIP.bam",
                            grs,
                            "Coverage of AVENIO genes.txt",
                            "PosDeduped-J-1663-input.bam",
                            healthy_reads,
                            "Plano 4")
fragment_active_vs_inactive_sub150("G-514-cfChIP.bam",
                            grs,
                            "Coverage of AVENIO genes.txt",
                            "PosDeduped-G-514-input.bam",
                            healthy_reads,
                            "SCLC 1")
fragment_active_vs_inactive_sub150("H-1169-cfChIP.bam",
                            grs,
                            "Coverage of AVENIO genes.txt",
                            "PosDeduped-H-1169-input.bam",
                            healthy_reads,
                            "SCLC 2")
fragment_active_vs_inactive_sub150("K-440-cfChIP.bam",
                            grs,
                            "Coverage of AVENIO genes.txt",
                            "PosDeduped-K-440-input.bam",
                            healthy_reads,
                            "SCLC 3")
fragment_active_vs_inactive_sub150("L-1100-cfChIP.bam",
                            grs,
                            "Coverage of AVENIO genes.txt",
                            "PosDeduped-K-440-input.bam",
                            healthy_reads,
                            "SCLC 4")

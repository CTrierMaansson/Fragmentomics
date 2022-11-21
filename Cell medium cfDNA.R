library(ggplot2)
th <- theme(
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

setwd("C:/Users/Christoffer/OneDrive/1PhD/RNA-seq/BGI/RNA-seq 10102022")
`%ni%` <- Negate(`%in%`)
TPM <- read.table("Gene abundance 03112022.txt", header = T)
colnames(TPM)[8:10] <- c("HCC827-MET_R1", "HCC827-MET_R2", "HCC827-MET_R3")
setwd("C:/Users/Christoffer/OneDrive/1PhD/RNA-seq/BGI/Data/Gene abundance files")

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

setwd("C:/Users/Christoffer/OneDrive/1PhD/R files/Important excel and text documents")
grs <- gr("AVENIO_genes.txt")

x #Data.frame with TPM data
y #vector of length 2 with root column names to be analyzed
dif_gene_express <- function(x,y){
    library(dplyr)
    library(parameters)
    df <- x
    group1_df <- df[grepl(y[1],colnames(df))]
    group2_df <- df[grepl(y[2],colnames(df))]
    rownames(group1_df) <- df$SYMBOL
    rownames(group2_df) <- df$SYMBOL
    group1_t <- as.data.frame(t(group1_df))
    group2_t <- as.data.frame(t(group2_df))
    p <- c()
    group1 <- c()
    group2 <- c()
    se_group1 <- c()
    se_group2 <- c()
    for (i in 1:(length(group1_t))){
        g1 <- group1_t[,i]
        g2 <- group2_t[,i]
        se_group1[i] <- parameters::standard_error(g1)
        se_group2[i] <- parameters::standard_error(g2)
        p[i] <- t.test(g2,g1, paired = F, alternative = "two.sided")$p.value
        group1[i] <- mean(g1)
        group2[i] <- mean(g2)
    }
    q <- p.adjust(p, method = "fdr")
    kk <- data.frame(gene = colnames(group1_t),
                     group1_name = group1,
                     group1_se = se_group1,
                     group2_name = group2,
                     group_2_se = se_group2,
                     p.value = p,
                     q.value = q,
                     Log2FC = log2(group1+1)-log2(group2+1))
    colnames(kk)[2] <- paste("mean_TPM_",strsplit(y[1],"_")[[1]], sep = "")
    colnames(kk)[3] <- paste("se_TPM_",strsplit(y[1],"_")[[1]], sep = "")
    colnames(kk)[4] <- paste("mean_TPM_",strsplit(y[2],"_")[[1]], sep = "")
    colnames(kk)[5] <- paste("se_TPM_",strsplit(y[2],"_")[[1]], sep = "")
    kk <- kk[order(-abs(kk$Log2FC)),]
    return(kk)
}
TPM_reduced <- TPM[TPM$SYMBOL %in% grs$SYMBOL,]
HCC827_vs_A549 <- dif_gene_express(TPM_reduced, c("HCC827_", "A549_"))
volcano_dif <- function(x,y){
    library(ggplot2)
    library(ggrepel)
    p <- -log10(x$q.value)
    x$Log10p <- p
    x <- x %>% dplyr::filter(!is.na(q.value))
    which_m <- apply(x, MARGIN =1, FUN = function(z){
        g1 <- as.numeric(z[2])
        g2 <- as.numeric(z[4])
        return(log2(max(g1,g2)+1))
    })
    x$max_log2_TPM <- which_m
    group1_high <- length(x[x$Log2FC>0,]$Log2FC)
    group2_high <- length(x[x$Log2FC<0,]$Log2FC)
    gg <- ggplot(data = x, aes(x = Log2FC, y = Log10p,
                               size = max_log2_TPM, fill = Log2FC))+
        geom_point(shape = 21,
                   stroke = 0.5,)+
        scale_size_continuous(name = expression(bold(log[2]("TPM+1"))),
                              limits = c(0,15))+
        scale_fill_gradient2(low = "#ffa10c", high = "#6a00fc",
                             mid = "grey",
                             name = expression(bold(log[2]("FC"))))+
        guides(colour = guide_colourbar(order = 1),
               size = guide_legend(order = 2))+
        xlab(expression(bold(log[2]("FC"))))+
        ylab(expression(bold(paste("-",log[10]("q-value"), sep = ""))))+
        labs(title = y)+
        geom_vline(xintercept = c(0,-1,1),
                   linetype = c("solid","dashed","dashed"),
                   colour = c("black", "black", "black"),
                   size = c(1,1,1))+
        geom_label(x = 4, y = 0, label = paste(strsplit(colnames(x)[2],"_")[[1]][3]),
                   color = "#6a00fc", label.size = 0, size = 5,
                   fill="white")+
        geom_label(x = -4, y = 0, label = paste(strsplit(colnames(x)[4],"_")[[1]][3]),
                   color = "#ffa10c", label.size = 0, size = 5,
                   fill="white")+
        theme_bw(base_size = 17)+
        scale_x_continuous(breaks = seq(-10, 10, by = 2), limits = c(-10,10))+
        geom_label_repel(
            aes(label=ifelse(Log2FC > 1, as.character(gene),"")),
            segment.color="#6a00fc",
            color="#6a00fc",
            nudge_x = 1,
            label.size = NA,
            size = 2.5,
            fill = alpha(c("white"),0),
            parse = F,
            max.overlaps = 100)+
        geom_label_repel(
            aes(label=ifelse(Log2FC < -1, as.character(gene),"")),
            segment.color="#ffa10c",
            color="#ffa10c",
            nudge_x = -1,
            label.size = NA,
            size = 2.5,
            fill = alpha(c("white"),0),
            parse = F,
            max.overlaps = 100)
    return(gg)
}
volcano_dif(HCC827_vs_A549, "HCC827 compared to A549")

setwd("D:/Cell medium DNA")
bam <- function(q,z){
    library(chromstaR)
    file <- chromstaR::readBamFileAsGRanges(
        q, bamindex = z,
        chromosomes = NULL,
        pairedEndReads = F,
        remove.duplicate.reads = T, 
        blacklist = NULL,
        what = c("isize", "mapq"))
    return(file)
}
A549_cfDNA_R1_bam <- bam("Deduped-A549_1_cell_free_medium.bam",
                         "Deduped-A549_1_cell_free_medium.bam.bai")
A549_cfDNA_R2_bam <- bam("Deduped-A549_2_cell_free_medium.bam",
                         "Deduped-A549_2_cell_free_medium.bam.bai")
A549_cfDNA_R3_bam <- bam("Deduped-A549_3_cell_free_medium.bam",
                         "Deduped-A549_3_cell_free_medium.bam.bai")
HCC827_cfDNA_R1_bam <- bam("Deduped-HCC827_1_cell_free_medium.bam",
                           "Deduped-HCC827_1_cell_free_medium.bam.bai")
HCC827_cfDNA_R2_bam <- bam("Deduped-HCC827_2_cell_free_medium.bam",
                           "Deduped-HCC827_2_cell_free_medium.bam.bai")
HCC827_cfDNA_R3_bam <- bam("Deduped-HCC827_3_cell_free_medium.bam",
                           "Deduped-HCC827_3_cell_free_medium.bam.bai")
mcols(A549_cfDNA_R1_bam)$size <- abs(mcols(A549_cfDNA_R1_bam)$isize)
mcols(A549_cfDNA_R2_bam)$size <- abs(mcols(A549_cfDNA_R2_bam)$isize)
mcols(A549_cfDNA_R3_bam)$size <- abs(mcols(A549_cfDNA_R3_bam)$isize)
mcols(HCC827_cfDNA_R1_bam)$size <- abs(mcols(HCC827_cfDNA_R1_bam)$isize)
mcols(HCC827_cfDNA_R2_bam)$size <- abs(mcols(HCC827_cfDNA_R2_bam)$isize)
mcols(HCC827_cfDNA_R3_bam)$size <- abs(mcols(HCC827_cfDNA_R3_bam)$isize)

size_df <- data.frame(cell = c(rep("A549", nrow(mcols(A549_cfDNA_R1_bam))),
                               rep("A549", nrow(mcols(A549_cfDNA_R2_bam))),
                               rep("A549", nrow(mcols(A549_cfDNA_R3_bam))),
                               rep("HCC827", nrow(mcols(HCC827_cfDNA_R1_bam))),
                               rep("HCC827", nrow(mcols(HCC827_cfDNA_R2_bam))),
                               rep("HCC827", nrow(mcols(HCC827_cfDNA_R3_bam)))),
                      size = c(mcols(A549_cfDNA_R1_bam)$size,
                               mcols(A549_cfDNA_R2_bam)$size,
                               mcols(A549_cfDNA_R3_bam)$size,
                               mcols(HCC827_cfDNA_R1_bam)$size,
                               mcols(HCC827_cfDNA_R2_bam)$size,
                               mcols(HCC827_cfDNA_R3_bam)$size))
library(ggplot2)

size_df %>% 
    ggplot(aes(fill = cell))+
    geom_histogram(aes(size),
                   binwidth = 5)+
    theme_bw()

x #granges object returned by bam and size is in positve values
y #cutoff for fragment length

sub_y_reads <- function(x,y){
    x <- x[!is.na(mcols(x)$size),]
    res <- x[mcols(x)$size <= y,]
    return(res)
}
A549_cfDNA_R1_sub150 <- sub_y_reads(A549_cfDNA_R1_bam, 150) 
A549_cfDNA_R2_sub150 <- sub_y_reads(A549_cfDNA_R2_bam, 150) 
A549_cfDNA_R3_sub150 <- sub_y_reads(A549_cfDNA_R3_bam, 150) 
HCC827_cfDNA_R1_sub150 <- sub_y_reads(HCC827_cfDNA_R1_bam, 150) 
HCC827_cfDNA_R2_sub150 <- sub_y_reads(HCC827_cfDNA_R2_bam, 150) 
HCC827_cfDNA_R3_sub150 <- sub_y_reads(HCC827_cfDNA_R3_bam, 150) 

size_df_sub150 <- data.frame(cell = c(rep("A549", nrow(mcols(A549_cfDNA_R1_sub150))),
                               rep("A549", nrow(mcols(A549_cfDNA_R2_sub150))),
                               rep("A549", nrow(mcols(A549_cfDNA_R3_sub150))),
                               rep("HCC827", nrow(mcols(HCC827_cfDNA_R1_sub150))),
                               rep("HCC827", nrow(mcols(HCC827_cfDNA_R2_sub150))),
                               rep("HCC827", nrow(mcols(HCC827_cfDNA_R3_sub150)))),
                      size = c(mcols(A549_cfDNA_R1_sub150)$size,
                               mcols(A549_cfDNA_R2_sub150)$size,
                               mcols(A549_cfDNA_R3_sub150)$size,
                               mcols(HCC827_cfDNA_R1_sub150)$size,
                               mcols(HCC827_cfDNA_R2_sub150)$size,
                               mcols(HCC827_cfDNA_R3_sub150)$size))
size_df_sub150 %>% 
    ggplot(aes(fill = cell))+
    geom_histogram(aes(size),
                   binwidth = 5)+
    theme_bw()


x #Granges object returned by grs()
y #Granges object returned by sub_y_reads() or bam()

AVENIO_overlaps <- function(x, y) {
    library(GenomicAlignments)
    res <- summarizeOverlaps(features = x,
                             reads = y,
                             mode = "Union",
                             ignor.strand = TRUE)
    df <- data.frame(genes = grs$SYMBOL,
                     counts = assays(res)$counts[,1])
    return(df)
}
A549_R1_sub150_counts <- AVENIO_overlaps(grs, A549_cfDNA_R1_sub150)
A549_R2_sub150_counts <- AVENIO_overlaps(grs, A549_cfDNA_R2_sub150)
A549_R3_sub150_counts <- AVENIO_overlaps(grs, A549_cfDNA_R3_sub150)
HCC827_R1_sub150_counts <- AVENIO_overlaps(grs, HCC827_cfDNA_R1_sub150)
HCC827_R2_sub150_counts <- AVENIO_overlaps(grs, HCC827_cfDNA_R2_sub150)
HCC827_R3_sub150_counts <- AVENIO_overlaps(grs, HCC827_cfDNA_R3_sub150)
A549_R1_counts <- AVENIO_overlaps(grs, A549_cfDNA_R1_bam)
A549_R2_counts <- AVENIO_overlaps(grs, A549_cfDNA_R2_bam)
A549_R3_counts <- AVENIO_overlaps(grs, A549_cfDNA_R3_bam)
HCC827_R1_counts <- AVENIO_overlaps(grs, HCC827_cfDNA_R1_bam)
HCC827_R2_counts <- AVENIO_overlaps(grs, HCC827_cfDNA_R2_bam)
HCC827_R3_counts <- AVENIO_overlaps(grs, HCC827_cfDNA_R3_bam)

setwd("C:/Users/Christoffer/OneDrive/1PhD/R files/Important excel and text documents")
x #Gene count table returned by gene_count()
y #.txt file of the AVENIO genes with the number of bases sequenced in each gene. based on data.frame returned by coverages()
enrichment <- function(x,y){
    library(GenomicAlignments)
    library(GenomicRanges)
    library(BiocParallel)
    library(dplyr)
    y <- read.table(y, header=T)
    y$coverage <- as.numeric(y$coverage)
    sub_reads <- x
    sub_reads <- sub_reads %>% filter(counts > 10)
    y <- y[match(sub_reads$genes, y$SYMBOL),]
    len <- sum(sub_reads$counts)
    RPKM <- c()
    for (i in 1:length(sub_reads$genes)){
        RPKM[i] <- (sub_reads$counts[i]*1000*1000000)/(len*y$coverage[i])
    }
    sub_reads$e <- RPKM
    res <- data.frame(genes = sub_reads$genes, enrichment=sub_reads$e)
    res <- res[order(res$enrichment),]
    return(res)
}
A_R1_e <- enrichment(A549_R1_sub150_counts, 
           "Coverage of AVENIO genes.txt")
TPM_reduced_A <- TPM_reduced %>% filter(SYMBOL %in% A_R1_e$genes)
try_df <- data.frame(TPM = TPM_reduced_A[order(A_R1_e$genes, TPM_reduced_A$SYMBOL),]$A549_R1, 
                     enrichment = A_R1_e$enrichment)
try_df %>% ggplot(aes(x = enrichment,
                      y = TPM))+
    geom_point()


A549_R1_counts
x #gene counts for all genes returned by AVENIO_overlaps()
y #gene counts for all genes of sub150 reads returned by AVENIO_overlaps()
p_sub_y <- function(x,y) {
    x <- x %>% filter(counts > 10)
    y <- y %>% filter(counts > 10)
    y <- y %>% filter(genes %in% x$genes)
    x <- x %>% filter(genes %in% y$genes)
    p <- y$counts/x$counts *100
    df <- data.frame(genes = y$genes,
                     sub_fraction = p)
    return(df)
}
fraction_A549_R1 <- p_sub_y(A549_R1_counts, A549_R1_sub150_counts)
fraction_A549_R2 <- p_sub_y(A549_R2_counts, A549_R2_sub150_counts)
fraction_A549_R3 <- p_sub_y(A549_R3_counts, A549_R3_sub150_counts)
fraction_HCC827_R1 <- p_sub_y(HCC827_R1_counts, HCC827_R1_sub150_counts)
fraction_HCC827_R2 <- p_sub_y(HCC827_R2_counts, HCC827_R2_sub150_counts)
fraction_HCC827_R3 <- p_sub_y(HCC827_R3_counts, HCC827_R3_sub150_counts)

x # List of sub fractions returned by p_sub_y()
mean_fraction <- function(x){
    df <- data.frame(genes = x[[1]]$genes)
    for (i in 1:length(x)){
        df[ , ncol(df) + 1] <- x[[i]]$sub_fraction
        colnames(df)[ncol(df)] <- paste0(i)
    }
    df_values <- df %>% dplyr::select(-genes)
    means <- rowMeans(df_values)
    res <- data.frame(genes = df$genes,
                      sub_fraction = means)
    return(res)
}

mean_A549_fraction <- mean_fraction(list(fraction_A549_R1, 
                   fraction_A549_R2, 
                   fraction_A549_R3))
mean_HCC827_fraction <- mean_fraction(list(fraction_HCC827_R1, 
                                           fraction_HCC827_R2, 
                                           fraction_HCC827_R3))
x #Name on combined TPM file
y #vector of names on samples
RNA_TPM <- function(x,y){
    library(dplyr)
    options(scipen=999)
    TPM <- read.table(x, header = T)
    df <- data.frame(SYMBOL = TPM$SYMBOL)
    for (i in 2:length(TPM)){
        res <- log2(TPM[,i]+1)
        df[ , ncol(df) + 1] <- res
        colnames(df)[ncol(df)] <- paste0(i)
    }
    setwd("C:/Users/Christoffer/OneDrive/1PhD/R files/Important excel and text documents")
    AVENIO_gr <- gr("AVENIO_genes.txt")
    AVENIO_genelist <- AVENIO_gr$SYMBOL
    TPM_AVENIO <- df %>% filter(SYMBOL %in% AVENIO_genelist)
    colnames(TPM_AVENIO) <- c("SYMBOL",y)
    setwd("C:/Users/Christoffer/OneDrive/1PhD/RNA-seq/BGI/RNA-seq 10102022")
    return(TPM_AVENIO)
}
setwd("C:/Users/Christoffer/OneDrive/1PhD/RNA-seq/BGI/RNA-seq 10102022")
TPM_AVENIO <- RNA_TPM("Gene abundance 03112022.txt",c("log2_A549_R1",
                                                      "log2_A549_R2",
                                                      "log2_A549_R3",
                                                      "log2_HCC827_R1",
                                                      "log2_HCC827_R2",
                                                      "log2_HCC827_R3"))
x #data.frame returned by RNA_TPM()
y #colname of cell line to gen the mean log2(TPM+1)
RNA_mean <- function(x,y){
    df <- x[grepl(y,colnames(x))]
    return(data.frame(SYMBOL = x$SYMBOL,
                      log2_TPM = rowMeans(df)))
}
mean_A549_TPM <- RNA_mean(TPM_AVENIO, "log2_A549_")
mean_HCC827_TPM <- RNA_mean(TPM_AVENIO, "log2_HCC827_")
mean_A549_TPM

x #mean log2TPM values for cell line returned by RNA_mean
y #mean sub fractions for cell line returned by mean_fraction()
z #cut-off for defintion of active gene
box_fraction <- function(x,y,z, c){
    library(ggpubr)
    x <- x %>% filter(SYMBOL %in% y$genes)
    status <- c()
    for (i in 1:nrow(x)){
        if(x$log2_TPM[i] > z){
            status[i] <- "Active"
        }
        else{
            status[i] <- "Inactive"
        }
    }
    y$status <- status
    gg1 <- y %>% group_by(status) %>% 
        ggplot(aes(x = status, y = sub_fraction, fill = status))+
        geom_boxplot(outlier.alpha = 0)+
        geom_jitter(alpha = 0.5, width = 0.1)+
        theme_bw()+
        stat_compare_means(method = "t.test",
                           label.x = 1.5, 
                           label.y = c(max(y$sub_fraction)+2,
                                       max(y$sub_fraction)+1,
                                       max(y$sub_fraction)+1.5),
                           comparisons = list(c("Active", "Inactive")))+
        labs(title = paste(c),
             caption = paste0("Active genes: n = ", sum(y$status == "Active"),
                              ", Inactive genes: n = ", sum(y$status == "Inactive")),
             x = "",
             y = "Sub 150 fraction (%)")+
        scale_fill_manual("Gene activity", values = c("Active" = "green4", 
                                                      "Inactive" = "firebrick"))+
        th
    return(gg1)

}
box_fraction(mean_A549_TPM, mean_A549_fraction, 1, "A549") 
box_fraction(mean_HCC827_TPM, mean_HCC827_fraction, 1, "HCC827") 

####Percent of input####

A549_pippin_R1_bam <- bam("Deduped-Pippin_A1P.bam",
                          "Deduped-Pippin_A1P.bam.bai")
A549_pippin_R2_bam <- bam("Deduped-Pippin_A2P.bam",
                          "Deduped-Pippin_A2P.bam.bai")
A549_pippin_R3_bam <- bam("Deduped-Pippin_A3P.bam",
                          "Deduped-Pippin_A3P.bam.bai")
HCC827_pippin_R1_bam <- bam("Deduped-Pippin_H1P.bam",
                            "Deduped-Pippin_H1P.bam.bai")
HCC827_pippin_R2_bam <- bam("Deduped-Pippin_H2P.bam",
                            "Deduped-Pippin_H2P.bam.bai")
HCC827_pippin_R3_bam <- bam("Deduped-Pippin_H3P.bam",
                            "Deduped-Pippin_H3P.bam.bai")
mcols(A549_pippin_R1_bam)$size <- abs(mcols(A549_pippin_R1_bam)$isize)
mcols(A549_pippin_R2_bam)$size <- abs(mcols(A549_pippin_R2_bam)$isize)
mcols(A549_pippin_R3_bam)$size <- abs(mcols(A549_pippin_R3_bam)$isize)
mcols(HCC827_pippin_R1_bam)$size <- abs(mcols(HCC827_pippin_R1_bam)$isize)
mcols(HCC827_pippin_R2_bam)$size <- abs(mcols(HCC827_pippin_R2_bam)$isize)
mcols(HCC827_pippin_R3_bam)$size <- abs(mcols(HCC827_pippin_R3_bam)$isize)

A549_pippin_R1_counts <- AVENIO_overlaps(grs, A549_pippin_R1_bam)
A549_pippin_R2_counts <- AVENIO_overlaps(grs, A549_pippin_R2_bam)
A549_pippin_R3_counts <- AVENIO_overlaps(grs, A549_pippin_R3_bam)
HCC827_pippin_R1_counts <- AVENIO_overlaps(grs, HCC827_pippin_R1_bam)
HCC827_pippin_R2_counts <- AVENIO_overlaps(grs, HCC827_pippin_R2_bam)
HCC827_pippin_R3_counts <- AVENIO_overlaps(grs, HCC827_pippin_R3_bam)

x # counts of input sample returned by AVENIO_overlaps()
y # counts of pippin sample returned by AVENIO_overlaps()
percent_of_input <- function(x,y){
    x <- x %>% filter(counts > 20)
    y <- y %>% filter(counts > 20)
    y <- y %>% filter(genes %in% x$genes)
    x <- x %>% filter(genes %in% y$genes)
    percent <- y$counts/x$counts *100
    df <- data.frame(genes = x$genes, 
                     percent = percent)
    return(df)
}
A549_percent_R1 <- percent_of_input(A549_R1_counts, A549_pippin_R1_counts)
A549_percent_R2 <- percent_of_input(A549_R2_counts, A549_pippin_R2_counts)
A549_percent_R3 <- percent_of_input(A549_R3_counts, A549_pippin_R3_counts)
HCC827_percent_R1 <- percent_of_input(HCC827_R1_counts, HCC827_pippin_R1_counts)
HCC827_percent_R2 <- percent_of_input(HCC827_R2_counts, HCC827_pippin_R2_counts)
HCC827_percent_R3 <- percent_of_input(HCC827_R3_counts, HCC827_pippin_R3_counts)

x # List of sub fractions returned by p_sub_y()
mean_percent <- function(x){
    df <- data.frame(genes = x[[1]]$genes)
    for (i in 1:length(x)){
        df[ , ncol(df) + 1] <- x[[i]]$percent
        colnames(df)[ncol(df)] <- paste0(i)
    }
    df_values <- df %>% dplyr::select(-genes)
    means <- rowMeans(df_values)
    res <- data.frame(genes = df$genes,
                      percent = means)
    return(res)
}

mean_A549_percent <- mean_percent(list(A549_percent_R1,
                          A549_percent_R2,
                          A549_percent_R3))
mean_HCC827_percent <- mean_percent(list(HCC827_percent_R1,
                          HCC827_percent_R2,
                          HCC827_percent_R3))
x #mean log2TPM values for cell line returned by RNA_mean
y #mean sub fractions for cell line returned by mean_fraction()
z #cut-off for defintion of active gene
c #name of cell line
box_percent <- function(x,y,z, c){
    library(ggpubr)
    x <- x %>% filter(SYMBOL %in% y$genes)
    status <- c()
    for (i in 1:nrow(x)){
        if(x$log2_TPM[i] > z){
            status[i] <- "Active"
        }
        else{
            status[i] <- "Inactive"
        }
    }
    y$status <- status
    gg1 <- y %>% group_by(status) %>% 
        ggplot(aes(x = status, y = percent, fill = status))+
        geom_boxplot(outlier.alpha = 0)+
        geom_jitter(alpha = 0.5, width = 0.1)+
        theme_bw()+
        stat_compare_means(method = "t.test",
                           label.x = 1.5, 
                           label.y = c(max(y$percent)+2,
                                       max(y$percent)+1,
                                       max(y$percent)+1.5),
                           comparisons = list(c("Active", "Inactive")))+
        labs(title = paste(c),
             caption = paste0("Active genes: n = ", sum(y$status == "Active"),
                              ", Inactive genes: n = ", sum(y$status == "Inactive")),
             x = "",
             y = "% of input")+
        scale_fill_manual("Gene activity", values = c("Active" = "green4", 
                                                      "Inactive" = "firebrick"))+
        th
    return(gg1)
    
}

box_percent(mean_A549_TPM, mean_A549_percent, 1, "A549")
box_percent(mean_HCC827_TPM, mean_HCC827_percent, 1, "HCC827")

####Top and low expressed genes####
mean_A549_TPM
A_TPM <- mean_A549_TPM %>% arrange(-log2_TPM)
A_TPM_zero <- mean_A549_TPM  %>% filter(log2_TPM == 0)

A_top <- A_TPM$SYMBOL[1:15]
A_zero <- A_TPM_zero$SYMBOL

H_TPM <- mean_HCC827_TPM %>% arrange(-log2_TPM)
H_TPM_zero <- mean_HCC827_TPM  %>% filter(log2_TPM == 0)

H_top <- H_TPM$SYMBOL[1:15]
H_zero <- H_TPM_zero$SYMBOL

mean_A_R1_top <- A549_percent_R1 %>% filter(genes %in% A_top) %>% 
    summarise(mean(percent)) 
mean_A_R2_top <- A549_percent_R2 %>% filter(genes %in% A_top) %>% 
    summarise(mean(percent)) 
mean_A_R3_top <- A549_percent_R3 %>% filter(genes %in% A_top) %>% 
    summarise(mean(percent)) 
mean_H_R1_top <- HCC827_percent_R1 %>% filter(genes %in% H_top) %>% 
    summarise(mean(percent)) 
mean_H_R2_top <- HCC827_percent_R2 %>% filter(genes %in% H_top) %>% 
    summarise(mean(percent)) 
mean_H_R3_top <- HCC827_percent_R3 %>% filter(genes %in% H_top) %>% 
    summarise(mean(percent)) 

mean_A_R1_zero <- A549_percent_R1 %>% filter(genes %in% A_zero) %>% 
    summarise(mean(percent)) 
mean_A_R2_zero <- A549_percent_R2 %>% filter(genes %in% A_zero) %>% 
    summarise(mean(percent)) 
mean_A_R3_zero <- A549_percent_R3 %>% filter(genes %in% A_zero) %>% 
    summarise(mean(percent)) 
mean_H_R1_zero <- HCC827_percent_R1 %>% filter(genes %in% H_zero) %>% 
    summarise(mean(percent)) 
mean_H_R2_zero <- HCC827_percent_R2 %>% filter(genes %in% H_zero) %>% 
    summarise(mean(percent)) 
mean_H_R3_zero <- HCC827_percent_R3 %>% filter(genes %in% H_zero) %>% 
    summarise(mean(percent)) 
fraction_A549_R1
mean_A_R1_top
mean_A_R1_zero
mean_A_R2_top
mean_A_R2_zero
mean_A_R3_top
mean_A_R3_zero

mean_H_R1_top
mean_H_R1_zero
mean_H_R2_top
mean_H_R2_zero
mean_H_R3_top
mean_H_R3_zero


df_top_zero <- data.frame(cell = c(rep("A549",3),
                                   rep("HCC827",3),
                                   rep("A549",3),
                                   rep("HCC827",3)),
                          Gene_activity = c(rep("Zero", 6),
                                            rep("High", 6)),
                          percent = c(as.numeric(mean_A_R1_zero),
                                      as.numeric(mean_A_R2_zero),
                                      as.numeric(mean_A_R3_zero),
                                      as.numeric(mean_H_R1_zero),
                                      as.numeric(mean_H_R2_zero),
                                      as.numeric(mean_H_R3_zero),
                                      as.numeric(mean_A_R1_top),
                                      as.numeric(mean_A_R2_top),
                                      as.numeric(mean_A_R3_top),
                                      as.numeric(mean_H_R1_top),
                                      as.numeric(mean_H_R2_top),
                                      as.numeric(mean_H_R3_top)),
                          pairing = rep(c(LETTERS[1],
                                      LETTERS[2],
                                      LETTERS[3],
                                      LETTERS[4],
                                      LETTERS[5],
                                      LETTERS[6]),2))
df_top_zero %>% group_by(Gene_activity) %>% 
    ggplot(aes(x = Gene_activity, y = sub_fraction,
               color = Gene_activity, shape = cell))+
    geom_line(aes(group = pairing),
              color = "black", size = 1)+
    geom_point(size = 4)+
    scale_color_manual("Gene activity", values = c("High" = "green4", 
                                                  "Zero" = "firebrick"))+
    theme_bw()+
    labs(title ="Pippin % of input in active\nand inactive genes",
         x = "",
         y = "% of input")+
    th

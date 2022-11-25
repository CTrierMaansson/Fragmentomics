
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
setwd("D:/Lung cancer cfChIP/PosDeduped")
bamfile <- function(x,z){
    library(chromstaR)
    file <- readBamFileAsGRanges(
        x, bamindex = z,
        chromosomes = NULL,
        pairedEndReads = F,
        remove.duplicate.reads = T, 
        blacklist = NULL,
        what = c("mapq", "isize", "seq"))
}    
NAC.1_cfChIP <- bamfile("PosDeduped-A-1279-cfChIP.bam",
                        "PosDeduped-A-1279-cfChIP.bam.bai")
NAC.2_cfChIP <- bamfile("PosDeduped-B-1288-cfChIP.bam",
                        "PosDeduped-B-1288-cfChIP.bam.bai")
NAC.3_cfChIP <- bamfile("PosDeduped-C-1475-cfChIP.bam",
                        "PosDeduped-C-1475-cfChIP.bam.bai")
NAC.4_cfChIP <- bamfile("PosDeduped-D-1578-cfChIP.bam",
                        "PosDeduped-D-1578-cfChIP.bam.bai")
NSC.1_cfChIP <- bamfile("PosDeduped-E-439-cfChIP.bam",
                        "PosDeduped-E-439-cfChIP.bam.bai")
NSC.2_cfChIP <- bamfile("PosDeduped-F-1449-cfChIP.bam",
                        "PosDeduped-F-1449-cfChIP.bam.bai")
NSC.3_cfChIP <- bamfile("PosDeduped-I-645-cfChIP.bam",
                        "PosDeduped-I-645-cfChIP.bam.bai")
NSC.4_cfChIP <- bamfile("PosDeduped-J-1663-cfChIP.bam",
                        "PosDeduped-J-1663-cfChIP.bam.bai")
SSC.1_cfChIP <- bamfile("PosDeduped-G-514-cfChIP.bam",
                        "PosDeduped-G-514-cfChIP.bam.bai")
SSC.2_cfChIP <- bamfile("PosDeduped-H-1169-cfChIP.bam",
                        "PosDeduped-H-1169-cfChIP.bam.bai")
SSC.3_cfChIP <- bamfile("PosDeduped-K-440-cfChIP.bam",
                        "PosDeduped-K-440-cfChIP.bam.bai")
SSC.4_cfChIP <- bamfile("PosDeduped-L-1100-cfChIP.bam",
                        "PosDeduped-L-1100-cfChIP.bam.bai")
setwd("D:/Healthy cfChIP/PosDeduped")
HC.1_cfChIP <- bamfile("PosDeduped-Rask_kontrol_1_cfChIP.bam",
                       "PosDeduped-Rask_kontrol_1_cfChIP.bam.bai")
HC.2_cfChIP <- bamfile("PosDeduped-Rask_kontrol_2_cfChIP.bam",
                       "PosDeduped-Rask_kontrol_2_cfChIP.bam.bai")
HC.3_cfChIP <- bamfile("PosDeduped-Rask_kontrol_3_cfChIP.bam",
                       "PosDeduped-Rask_kontrol_3_cfChIP.bam.bai")
HC.4_cfChIP <- bamfile("PosDeduped-Rask_kontrol_4_cfChIP.bam",
                       "PosDeduped-Rask_kontrol_4_cfChIP.bam.bai")

setwd("D:/Lung cancer input/PosDeduped")
NAC.1_input <- bamfile("PosDeduped-A-1279-input.bam",
                       "PosDeduped-A-1279-input.bam.bai")
NAC.2_input <- bamfile("PosDeduped-B-1288-input.bam",
                       "PosDeduped-B-1288-input.bam.bai")
NAC.3_input <- bamfile("PosDeduped-C-1475-input.bam",
                       "PosDeduped-C-1475-input.bam.bai")
NAC.4_input <- bamfile("PosDeduped-D-1578-input.bam",
                       "PosDeduped-D-1578-input.bam.bai")
NSC.1_input <- bamfile("PosDeduped-E-439-input.bam",
                       "PosDeduped-E-439-input.bam.bai")
NSC.2_input <- bamfile("PosDeduped-F-1449-input.bam",
                       "PosDeduped-F-1449-input.bam.bai")
NSC.3_input <- bamfile("PosDeduped-I-645-input.bam",
                       "PosDeduped-I-645-input.bam.bai")
NSC.4_input <- bamfile("PosDeduped-J-1663-input.bam",
                       "PosDeduped-J-1663-input.bam.bai")
SSC.1_input <- bamfile("PosDeduped-G-514-input.bam",
                       "PosDeduped-G-514-input.bam.bai")
SSC.2_input <- bamfile("PosDeduped-H-1169-input.bam",
                       "PosDeduped-H-1169-input.bam.bai")
SSC.3_input <- bamfile("PosDeduped-K-440-input.bam",
                       "PosDeduped-K-440-input.bam.bai")
SSC.4_input <- bamfile("PosDeduped-L-1100-input.bam",
                       "PosDeduped-L-1100-input.bam.bai")

setwd("D:/Healthy input/PosDeduped")
HC.1_input <- bamfile("PosDeduped-Rask_kontrol_1_input.bam",
                      "PosDeduped-Rask_kontrol_1_input.bam.bai")
HC.2_input <- bamfile("PosDeduped-Rask_kontrol_2_input.bam",
                      "PosDeduped-Rask_kontrol_2_input.bam.bai")
HC.3_input <- bamfile("PosDeduped-Rask_kontrol_3_input.bam",
                      "PosDeduped-Rask_kontrol_3_input.bam.bai")
HC.4_input <- bamfile("PosDeduped-Rask_kontrol_4_input.bam",
                      "PosDeduped-Rask_kontrol_4_input.bam.bai")

x # BAM file returned by bamfile()
m #number of bases to be evaluated in each fragment end
mer_count <- function(x,m){
    df <- data.frame(seq = mcols(x)$seq,
                     strand = strand(x))
    mer_fun <- function(def){
        if (def[2] == "+"){
            seqs <- substring(def[1], 1, m)
        }
        if (def[2] == "-"){
            seqs <- substring(def[1],(nchar(def[1])-m+1) , nchar(def[1]))
        }
        if (grepl("N", seqs)){
            seqs <- NA
        }
        return(seqs)
    }
    res <- apply(df, MARGIN = 1, FUN = mer_fun)
    res <- res[!is.na(res)]
    return(table(res))
    
}
warnings()
mer_count_NAC.1_cfChIP <- mer_count(NAC.1_cfChIP,3)
mer_count_NAC.2_cfChIP <- mer_count(NAC.2_cfChIP,3)
mer_count_NAC.3_cfChIP <- mer_count(NAC.3_cfChIP,3)
mer_count_NAC.4_cfChIP <- mer_count(NAC.4_cfChIP,3)
mer_count_NSC.1_cfChIP <- mer_count(NSC.1_cfChIP,3)
mer_count_NSC.2_cfChIP <- mer_count(NSC.2_cfChIP,3)
mer_count_NSC.3_cfChIP <- mer_count(NSC.3_cfChIP,3)
mer_count_NSC.4_cfChIP <- mer_count(NSC.4_cfChIP,3)
mer_count_SSC.1_cfChIP <- mer_count(SSC.1_cfChIP,3)
mer_count_SSC.2_cfChIP <- mer_count(SSC.2_cfChIP,3)
mer_count_SSC.3_cfChIP <- mer_count(SSC.3_cfChIP,3)
mer_count_SSC.4_cfChIP <- mer_count(SSC.4_cfChIP,3)
mer_count_HC.1_cfChIP <- mer_count(HC.1_cfChIP, 3)
mer_count_HC.2_cfChIP <- mer_count(HC.2_cfChIP, 3)
mer_count_HC.3_cfChIP <- mer_count(HC.3_cfChIP, 3)
mer_count_HC.4_cfChIP <- mer_count(HC.4_cfChIP, 3)

mer_count_NAC.1_input <- mer_count(NAC.1_input,3)
mer_count_NAC.2_input <- mer_count(NAC.2_input,3)
mer_count_NAC.3_input <- mer_count(NAC.3_input,3)
mer_count_NAC.4_input <- mer_count(NAC.4_input,3)
mer_count_NSC.1_input <- mer_count(NSC.1_input,3)
mer_count_NSC.2_input <- mer_count(NSC.2_input,3)
mer_count_NSC.3_input <- mer_count(NSC.3_input,3)
mer_count_NSC.4_input <- mer_count(NSC.4_input,3)
mer_count_SSC.1_input <- mer_count(SSC.1_input,3)
mer_count_SSC.2_input <- mer_count(SSC.2_input,3)
mer_count_SSC.3_input <- mer_count(SSC.3_input,3)
mer_count_SSC.4_input <- mer_count(SSC.4_input,3)
mer_count_HC.1_input <- mer_count(HC.1_input, 3)
mer_count_HC.2_input <- mer_count(HC.2_input, 3)
mer_count_HC.3_input <- mer_count(HC.3_input, 3)
mer_count_HC.4_input <- mer_count(HC.4_input, 3)

prop_df <- function(x){
    df <- as.data.frame.table(proportions(x))
    return(df)
}
gc()
prop_mer_NAC.1_cfChIP <- prop_df(mer_count_NAC.1_cfChIP)
prop_mer_NAC.2_cfChIP <- prop_df(mer_count_NAC.2_cfChIP)
prop_mer_NAC.3_cfChIP <- prop_df(mer_count_NAC.3_cfChIP)
prop_mer_NAC.4_cfChIP <- prop_df(mer_count_NAC.4_cfChIP)
prop_mer_NSC.1_cfChIP <- prop_df(mer_count_NSC.1_cfChIP)
prop_mer_NSC.2_cfChIP <- prop_df(mer_count_NSC.2_cfChIP)
prop_mer_NSC.3_cfChIP <- prop_df(mer_count_NSC.3_cfChIP)
prop_mer_NSC.4_cfChIP <- prop_df(mer_count_NSC.4_cfChIP)
prop_mer_SSC.1_cfChIP <- prop_df(mer_count_SSC.1_cfChIP)
prop_mer_SSC.2_cfChIP <- prop_df(mer_count_SSC.2_cfChIP)
prop_mer_SSC.3_cfChIP <- prop_df(mer_count_SSC.3_cfChIP)
prop_mer_SSC.4_cfChIP <- prop_df(mer_count_SSC.4_cfChIP)
prop_mer_HC.1_cfChIP <- prop_df(mer_count_HC.1_cfChIP)
prop_mer_HC.2_cfChIP <- prop_df(mer_count_HC.2_cfChIP)
prop_mer_HC.3_cfChIP <- prop_df(mer_count_HC.3_cfChIP)
prop_mer_HC.4_cfChIP <- prop_df(mer_count_HC.4_cfChIP)

prop_mer_NAC.1_input <- prop_df(mer_count_NAC.1_input)
prop_mer_NAC.2_input <- prop_df(mer_count_NAC.2_input)
prop_mer_NAC.3_input <- prop_df(mer_count_NAC.3_input)
prop_mer_NAC.4_input <- prop_df(mer_count_NAC.4_input)
prop_mer_NSC.1_input <- prop_df(mer_count_NSC.1_input)
prop_mer_NSC.2_input <- prop_df(mer_count_NSC.2_input)
prop_mer_NSC.3_input <- prop_df(mer_count_NSC.3_input)
prop_mer_NSC.4_input <- prop_df(mer_count_NSC.4_input)
prop_mer_SSC.1_input <- prop_df(mer_count_SSC.1_input)
prop_mer_SSC.2_input <- prop_df(mer_count_SSC.2_input)
prop_mer_SSC.3_input <- prop_df(mer_count_SSC.3_input)
prop_mer_SSC.4_input <- prop_df(mer_count_SSC.4_input)
prop_mer_HC.1_input <- prop_df(mer_count_HC.1_input)
prop_mer_HC.2_input <- prop_df(mer_count_HC.2_input)
prop_mer_HC.3_input <- prop_df(mer_count_HC.3_input)
prop_mer_HC.4_input <- prop_df(mer_count_HC.4_input)
library(dplyr)
prop_mer_cfChIP <- prop_mer_NAC.1_cfChIP %>% 
    left_join(prop_mer_NAC.2_cfChIP, by = "res") %>% 
    left_join(prop_mer_NAC.3_cfChIP, by = "res") %>% 
    left_join(prop_mer_NAC.4_cfChIP, by = "res") %>% 
    left_join(prop_mer_NSC.1_cfChIP, by = "res") %>% 
    left_join(prop_mer_NSC.2_cfChIP, by = "res") %>% 
    left_join(prop_mer_NSC.3_cfChIP, by = "res") %>% 
    left_join(prop_mer_NSC.4_cfChIP, by = "res") %>% 
    left_join(prop_mer_SSC.1_cfChIP, by = "res") %>% 
    left_join(prop_mer_SSC.2_cfChIP, by = "res") %>% 
    left_join(prop_mer_SSC.3_cfChIP, by = "res") %>% 
    left_join(prop_mer_SSC.4_cfChIP, by = "res") %>% 
    left_join(prop_mer_HC.1_cfChIP, by = "res") %>%
    left_join(prop_mer_HC.2_cfChIP, by = "res") %>%
    left_join(prop_mer_HC.3_cfChIP, by = "res") %>%
    left_join(prop_mer_HC.4_cfChIP, by = "res")
prop_mer_input <- prop_mer_NAC.1_input %>% 
    left_join(prop_mer_NAC.2_input, by = "res") %>% 
    left_join(prop_mer_NAC.3_input, by = "res") %>% 
    left_join(prop_mer_NAC.4_input, by = "res") %>% 
    left_join(prop_mer_NSC.1_input, by = "res") %>% 
    left_join(prop_mer_NSC.2_input, by = "res") %>% 
    left_join(prop_mer_NSC.3_input, by = "res") %>% 
    left_join(prop_mer_NSC.4_input, by = "res") %>% 
    left_join(prop_mer_SSC.1_input, by = "res") %>% 
    left_join(prop_mer_SSC.2_input, by = "res") %>% 
    left_join(prop_mer_SSC.3_input, by = "res") %>% 
    left_join(prop_mer_SSC.4_input, by = "res") %>% 
    left_join(prop_mer_HC.1_input, by = "res") %>%
    left_join(prop_mer_HC.2_input, by = "res") %>%
    left_join(prop_mer_HC.3_input, by = "res") %>%
    left_join(prop_mer_HC.4_input, by = "res")

colnames(prop_mer_cfChIP) <- c("motif", 
                               "NAC.1_cfChIP",
                               "NAC.2_cfChIP",
                               "NAC.3_cfChIP",
                               "NAC.4_cfChIP",
                               "NSC.1_cfChIP",
                               "NSC.2_cfChIP",
                               "NSC.3_cfChIP",
                               "NSC.4_cfChIP",
                               "SSC.1_cfChIP",
                               "SSC.2_cfChIP",
                               "SSC.3_cfChIP",
                               "SSC.4_cfChIP",
                               "HC.1_cfChIP",
                               "HC.2_cfChIP",
                               "HC.3_cfChIP",
                               "HC.4_cfChIP")
colnames(prop_mer_input) <- c("motif", 
                               "NAC.1_input",
                               "NAC.2_input",
                               "NAC.3_input",
                               "NAC.4_input",
                               "NSC.1_input",
                               "NSC.2_input",
                               "NSC.3_input",
                               "NSC.4_input",
                               "SSC.1_input",
                               "SSC.2_input",
                               "SSC.3_input",
                               "SSC.4_input",
                               "HC.1_input",
                               "HC.2_input",
                               "HC.3_input",
                               "HC.4_input")
library(dplyr)
setwd("C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver")
prop_mer_input <- read.table("fragment end motif proportions input.txt", header = T) %>% select(-X)
prop_mer_cfChIP <- read.table("fragment end motif proportions cfChIP.txt", header = T) %>% select(-X)

count_df <- function(x){
    df <- as.data.frame.table(x)
    return(df) 
}
count_mer_NAC.1_cfChIP <- count_df(mer_count_NAC.1_cfChIP)
count_mer_NAC.2_cfChIP <- count_df(mer_count_NAC.2_cfChIP)
count_mer_NAC.3_cfChIP <- count_df(mer_count_NAC.3_cfChIP)
count_mer_NAC.4_cfChIP <- count_df(mer_count_NAC.4_cfChIP)
count_mer_NSC.1_cfChIP <- count_df(mer_count_NSC.1_cfChIP)
count_mer_NSC.2_cfChIP <- count_df(mer_count_NSC.2_cfChIP)
count_mer_NSC.3_cfChIP <- count_df(mer_count_NSC.3_cfChIP)
count_mer_NSC.4_cfChIP <- count_df(mer_count_NSC.4_cfChIP)
count_mer_SSC.1_cfChIP <- count_df(mer_count_SSC.1_cfChIP)
count_mer_SSC.2_cfChIP <- count_df(mer_count_SSC.2_cfChIP)
count_mer_SSC.3_cfChIP <- count_df(mer_count_SSC.3_cfChIP)
count_mer_SSC.4_cfChIP <- count_df(mer_count_SSC.4_cfChIP)
count_mer_HC.1_cfChIP <- count_df(mer_count_HC.1_cfChIP)
count_mer_HC.2_cfChIP <- count_df(mer_count_HC.2_cfChIP)
count_mer_HC.3_cfChIP <- count_df(mer_count_HC.3_cfChIP)
count_mer_HC.4_cfChIP <- count_df(mer_count_HC.4_cfChIP)

count_mer_NAC.1_input <- count_df(mer_count_NAC.1_input)
count_mer_NAC.2_input <- count_df(mer_count_NAC.2_input)
count_mer_NAC.3_input <- count_df(mer_count_NAC.3_input)
count_mer_NAC.4_input <- count_df(mer_count_NAC.4_input)
count_mer_NSC.1_input <- count_df(mer_count_NSC.1_input)
count_mer_NSC.2_input <- count_df(mer_count_NSC.2_input)
count_mer_NSC.3_input <- count_df(mer_count_NSC.3_input)
count_mer_NSC.4_input <- count_df(mer_count_NSC.4_input)
count_mer_SSC.1_input <- count_df(mer_count_SSC.1_input)
count_mer_SSC.2_input <- count_df(mer_count_SSC.2_input)
count_mer_SSC.3_input <- count_df(mer_count_SSC.3_input)
count_mer_SSC.4_input <- count_df(mer_count_SSC.4_input)
count_mer_HC.1_input <- count_df(mer_count_HC.1_input)
count_mer_HC.2_input <- count_df(mer_count_HC.2_input)
count_mer_HC.3_input <- count_df(mer_count_HC.3_input)
count_mer_HC.4_input <- count_df(mer_count_HC.4_input)

count_mer_cfChIP <- count_mer_NAC.1_cfChIP %>% 
    left_join(count_mer_NAC.2_cfChIP, by = "res") %>% 
    left_join(count_mer_NAC.3_cfChIP, by = "res") %>% 
    left_join(count_mer_NAC.4_cfChIP, by = "res") %>% 
    left_join(count_mer_NSC.1_cfChIP, by = "res") %>% 
    left_join(count_mer_NSC.2_cfChIP, by = "res") %>% 
    left_join(count_mer_NSC.3_cfChIP, by = "res") %>% 
    left_join(count_mer_NSC.4_cfChIP, by = "res") %>% 
    left_join(count_mer_SSC.1_cfChIP, by = "res") %>% 
    left_join(count_mer_SSC.2_cfChIP, by = "res") %>% 
    left_join(count_mer_SSC.3_cfChIP, by = "res") %>% 
    left_join(count_mer_SSC.4_cfChIP, by = "res") %>% 
    left_join(count_mer_HC.1_cfChIP, by = "res") %>%
    left_join(count_mer_HC.2_cfChIP, by = "res") %>%
    left_join(count_mer_HC.3_cfChIP, by = "res") %>%
    left_join(count_mer_HC.4_cfChIP, by = "res")
count_mer_input <- count_mer_NAC.1_input %>% 
    left_join(count_mer_NAC.2_input, by = "res") %>% 
    left_join(count_mer_NAC.3_input, by = "res") %>% 
    left_join(count_mer_NAC.4_input, by = "res") %>% 
    left_join(count_mer_NSC.1_input, by = "res") %>% 
    left_join(count_mer_NSC.2_input, by = "res") %>% 
    left_join(count_mer_NSC.3_input, by = "res") %>% 
    left_join(count_mer_NSC.4_input, by = "res") %>% 
    left_join(count_mer_SSC.1_input, by = "res") %>% 
    left_join(count_mer_SSC.2_input, by = "res") %>% 
    left_join(count_mer_SSC.3_input, by = "res") %>% 
    left_join(count_mer_SSC.4_input, by = "res") %>% 
    left_join(count_mer_HC.1_input, by = "res") %>%
    left_join(count_mer_HC.2_input, by = "res") %>%
    left_join(count_mer_HC.3_input, by = "res") %>%
    left_join(count_mer_HC.4_input, by = "res")
colnames(count_mer_cfChIP) <- c("motif", 
                               "NAC.1_cfChIP",
                               "NAC.2_cfChIP",
                               "NAC.3_cfChIP",
                               "NAC.4_cfChIP",
                               "NSC.1_cfChIP",
                               "NSC.2_cfChIP",
                               "NSC.3_cfChIP",
                               "NSC.4_cfChIP",
                               "SSC.1_cfChIP",
                               "SSC.2_cfChIP",
                               "SSC.3_cfChIP",
                               "SSC.4_cfChIP",
                               "HC.1_cfChIP",
                               "HC.2_cfChIP",
                               "HC.3_cfChIP",
                               "HC.4_cfChIP")
colnames(count_mer_input) <- c("motif", 
                              "NAC.1_input",
                              "NAC.2_input",
                              "NAC.3_input",
                              "NAC.4_input",
                              "NSC.1_input",
                              "NSC.2_input",
                              "NSC.3_input",
                              "NSC.4_input",
                              "SSC.1_input",
                              "SSC.2_input",
                              "SSC.3_input",
                              "SSC.4_input",
                              "HC.1_input",
                              "HC.2_input",
                              "HC.3_input",
                              "HC.4_input")


count_mer_input <- read.table("fragment end motif counts input.txt", header = T) %>% select(-X)
count_mer_cfChIP <- read.table("fragment end motif counts cfChIP.txt", header = T) %>% select(-X)

dif_motif_prop <- function(x,y){
    library(dplyr)
    rownames(x) <- x$motif
    rownames(y) <- y$motif
    x1 <- x %>% dplyr::select(-motif)
    y1 <- y %>% dplyr::select(-motif)
    x2 <- as.data.frame(t(x1))
    y2 <- as.data.frame(t(y1))
    p <- c()
    se_x <- c()
    se_y <- c()
    avg_x <- c()
    avg_y <- c()
    diff <- c()
    fold <- c()
    for (i in 1:length(x2)){
        g1 <- x2[,i]
        g2 <- y2[,i]
        se_x[i] <- parameters::standard_error(g1)
        se_y[i] <- parameters::standard_error(g2)
        avg_x[i] <- mean(g1)
        avg_y[i] <- mean(g2)
        p[i] <- t.test(g1,g2, paired = T, alternative = "two.sided")$p.value
        diff[i] <- mean(g2-g1)
        fold[i] <- mean(g2/g1)
    }
    q <- p.adjust(p, method = "fdr")
    kk <- data.frame(motif = colnames(x2),
                     avg_x = avg_x,
                     se_x = se_x,
                     avg_y = avg_y,
                     se_y = se_y,
                     p.value = p,
                     q.value = q,
                     mean_diff = diff,
                     FC = fold)
    kk <- kk[order(-abs(kk$mean_diff)),]
    return(kk)
}
prop_mer_input_cancer <- prop_mer_input[!grepl("HC.",colnames(prop_mer_input))]
prop_mer_cfChIP_cancer <- prop_mer_cfChIP[!grepl("HC.",colnames(prop_mer_cfChIP))]
prop_mer_input_healthy <- prop_mer_input[grepl("HC.",colnames(prop_mer_input))]
prop_mer_cfChIP_healthy <- prop_mer_cfChIP[grepl("HC.",colnames(prop_mer_cfChIP))]
prop_mer_input_healthy$motif <- prop_mer_input$motif
prop_mer_cfChIP_healthy$motif <- prop_mer_cfChIP$motif
dif_motif_prop(prop_mer_input_cancer, prop_mer_cfChIP_cancer)

volcano_motif <- function(x,y,z){
    library(ggplot2)
    library(ggrepel)
    p <- -log10(z$q.value)
    z$Log10p <- p
    z <- z %>% dplyr::filter(!is.na(q.value))
    tit <- paste("End motif of", x, "compared to",y)
    gg <- ggplot(data = z, aes(x = FC, y = Log10p,
                               fill = FC))+
        geom_point(shape = 21,
                   stroke = 0.5,
                   size = 4)+
        scale_fill_gradient2(low = "#ffa10c", high = "#6a00fc",
                             mid = "grey",
                             name = "FC", midpoint = 1)+
        xlab("FC")+
        ylab(expression(paste("-",log[10]("q-value"), sep = "")))+
        labs(title = tit)+
        geom_label(x = 1.2, y = 0.1, label = paste(y),
                   color = "#6a00fc", label.size = 0, size = 5,
                   fill="white")+
        geom_label(x = (min(z$FC)+0.1), y = 0.1, label = paste(x),
                   color = "#ffa10c", label.size = 0, size = 5,
                   fill="white")+
        geom_hline(yintercept = c(-log10(0.05),-log10(0.01),-log10(0.001)),
                   linetype = c("dashed","dashed","dashed"),
                   colour = c("black", "black", "black"),
                   size = c(1,1,1))+
        geom_vline(xintercept = 1,
                   linetype = "solid",
                   colour = "black",
                   size = 1.5)+
        theme_bw(base_size = 17)+
        geom_label_repel(
            aes(label=ifelse(Log10p > -log10(0.05), as.character(motif),"")),
            segment.color="black",
            color="black",
            nudge_x = 0,
            label.size = NA,
            size = 2.5,
            fill = alpha(c("white"),0),
            parse = F,
            max.overlaps = 100)+
        th
    return(gg)
}
volcano_motif("Input cancer", "cfChIP cancer", dif_motif_prop(prop_mer_input_cancer, prop_mer_cfChIP_cancer))
volcano_motif("Input healthy", "cfChIP healthy", dif_motif_prop(prop_mer_input_healthy, prop_mer_cfChIP_healthy))
volcano_motif("Input", "cfChIP", dif_motif_prop(prop_mer_input, prop_mer_cfChIP))
x #Proportion of fragment end motifs for sample 1
y #Proportion of fragment end motifs for sample 2
p #Name of sample 1
q #Name of sample 2
s #seed for set.seed()
n #Number of neighbors
umap_moitf <- function(x,y, p, q, s, n){
    library(umap)
    set.seed(s)
    rownames(x) <- x$motif
    rownames(y) <- y$motif
    x1 <- x %>% dplyr::select(-motif)
    y1 <- y %>% dplyr::select(-motif)
    x2 <- as.data.frame(t(x1))
    y2 <- as.data.frame(t(y1))
    df <- rbind(x2,y2)
    umap_res <- umap(df, n_neighbors = n)
    umap_df <- as.data.frame(umap_res$layout)
    umap_df[,3] <- c(rep(p, nrow(x2)), rep(q, nrow(y2)))
    umap_df[,4] <- c(rep(c(rep("Cancer", 12), rep("Healthy", 4)),2))
    umap_df[,5] <- paste(umap_df$V4, umap_df$V3)
    gg <- ggplot(data = umap_df, aes(x = V1, y = V2, color = V5,
                                     shape = V5, group = V5))+
        geom_point(size = 5)+
        labs(title = "Input and cfChIP",
             x = "UMAP-1", y = "UMAP-2")+
        scale_color_manual(name = "Sample",
                           labels = c("Cancer input", "Cancer cfChIP",
                                      "Healthy input", "Healthy cfChIP"),
                           values = c("#6a00fc","#ffa10c","#6a00fc","#ffa10c"))+
        scale_shape_manual(name = "Sample",
                           labels = c("Cancer input", "Cancer cfChIP",
                                      "Healthy input","Healthy cfChIP"),
                           values = c(19,19,17,17))+
        theme_bw()+
        th
    return(gg)
}
umap_moitf(prop_mer_input, prop_mer_cfChIP,
           "input", "cfChIP", 15, 15)
gc()

x #name of cfChIP BAM file
y #Granges object returned by gr() for the 197 AVENIO target genes
z #.txt file of the AVENIO genes with the number of bases sequenced in each gene. based on data.frame returned by coverages()
i #name of input BAM file
h #data.frame object returned by healthy() for normilization to healthy. Otherwise NULL
m #number of bases to be evaluated in each fragment end

end_motif_active_inactive <- function(x,y,z,i,h,m){
    or_wd <- getwd()
    setwd("D:/Healthy cfChIP/Deduped")
    c <- gene_count(x,y)
    setwd("C:/Users/Christoffer/OneDrive/1PhD/R files/Important excel and text documents")
    e <- e.score(c,z,y,h)
    top15 <- e$genes[(length(e$genes)-14):length(e$genes)]
    bottom15 <- e$genes[1:15]
    which <- grs[elementMetadata(y)[,1] %in% c(top15)]
    index1 <- match(top15,elementMetadata(which)[,1])
    which <- which[index1]
    what <- c("seq", "strand")
    param <- ScanBamParam(which = which, what = what)
    setwd("D:/Healthy input/PosDeduped")
    bam <- scanBam(i, param = param)
    seqs_top15 <- unname(unlist(lapply(bam,function(x) as.data.frame(x$seq))))
    strand_top15 <- unname(unlist(lapply(bam,function(x) as.data.frame(x$strand))))
    df_top15 <- data.frame(seq = seqs_top15, 
                           strand = strand_top15)
    df_top15 <- na.omit(df_top15)
    which <- grs[elementMetadata(y)[,1] %in% c(bottom15)]
    index1 <- match(bottom15,elementMetadata(which)[,1])
    which <- which[index1]
    what <- c("seq", "strand")
    param <- ScanBamParam(which = which, what = what)
    bam <- scanBam(i, param = param)
    seqs_bottom15 <- unname(unlist(lapply(bam,function(x) as.data.frame(x$seq))))
    strand_bottom15 <- unname(unlist(lapply(bam,function(x) as.data.frame(x$strand))))
    df_bottom15 <- data.frame(seq = seqs_bottom15, 
                              strand = strand_bottom15)
    df_bottom15 <- na.omit(df_bottom15)
    mer_fun <- function(def){
        if (def[2] == "+"){
            seqs <- substring(def[1], 1, m)
        }
        if (def[2] == "-"){
            seqs <- substring(def[1],(nchar(def[1])-m+1) , nchar(def[1]))
        }
        if (grepl("N", seqs)){
            seqs <- NA
        }
        return(seqs)
    }
    top15_mers <- apply(df_top15, MARGIN = 1, FUN = mer_fun)
    top15_mers <- top15_mers[!is.na(top15_mers)]
    top15_mer_df <- as.data.frame.table(proportions(table(top15_mers)))
    bottom15_mers <- apply(df_bottom15, MARGIN = 1, FUN = mer_fun)
    bottom15_mers <- bottom15_mers[!is.na(bottom15_mers)]
    bottom15_mer_df <- as.data.frame.table(proportions(table(bottom15_mers)))
    df_res <- top15_mer_df %>% left_join(bottom15_mer_df, by = c("top15_mers" = "bottom15_mers"))
    colnames(df_res) <- c("motif", "Active", "Inactive")
    setwd(or_wd)
    return(df_res)
}
?unlist
NAC.1_active_inactive_df <- end_motif_active_inactive("A-1279-cfChIP.bam",
                                                      grs,
                                                      "Coverage of AVENIO genes.txt",
                                                      "PosDeduped-A-1279-input.bam",
                                                      healthy_reads,3)
NAC.2_active_inactive_df <- end_motif_active_inactive("B-1288-cfChIP.bam",
                                                      grs,
                                                      "Coverage of AVENIO genes.txt",
                                                      "PosDeduped-B-1288-input.bam",
                                                      healthy_reads,3)
NAC.3_active_inactive_df <- end_motif_active_inactive("C-1475-cfChIP.bam",
                                                      grs,
                                                      "Coverage of AVENIO genes.txt",
                                                      "PosDeduped-C-1475-input.bam",
                                                      healthy_reads,3)

NAC.4_active_inactive_df <- end_motif_active_inactive("D-1578-cfChIP.bam",
                                                      grs,
                                                      "Coverage of AVENIO genes.txt",
                                                      "PosDeduped-D-1578-input.bam",
                                                      healthy_reads,3)

NSC.1_active_inactive_df <- end_motif_active_inactive("E-439-cfChIP.bam",
                                                      grs,
                                                      "Coverage of AVENIO genes.txt",
                                                      "PosDeduped-E-439-input.bam",
                                                      healthy_reads,3)

NSC.2_active_inactive_df <- end_motif_active_inactive("F-1449-cfChIP.bam",
                                                      grs,
                                                      "Coverage of AVENIO genes.txt",
                                                      "PosDeduped-F-1449-input.bam",
                                                      healthy_reads,3)

NSC.3_active_inactive_df <- end_motif_active_inactive("I-645-cfChIP.bam",
                                                      grs,
                                                      "Coverage of AVENIO genes.txt",
                                                      "PosDeduped-I-645-input.bam",
                                                      healthy_reads,3)

NSC.4_active_inactive_df <- end_motif_active_inactive("J-1663-cfChIP.bam",
                                                      grs,
                                                      "Coverage of AVENIO genes.txt",
                                                      "PosDeduped-J-1663-input.bam",
                                                      healthy_reads,3)

SSC.1_active_inactive_df <- end_motif_active_inactive("G-514-cfChIP.bam",
                                                      grs,
                                                      "Coverage of AVENIO genes.txt",
                                                      "PosDeduped-G-514-input.bam",
                                                      healthy_reads,3)

SSC.2_active_inactive_df <- end_motif_active_inactive("H-1169-cfChIP.bam",
                                                      grs,
                                                      "Coverage of AVENIO genes.txt",
                                                      "PosDeduped-H-1169-input.bam",
                                                      healthy_reads,3)

SSC.3_active_inactive_df <- end_motif_active_inactive("K-440-cfChIP.bam",
                                                      grs,
                                                      "Coverage of AVENIO genes.txt",
                                                      "PosDeduped-K-440-input.bam",
                                                      healthy_reads,3)

SSC.4_active_inactive_df <- end_motif_active_inactive("L-1100-cfChIP.bam",
                                                      grs,
                                                      "Coverage of AVENIO genes.txt",
                                                      "PosDeduped-L-1100-input.bam",
                                                      healthy_reads,3)
HC.1_active_inactive_df <- end_motif_active_inactive(x = "Rask_kontrol1_cfChIP.bam",
                                                      y = grs,
                                                      z = "Coverage of AVENIO genes.txt",
                                                      i = "PosDeduped-Rask_kontrol_1_input.bam",
                                                      h = NULL,
                                                      m = 3)
HC.2_active_inactive_df <- end_motif_active_inactive(x = "Rask_kontrol2_cfChIP.bam",
                                                     y = grs,
                                                     z = "Coverage of AVENIO genes.txt",
                                                     i = "PosDeduped-Rask_kontrol_2_input.bam",
                                                     h = NULL,
                                                     m = 3)
HC.3_active_inactive_df <- end_motif_active_inactive(x = "Rask_kontrol3_cfChIP.bam",
                                                     y = grs,
                                                     z = "Coverage of AVENIO genes.txt",
                                                     i = "PosDeduped-Rask_kontrol_3_input.bam",
                                                     h = NULL,
                                                     m = 3)
HC.4_active_inactive_df <- end_motif_active_inactive(x = "Rask_kontrol4_cfChIP.bam",
                                                     y = grs,
                                                     z = "Coverage of AVENIO genes.txt",
                                                     i = "PosDeduped-Rask_kontrol_4_input.bam",
                                                     h = NULL,
                                                     m = 3)

collected_active_inactive_df <- NAC.1_active_inactive_df %>% 
    left_join(NAC.2_active_inactive_df, by = "motif") %>%
    left_join(NAC.3_active_inactive_df, by = "motif") %>%
    left_join(NAC.4_active_inactive_df, by = "motif") %>%
    left_join(NSC.1_active_inactive_df, by = "motif") %>%
    left_join(NSC.2_active_inactive_df, by = "motif") %>%
    left_join(NSC.3_active_inactive_df, by = "motif") %>%
    left_join(NSC.4_active_inactive_df, by = "motif") %>%
    left_join(SSC.1_active_inactive_df, by = "motif") %>%
    left_join(SSC.2_active_inactive_df, by = "motif") %>%
    left_join(SSC.3_active_inactive_df, by = "motif") %>%
    left_join(SSC.4_active_inactive_df, by = "motif")

colnames(collected_active_inactive_df) <- c("motif",
                                            "NAC.1_active", "NAC.1_inactive",
                                            "NAC.2_active", "NAC.2_inactive",
                                            "NAC.3_active", "NAC.3_inactive",
                                            "NAC.4_active", "NAC.4_inactive",
                                            "NSC.1_active", "NSC.1_inactive",
                                            "NSC.2_active", "NSC.2_inactive",
                                            "NSC.3_active", "NSC.3_inactive",
                                            "NSC.4_active", "NSC.4_inactive",
                                            "SSC.1_active", "SSC.1_inactive",
                                            "SSC.2_active", "SSC.2_inactive",
                                            "SSC.3_active", "SSC.3_inactive",
                                            "SSC.4_active", "SSC.4_inactive",
                                            "HC.1_active", "HC.1_inactive",
                                            "HC.2_active", "HC.2_inactive",
                                            "HC.3_active", "HC.3_inactive",
                                            "HC.4_active", "HC.4_inactive")

write.table(collected_active_inactive_df, file = "fragment end motif proportion of active and inactive fragments.txt",
            sep = "\t", col.names= T)
setwd("C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver")
collected_active_inactive_df <- read.table("fragment end motif proportion of active and inactive fragments.txt", header = T)


diff_active_inactive_end_motif <- function(x){
    library(dplyr)
    rownames(x) <- x$motif
    x1 <- x %>% dplyr::select(-motif)
    x_active <- x1[grepl("_active", colnames(x1))]
    x_inactive <- x1[grepl("_inactive", colnames(x1))]
    y2 <- as.data.frame(t(x_active))
    x2 <- as.data.frame(t(x_inactive))
    p <- c()
    se_x <- c()
    se_y <- c()
    avg_x <- c()
    avg_y <- c()
    diff <- c()
    fold <- c()
    for (i in 1:length(x2)){
        g1 <- x2[,i]
        g2 <- y2[,i]
        se_x[i] <- parameters::standard_error(g1)
        se_y[i] <- parameters::standard_error(g2)
        avg_x[i] <- mean(g1)
        avg_y[i] <- mean(g2)
        p[i] <- t.test(g1,g2, paired = T, alternative = "two.sided")$p.value
        diff[i] <- mean(g2-g1)
        fold[i] <- mean(g2/g1)
    }
    q <- p.adjust(p, method = "fdr")
    kk <- data.frame(motif = colnames(x2),
                     avg_x = avg_x,
                     se_x = se_x,
                     avg_y = avg_y,
                     se_y = se_y,
                     p.value = p,
                     q.value = q,
                     mean_diff = diff,
                     FC = fold)
    kk <- kk[order(-abs(kk$mean_diff)),]
    return(kk)
}

volcano_motif("Inactive genes", "Active genes", diff_active_inactive_end_motif(collected_active_inactive_df))


x #data.frame returned by diff_active_inactive_end_motif() eller dif_motif_prop()
q #value cutoff
f #Can be "negative" eller "positive" for om ente FC < 1 eller FC > 1 endemotiver skal undersøge
sequence_end_motif <- function(x, q, f, tit){
    library(ggseqlogo)
    sig_motif <- x %>% filter(q.value < q)
    if(f == "positive"){
        sig_motif <- sig_motif %>%  filter(FC > 1)
        end_motif <- sig_motif$motif
    }
    if(f == "negative"){
        sig_motif <- sig_motif %>%  filter(FC < 1)
        end_motif <- sig_motif$motif
    }
    return(end_motif)
}
sequence_list <- list(sequence_end_motif(dif_motif_prop(prop_mer_input, prop_mer_cfChIP), 0.05, "positive"),
                      sequence_end_motif(dif_motif_prop(prop_mer_input, prop_mer_cfChIP), 0.05, "negative"),
                      sequence_end_motif(diff_active_inactive_end_motif(collected_active_inactive_df), 0.05, "positive"),
                      sequence_end_motif(diff_active_inactive_end_motif(collected_active_inactive_df), 0.05, "negative"))
x #list of character vectors with significant sequences
y #character vector of names of samples
sequence_end_motif_plot <- function(x,y){
    library(ggseqlogo)
    names(x) <- y
    cs <- make_col_scheme(chars=c('A', 'T', 'C', 'G'), 
                          groups=c('Weak', 'Weak', 'Strong', 'Strong'), 
                          cols=c("#ffa10c", "#ffa10c", "#6a00fc", "#6a00fc"),
                          name = "Bond strengh")
    gg <- ggseqlogo(x, ncol=2, method = "prob",
                    col_scheme = cs)+
        theme_logo()+
        th+
        theme(axis.text.x = element_blank(),
              strip.text = element_text(face = "bold",
                                        size = 16))
    return(gg)
}
sequence_end_motif_plot(sequence_list, c("cfChIP", "Input",
                                         "Active genes", "Inactive genes"))
   
sequence_list

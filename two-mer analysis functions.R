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
x # BAM file returned by bamfile()
m #number of bases to be evaluated in each fragment end
mer_count <- function(x,m){
    positve <- x[strand(x) == "+"]
    negative <- x[strand(x) == "-"]
    rev_strands <- DNAStringSet(stringi::stri_reverse(mcols(negative)$seq))
    positve_counts <- nucleotideFrequencyAt(
        x = mcols(positve)$seq, at = c(1,2,3), as.array = F)
    negative_counts <- nucleotideFrequencyAt(
        x = rev_strands,
        at = c(1,2,3), as.array = F)
    return(positve_counts+negative_counts)
    
}

x # table returned by mer_count()
prop_df <- function(x){
    prop <- proportions(x)
    df <- data.frame(res = names(prop),
                     Freq = unname(prop))
    return(df)
}

x # table returned by mer_count()
count_df <- function(x){
    df <- data.frame(res = names(x),
                     Freq = unname(x))
    return(df) 
}
x #prop_mer data.frame for group 1
x #prop_mer data.frame for group 2

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

x #Name of group 1
y #Name of group 2
z #data.frame returned by dif_motif_prop() eller diff_active_inactive_end_motif()
volcano_motif <- function(x,y,z){
    library(ggplot2)
    library(ggrepel)
    p <- -log10(z$q.value)
    z$Log10p <- p
    z <- z %>% dplyr::filter(!is.na(q.value))
    tit <- paste("End motif of", x, "compared to",y)
    if(min(z$FC)<0.8){
        xmin = log2(min(z$FC)-0.1)
        x_lab1 = log2(1-((1-min(z$FC))/2))
    }
    else{
        xmin = log2(0.8)
        x_lab1 = log2(0.9)
    }
    if(max(z$FC)>1.2){
        xmax = log2(max(z$FC)+0.1)
        x_lab2 = log2(1+((max(z$FC)-1)/2))
    }
    else{
        xmax = log2(1.2)
        x_lab2 = log2(1.1)
    }
    gg <- ggplot(data = z, aes(x = log2(FC), y = Log10p,
                               fill = log2(FC)))+
        geom_label(x = x_lab2, y = 0.1, label = paste(y),
                   color = "#6a00fc", label.size = 0, size = 5,
                   fill="white")+
        geom_label(x = x_lab1, y = 0.1, label = paste(x),
                   color = "#ffa10c", label.size = 0, size = 5,
                   fill="white")+
        geom_point(shape = 21,
                   stroke = 0.5,
                   size = 4)+
        scale_fill_gradient2(low = "#ffa10c", high = "#6a00fc",
                             mid = "grey",
                             name = expression(bold(paste(log[2]("FC"), sep = ""))),
                             midpoint = 0)+
        xlab(expression(bold(paste(log[2]("FC"), sep = ""))))+
        ylab(expression(bold(paste("-",log[10]("q-value"), sep = ""))))+
        labs(title = tit)+
        geom_hline(yintercept = c(-log10(0.05),-log10(0.01),-log10(0.001)),
                   linetype = c("dashed","dashed","dashed"),
                   colour = c("black", "black", "black"),
                   size = c(1,1,1))+
        geom_vline(xintercept = 0,
                   linetype = "solid",
                   colour = "black",
                   size = 1.5)+
        scale_x_continuous(limits = c(xmin, xmax))+
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

x #Proportion of fragment end motifs for sample 1
y #Proportion of fragment end motifs for sample 2
p #Name of sample 1
q #Name of sample 2
s #seed for set.seed()
n #Number of neighbors
t #Type of plot. For healthy and cancer type "both" (default), for cancer = "cancer", for healthy = "healthy"
umap_moitf <- function(x,y, p, q, s, n, t = "both"){
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
    if(t == "both"){
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
    }
    if(t == "cancer"){
        gg <- ggplot(data = umap_df, aes(x = V1, y = V2, color = V3))+
            geom_point(size = 5)+
            labs(title = "Input and cfChIP",
                 x = "UMAP-1", y = "UMAP-2")+
            scale_color_manual(name = "Sample",
                               labels = c("Cancer input", "Cancer cfChIP"),
                               values = c("#6a00fc","#ffa10c"))+
            theme_bw()+
            th
    }
    if(t == "healthy"){
        gg <- ggplot(data = umap_df, aes(x = V1, y = V2, color = V3))+
            geom_point(size = 5)+
            labs(title = "Input and cfChIP",
                 x = "UMAP-1", y = "UMAP-2")+
            scale_color_manual(name = "Sample",
                               labels = c("Healthy input", "Healthy cfChIP"),
                               values = c("#6a00fc","#ffa10c"))+
            theme_bw()+
            th
    }
    
    return(gg)
}

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

x #prop_mer data.frame with paired active and inactive gene fragment end motifs
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
        theme_bw()+
        th+
        theme(axis.text.x = element_blank(),
              strip.text = element_text(face = "bold",
                                        size = 16),
              legend.position = "bottom",
              strip.background = element_rect(fill = "white"),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    return(gg)
}

x #bam.file returned by bamfile()
y #genomic ranges of AVENIO genes returned by gr()
z #character vector of motifs to get fraction of
motif_fraction <- function(x, y, z){
    fractions <- c()
    y <- y[mcols(y)$SYMBOL != "DLGAP2"]
    for (i in 1:length(y)){
        print(round(i/length(y)*100,2))
        sub <- subsetByOverlaps(x, y[i],ignore.strand=TRUE)
        tab <- mer_count(sub, 3)
        df <- as.data.frame.table(tab)
        sig_df <- df %>% filter(res %in% z)
        res <- sum(sig_df$Freq)/sum(df$Freq) * 100
        fractions[i] <- res
    }
    res_df <- data.frame(genes = mcols(y)$SYMBOL,
                         fractions = fractions)
    return(res_df)
}
                                        
write.table(NAC.1_active_fraction, file = "NAC.1 active motif fraction.txt", 
            sep = "\t",col.names = T)

x #name on .txt file containing each gene, the number of times each gene is sequenced and the distance each target is from TSS
y #cutoff of distance from TSS which is too long
badgene <- function(x,y = 25){
    library(dplyr)
    df <- read.table(x,header = T)
    df$Relative_distance <- as.numeric(df$Relative_distance)
    sr <- df %>% filter(Regions == 1) %>% filter(Relative_distance < y)
    gene_names <- sr$SYMBOL
    return(gene_names)
}

multiple_region <- function(x){
    df <- read.table(x, header = T)
    df <- df %>% filter(Regions > 1)
    return(unique(df$SYMBOL))
}
x #data.frame returned by motif_fraction()
y #data.frame returned by e.score()
z #Name of sample
b #optional: if not NULL then use character vector returned by badgene()
r #optional: if not NULL then use character vector returned by multiple_region()
motif_fraction_vs_enrichment <- function(x,y,z,b = NULL, r = NULL){
    `%ni%` <- Negate(`%in%`)
    if (is.null(b)){
        x1 <- x
        y1 <- y
    }
    else{
        x1 <- x %>% filter(genes %ni% b)
        y1 <- y %>% filter(genes %ni% b)
    }
    if (is.null(r)){
        x2 <- x1
        y2 <- y1
        tit = "all genes"
        
    }
    else{
        x2 <- x1 %>% filter(genes %in% r)
        y2 <- y1 %>% filter(genes %in% r)
        tit = "multi region genes"
    }
    df <- x2 %>% left_join(y2, by = "genes")
    colnames(df) <- c("genes", "fractions", "enrichment")
    linear <- lm(enrichment~fractions, data=df)
    coef <- as.numeric(linear$coefficients[2])
    res <- cor.test(df$enrichment, df$fractions, method = "spearman")
    linear <- lm(enrichment~fractions, data=df)
    p_values <- res$p.value
    if (p_values < 0.0001){
        p <- ", P < 0.0001"
    }
    else{
        p <- paste(", P =", round(p_values,4))
    }
    rhos <- as.numeric(res$estimate)
    gg <- df %>% ggplot(aes(x = fractions,
                      y = enrichment))+
        geom_point(size = 4, color = "#ffa10c")+
        theme_bw()+
        geom_smooth(data = df,
                    aes(x=fractions, y = enrichment),
                    method = "lm", se = F, color = "black")+
        labs(title = paste(z, "for", tit), 
             x = "Fraction of active motifs", y = "cfChIP enrichment",
             subtitle = paste("Spearman's rho =", round(rhos,3), 
                              as.character(p),", alpha = ",round(coef,2)))+
        th
    
    return(gg)
}
x #correlation matrix with samples, rhos, p and alpha values
z #Title of plot
matrix_plot <- function(x,z){
    library(tidyr)
    options(scipen = 999)
    plot.dat <- gather(x, key = "type", value = "values", -sample)
    gg <- ggplot(data = plot.dat,aes(x = factor(type,level = c("rho", "p.value","alpha")),
                                                y = sample, label = round(values,4)))+
        geom_tile(data = plot.dat %>% filter(type == "rho"),
                  aes(fill = values))+
        geom_text(col = "black")+
        scale_fill_gradient("Rho", limits = c(0, 1), 
                             low = "white", high = "#1B7837")+
        scale_x_discrete(limits = c("rho", "p.value","alpha"))+
        ggnewscale::new_scale_fill()+
        geom_tile(data = plot.dat %>% filter(type == "p.value"),
                  aes(fill = values))+
        geom_text(col = "black")+
        scale_fill_gradient("p.value", limits = c(0,1), 
                             low = "#762A83", high = "white")+
        ggnewscale::new_scale_fill()+
        geom_tile(data = plot.dat %>% filter(type == "alpha"),
                  aes(fill = values))+
        geom_text(col = "black")+
        scale_fill_gradient("Alpha", 
                            low = "white", high = "tomato3")+
        theme_bw()+
        labs(x = "", y = "Cancer sample",
             title = z)+
        th
        return(gg) 
 
    
}
x #BAM file returned by bamfile()
y #fragment length cut-off, default = 150
z #TRUE to get fragments shorter than or equal to cutoff, whereas FALSE gives fragments longer than cutoff. Default is TRUE
short_long_fragments <- function(x, y = 150,z = TRUE){
    mcols(x)$abs_size <- abs(mcols(x)$isize)
    x <- x[!is.na(mcols(x)$abs_size)]
    if(isTRUE(z)){
        res <- x[mcols(x)$abs_size <= y]
        return(res)
    }
    if(isFALSE(z)){
        res <- x[mcols(x)$abs_size > y]
        return(res)
    }
}


x #BAM file returned by bamfile()
y #fragment length cut-off, default = 150
m #number of bases to be evaluated in each fragment end
short_long_fragments_motifs <- function(x,y,m){
    short <- short_long_fragments(x,y,z = TRUE)
    long <- short_long_fragments(x,y,z = FALSE)
    short_table <- mer_count(short, m = m)
    long_table <- mer_count(long, m = m)
    short_df <- prop_df(short_table)
    long_df <- prop_df(long_table)
    df <- long_df %>% left_join(short_df, by = "res")
    colnames(df) <- c("motif", "long", "short")
    return(df)
}

x #list of character vectors with significant sequences
y #character vector of names of samples
tit #title of plot
motif_venn <- function(x,y, tit){
    library(ggVennDiagram)
    library(ggpubr)
    gg <- ggVennDiagram(x, set_color = "black", 
                             edge_lty = 1, edge_size = 2, 
                             label_percent_digit = 2,
                             category.names = y,
                             color = "black", label_size = 5, label = "count")+
        scale_fill_gradient(low = "white", high = "#1B7837", name = "Motifs")+
        theme_bw(base_size = 12)+
        labs(title = tit,
             x = "",
             y = "")+
        th+
        theme(
            axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.line = element_blank(),
            axis.title = element_blank(),
            legend.background = element_rect(color = NA),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            axis.ticks = element_blank())
    return(gg)
}
e #enrichment dataframe
n #Name of sample
x #input BAM file returned by bamfile()
g #Granges object returned by gr()
y #fragment length cut-off, default = 150
m #number of bases to be evaluated in each fragment end. Default = 3
active_short_long_fragments_motifs <- function(e, n, x, g, y = 150, m = 3){
    e_sample <- e %>% select(c("genes", n))
    bottom_15 <- e_sample$genes[order(e_sample[,2])][1:15]
    top_15 <- e_sample$genes[order(-e_sample[,2])][1:15]
    g_bottom15 <- g[mcols(g)$SYMBOL %in% bottom_15]
    g_top15 <- g[mcols(g)$SYMBOL %in% top_15]
    sub_bottom <- subsetByOverlaps(x, g_bottom15,ignore.strand=TRUE)
    sub_top <- subsetByOverlaps(x, g_top15,ignore.strand=TRUE)
    bottom_fractions <- short_long_fragments_motifs(sub_bottom, y = y, m=m)
    top_fractons <- short_long_fragments_motifs(sub_top, y = y, m=m)
    df <- bottom_fractions %>% left_join(top_fractons, by = "motif")
    colnames(df) <- c("motif", "bottom_long", "bottom_short", "top_long", "top_short")
    return(df)
}

x #list of data.frames returned by e.score()
y #vector of names for each data.frame.
z #vector of names for groups
bind_samples <- function(x,y,z){
    library(dplyr)
    len <- length(x)
    lens <- c()
    for (i in 1:len){
        tran <- transpose(x[[i]], y[i])
        x[[i]] <- tran
        lens[i] <- length(tran)
    }
    short <- which.min(lens)
    std <- colnames(x[[short]])
    for (i in 1:len){
        x[[i]] <- x[[i]] %>% dplyr::select(all_of(std))
    }
    res <- bind_rows(x)
    res <- as.data.frame(sapply(res, as.numeric))
    rownames(res) <- y
    res$groups <- z
    return(res)
}

x #data.frame returned by e.score()
y #name of sample
transpose <- function(x,y){
    x1 <- as.data.frame(t(x))
    colnames(x1) <- x1[1,1:length(x1)]
    x1 <- x1[-1,]
    rownames(x1) <- y
    return(x1)
}

x #list of data.frames returned by e.score()
average_enrichment <- function(x){
    len <- length(x)
    std <- x[[1]]$genes
    df <- data.frame(genes = std)
    for (i in 1:len){
        x[[i]] <- x[[i]][match(std,x[[i]]$genes),]
        df[ , ncol(df) + 1] <- x[[i]]$enrichment
        colnames(df)[ncol(df)] <- paste0(i)
    }
    avg <- c()
    for (i in 1:length(df$genes)){
        avg[i] <- mean(as.numeric(df[i,2:length(df)]))
    }
    df$enrichment <- avg
    df <- df[order(df$enrichment),]
    return(df)
}

x #Average Enrichment data.frame return from average_enrichmet() for sample 1
y #Average Enrichment data.frame return from average_enrichmet() for sample 2
conf <- function(x,y,b = NULL,v,w){
    library(dplyr)
    `%ni%` <- Negate(`%in%`)
    if(length(x$genes) == length(y$genes)){
        y = y
        x = x
    }
    else{
        if (length(x$genes) > length(y$genes)){
            x <- x %>% filter(genes %in% y$genes)
        }
        if (length(y$genes) > length(x$genes)){
            y <- y %>% filter(genes %in% x$genes)
        }
    }
    y <- y[(match(x$genes, y$genes)),]
    x_df <- x
    y_df <- y
    rownames(x_df) <- x_df$genes
    rownames(y_df) <- y_df$genes
    x_df <- x_df %>% select(-genes)
    y_df <- y_df %>% select(-genes)
    group1_t <- as.data.frame(t(x_df))
    group2_t <- as.data.frame(t(y_df))
    p <- c()
    group1 <- c()
    group2 <- c()
    se_group1 <- c()
    se_group2 <- c()
    conf_lower <- c()
    conf_upper <- c()
    clean <- c()
    for (i in 1:(length(group1_t))){
        g1 <- group1_t[,i]
        g2 <- group2_t[,i]
        res <- t.test(log2(g1),log2(g2), paired = F, alternative = "two.sided")
        p[i] <- res$p.value
        group1[i] <- mean(g1)
        group2[i] <- mean(g2)
        conf_lower[i] <- res$conf.int[1]
        conf_upper[i] <- res$conf.int[2]
        clean[i] <- paste0("[",round(res$conf.int[1],4), " - ", round(res$conf.int[2],4),"]")
    }
    q <- p.adjust(p, method = "fdr")
    kk <- data.frame(gene = colnames(group1_t),
                     group1_name = round(group1,0),
                     group2_name = round(group2,0),
                     conf_lower = conf_lower,
                     conf_upper = conf_upper,
                     Log2FC = round(log2(group1/group2),4),
                     clean = clean)
    kk <- kk %>% filter(conf_upper > 0 & conf_lower > 0 | conf_upper < 0 & conf_lower < 0)
    kk$output <- paste(kk$Log2FC, kk$clean)
    if(!is.null(b)){
        kk <- kk[kk$gene %ni% b,]
    }
    kk <- kk %>% select(gene, group1_name, group2_name, Log2FC, output)
    colnames(kk) <- c("genes", paste(v), paste(w), "Log2FC", "output")
    return(kk)
}

x #data.frame with genes of interest returned by conf()
g #Granges object returned by gr()
i #input BAM file to be analyzed, returned by bamfile()
diff_gene_reads <- function(x,g,i){
    gene_g <- g[mcols(g)$SYMBOL %in% x$genes]
    sub <- subsetByOverlaps(i, gene_g,ignore.strand=TRUE)
    return(sub)
}

x #input BAM file to be analyzed, returned by bamfile()
y #cutoff for short fragments
sub150_fraction <- function(x,y){
    mcols(x)$abs_size <- abs(mcols(x)$isize)
    res <- mcols(x)$abs_size<y
    res <- res[!is.na(res)]
    return(mean(res))
}
x #data.frame with genes of interest returned by conf()
g #Granges object returned by gr()
l #list of input BAM file to be analyzed, returned by bamfile()
y #cutoff for short fragments
m #number of bases to be evaluated in each fragment end
diff_gene_epigenetic_analysis <- function(x,g,l,y,m){
    fragment_res <- c()
    for(j in 1:length(l)){
        print(j)
        sub <- diff_gene_reads(x = x, g = g, i = l[[j]])
        fragment_res[j] <- sub150_fraction(sub, y = y)
        tab <- mer_count(sub,m=m)
        if(j == 1){
            df <- as.data.frame.table(proportions(tab))
        }
        else{
            df <- df %>% 
                left_join(as.data.frame.table(proportions(tab)), by = "res")
        }
        
    }
    fragment_df <- data.frame(sub150 = fragment_res)
    return(list(fragment_df,df))
}
cancer_input <- list(NAC.1_input, NAC.2_input,
                     NAC.3_input, NAC.4_input,
                     NSC.1_input, NSC.2_input,
                     NSC.3_input, NSC.4_input,
                     SSC.1_input, SSC.2_input,
                     SSC.3_input, SSC.4_input)
healthy_input <- list(HC.1_input, HC.2_input,
                      HC.3_input, HC.4_input)

healthy_upregulated_cancer_files <- diff_gene_epigenetic_analysis(healthy_upregulated, 
                                                                  grs, 
                                                                  cancer_input,
                                                                  150,
                                                                  3)
healthy_upregulated_healthy_files <- diff_gene_epigenetic_analysis(healthy_upregulated, 
                                                                   grs, 
                                                                   healthy_input,
                                                                   150,
                                                                   3)
cancer_upregulated_cancer_files <- diff_gene_epigenetic_analysis(cancer_upregulated, 
                                                                 grs, 
                                                                 cancer_input,
                                                                 150,
                                                                 3)
cancer_upregulated_healthy_files <- diff_gene_epigenetic_analysis(cancer_upregulated, 
                                                                  grs, 
                                                                  healthy_input,
                                                                  150,
                                                                  3)
colnames(cancer_upregulated_cancer_files[[2]])[1] <- "motif"
colnames(cancer_upregulated_healthy_files[[2]])[1] <- "motif"
colnames(healthy_upregulated_cancer_files[[2]])[1] <- "motif"
colnames(healthy_upregulated_healthy_files[[2]])[1] <- "motif"
x #data.frame returned by diff_gene_epigenetic_analysis()[[2]] for sample 1
y #data.frame returned by diff_gene_epigenetic_analysis()[[2]] for sample 2
diff_gene_motiv_analysis <- function(x, y){
    library(dplyr)
    rownames(x) <- x$motif
    rownames(y) <- y$motif
    x1 <- x %>% dplyr::select(-"motif")
    y1 <- y %>% dplyr::select(-"motif")
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
        p[i] <- t.test(g1,g2, paired = F, alternative = "two.sided")$p.value
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
volcano_motif("Healthy","Cancer",
              diff_gene_motiv_analysis(cancer_upregulated_healthy_files[[2]],cancer_upregulated_cancer_files[[2]]))
volcano_motif("Healthy","Cancer",
              diff_gene_motiv_analysis(healthy_upregulated_healthy_files[[2]],healthy_upregulated_cancer_files[[2]]))

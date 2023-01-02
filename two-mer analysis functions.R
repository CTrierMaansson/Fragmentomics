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
        what = c("mapq","isize","pos","cigar","seq"))
}
?readBamFileAsGRanges
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
    library(ggrepel)
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
            labs(title = paste(p, "and", q),
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
            labs(title = paste(p, "and", q),
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
            labs(title = paste(p, "and", q),
                 x = "UMAP-1", y = "UMAP-2")+
            scale_color_manual(name = "Sample",
                               labels = c("Healthy input", "Healthy cfChIP"),
                               values = c("#6a00fc","#ffa10c"))+
            theme_bw()+
            th
    }
    if(t == "Activity"){
        umap_df[,4] <- c(rep(c(rep("Cancer", 12), rep("Healthy", 4)),2))
        umap_df[,5] <- paste(umap_df$V4, umap_df$V3)
        umap_df[,6] <- rep(c(colnames(x)[-1],colnames(y)[-17]))
        gg <- ggplot(data = umap_df, aes(x = V1, y = V2, color = V5,
                                         shape = V5, group = V5))+
            geom_point(size = 5)+
            labs(title = paste(p, "and", q),
                 x = "UMAP-1", y = "UMAP-2")+
            scale_color_manual(name = "Sample",
                               labels = c("Cancer active", "Cancer inactive",
                                          "Healthy active", "Healthy inactive"),
                               values = c("#6a00fc","#ffa10c","#6a00fc","#ffa10c"))+
            scale_shape_manual(name = "Sample",
                               labels = c("Cancer active", "Cancer inactive",
                                          "Healthy active","Healthy inactive"),
                               values = c(19,19,17,17))+
            theme_bw()+
            geom_label_repel(
                aes(label=V6),
                segment.color="black",
                color="black",
                nudge_x = 0,
                label.size = NA,
                size = 2.5,
                fill = alpha(c("white"),0),
                parse = F,
                max.overlaps = 100)+
            th
    }
    
    return(gg)
}

x #enrichment data.frame 
y #Granges object returned by gr() for the 197 AVENIO target genes
z #name of sample to isolate from enrichment data.frame
i #name of input BAM file
m #number of bases to be evaluated in each fragment end
end_motif_active_inactive <- function(x,y,z,i,m){
    library(Rsamtools)
    df_sample <- x[c("genes",z)]
    df_sample <- df_sample[order(df_sample[,2]),]
    e <- df_sample
    return(e)
    top15 <- e$genes[(length(e$genes)-14):length(e$genes)]
    bottom15 <- e$genes[1:15]
    return(list(top15, bottom15))
    which <- grs[elementMetadata(y)[,1] %in% c(top15)]
    index1 <- match(top15,elementMetadata(which)[,1])
    which <- which[index1]
    what <- c("seq", "strand")
    param <- ScanBamParam(which = which, what = what)
    bam <- scanBam(i, param = param)
    seqs_top15 <- unname(unlist(lapply(bam,function(x) x$seq)))
    strand_top15 <- unname(unlist(lapply(bam,function(x) x$strand)))
    grang_top15 <- GRanges(seqnames = "chr1", ranges = IRanges(start = 1,
                                                               end = 2),
                           strand = strand_top15)
    mcols(grang_top15)$seq <- DNAStringSet(do.call(c,unlist(seqs_top15)))
    which <- grs[elementMetadata(y)[,1] %in% c(bottom15)]
    index1 <- match(bottom15,elementMetadata(which)[,1])
    which <- which[index1]
    what <- c("seq", "strand")
    param <- ScanBamParam(which = which, what = what)
    bam <- scanBam(i, param = param)
    seqs_bottom15 <- unname(unlist(lapply(bam,function(x) x$seq)))
    strand_bottom15 <- unname(unlist(lapply(bam,function(x) x$strand)))
    grang_bottom15 <- GRanges(seqnames = "chr1", ranges = IRanges(start = 1,
                                                               end = 2),
                           strand = strand_bottom15)
    mcols(grang_bottom15)$seq <- DNAStringSet(do.call(c,unlist(seqs_bottom15)))
    
    top15_mers <- mer_count(grang_top15)
    top15_mer_df <- prop_df(top15_mers)
    bottom15_mers <- mer_count(grang_bottom15)
    bottom15_mer_df <- prop_df(bottom15_mers)
    df_res <- top15_mer_df %>% left_join(bottom15_mer_df, by = "res")
    colnames(df_res) <- c("motif", "Active", "Inactive")
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
    library(Biostrings)
    fractions <- c()
    y <- y[mcols(y)$SYMBOL != "DLGAP2"]
    for (i in 1:length(y)){
        print(round(i/length(y)*100,2))
        sub <- subsetByOverlaps(x, y[i],ignore.strand=TRUE)
        tab <- mer_count(sub, 3)
        df <- prop_df(tab)
        sig_df <- df %>% filter(res %in% z)
        res <- sum(sig_df$Freq)
        fractions[i] <- res
    }
    res_df <- data.frame(genes = mcols(y)$SYMBOL,
                         fractions = fractions)
    return(res_df)
}
                                        

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
    df$fractions <- df$fractions * 100
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
             x = "Fraction of active motifs (%)", y = "cfChIP enrichment",
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

x #Data.frame of sample 1
y #data.frame of sample 2
z #Character vecotr of sample names
clusters <- function(x,y,z){
    library(plotly)
    rownames(x) <- x$motif
    rownames(y) <- y$motif
    x1 <- x %>% dplyr::select(-"motif")
    y1 <- y %>% dplyr::select(-"motif")
    x2 <- as.data.frame(t(x1))
    y2 <- as.data.frame(t(y1))
    df <- rbind(x2,y2)
    mat1 <- as.matrix(x2)
    mat2 <- as.matrix(y2)
    lgd_ = rep(NA, 11)
    lgd_[c(1,3,5,7,9,11)] = c(0.0,0.2,0.4,0.6,0.8,1.0)
    col <- colorRampPalette(c("black","#ffa10c"))
    heatmap(cor(mat1,mat2, method = "spearman"), 
                  col = col(11), symm = T, cexRow = 0.5,
                  cexCol = 0.5, xlab = z[1], ylab = z[2])
    legend(x =0.85, y = 0.9,legend=rev(lgd_),
           fill=rev(col(11)),
           title = "Correlation \n",
           horiz = F,
           bty = "n",
           cex = 0.75,
           title.adj = 2,
           y.intersp = 0.5)
   res <- heatmap(cor(mat1,mat2, method = "spearman"), 
            col = col(11), symm = T, cexRow = 0.5,
            cexCol = 0.5, xlab = z[1], ylab = z[2])
   return(colnames(mat1)[res$colInd])
}

x #list of bamfiles returned by bamfile()
y #list of character vectors with motifs to group fragments by
z #Name of groups
t #title of plot
fragment_length_motif_plot <- function(x, y, z, t){
    library(Biostrings)
    for(i in 1:length(x)){
        print(i)
        positive <- x[[i]][strand(x[[i]]) == "+"]
        positive_g1 <- positive[as.character(subseq(mcols(positive)$seq,1,3)) %in% y[[1]]]
        positive_g2 <- positive[as.character(subseq(mcols(positive)$seq,1,3)) %in% y[[2]]]
        negative <- x[[i]][strand(x[[i]]) == "-"]
        rev_strands <- DNAStringSet(stringi::stri_reverse(mcols(negative)$seq))
        negative_g1 <- negative[as.character(subseq(rev_strands,1,3)) %in% y[[1]]]
        negative_g2 <- negative[as.character(subseq(rev_strands,1,3)) %in% y[[2]]]
        g1_length_positive <- mcols(positive_g1)$isize
        g2_length_positive <- mcols(positive_g2)$isize
        g1_length_negative <- abs(mcols(negative_g1)$isize)
        g2_length_negative <- abs(mcols(negative_g2)$isize)
        g1_length <- c(g1_length_positive, g1_length_negative)
        g2_length <- c(g2_length_positive, g2_length_negative)
        g1_length <- g1_length[g1_length > 0]
        g1_length <- g1_length[!is.na(g1_length)]
        g2_length <- g2_length[g2_length > 0]
        g2_length <- g2_length[!is.na(g2_length)]
        df <- data.frame(len = c(g1_length,g2_length),
                         Sample = c(rep(z[1],length(g1_length)),rep(z[2],length(g2_length))))
        if(i == 1){
            df1 <- df
        }
        else{
            df1 <- rbind(df1,df)
        }
    }
    g1 <- df1 %>% filter(Sample == z[1])
    g2 <- df1 %>% filter(Sample == z[2])
    g2_median <- median(g2$len)
    g1_median <- median(g1$len)
    print("Analysis done, printing plot")
    gg <- ggplot(data = df1, aes(x = len, color = Sample))+
        geom_density(size = 1)+
        geom_vline(xintercept = g1_median, color="#ffa10c",
                   linetype="dashed", size = 1)+
        geom_vline(xintercept = g2_median, color="#6a00fc",
                   linetype="dashed", size = 1)+
        theme_bw(base_size = 15)+
        scale_color_manual(values = c("#6a00fc","#ffa10c"),
                           name = "Group")+
        labs(title = t, y = "Density", x = "Fragment length (bp)")+
        scale_x_continuous(breaks = seq(0,400,50), limits = c(0,400))+
        th
    return(gg)
}
x #Numeric vector of the motif proportions
entrop <- function(x){
    res <- -sum(x * log2(x))
    return(res)
}
x #Proportion of fragment end motifs for sample 1
y #Proportion of fragment end motifs for sample 2
shannon_entropy_plot <- function(x,y){
    library(dplyr)
    library(tidyr)
    library(ggpubr)
    library(ggplot2)
    x <- x %>% dplyr::select(-motif)
    y <- y %>% dplyr::select(-motif)
    e_1 <- c()
    e_2 <- c()
    for (i in 1:length(x)){
        e_1[i] <- entrop(x[,i])
        e_2[i] <- entrop(y[,i])
    }
    name_1 <- strsplit(colnames(x)[1], "_")[[1]][2]
    name_2 <- strsplit(colnames(y)[1], "_")[[1]][2]
    df <- data.frame(type = c(rep(name_1,length(x)),
                              rep(name_2,length(y))),
                     entropy = c(e_1, e_2),
                     pairs = rep(seq(1:length(x)),2))
    gg <- ggplot(df,aes(x = type, y = entropy, fill = type))+
        geom_boxplot(outlier.alpha = 0)+
        geom_point(alpha = 0.5)+
        labs(x = "",
             y = "Shannon entropy (H)")+
        geom_line(aes(group=pairs)) +
        scale_fill_manual(values = c("#6a00fc","#ffa10c"),
                          name = "Sample")+
        stat_compare_means(aes(x = type, y = entropy),
                           method = "t.test",
                           label = "p.format")+
        theme_bw()+
        th
    return(gg)
    
}
x #Proportion of fragment end motifs for sample 1
y #Proportion of fragment end motifs for sample 2
transposed_shannon_entropy_plot <- function(x,y){
    library(ggrepel)
    rownames(x) <- x$motif
    rownames(y) <- y$motif
    x <- x %>% dplyr::select(-motif)
    y <- y %>% dplyr::select(-motif)
    x1 <- as.data.frame(t(x))
    y1 <- as.data.frame(t(y))
    e_1 <- c()
    e_2 <- c()
    fc <- c()
    for (i in 1:length(x1)){
        e_1[i] <- entrop(x1[,i])
        e_2[i] <- entrop(y1[,i])
        fc[i] <- log2(entrop(x1[,i])/entrop(y1[,i]))
    }
    name_1 <- strsplit(rownames(x1)[1], "_")[[1]][2]
    name_2 <- strsplit(rownames(y1)[1], "_")[[1]][2]
    df <- data.frame(type = c(rep(name_1,length(x1)),
                              rep(name_2,length(y1))),
                     entropy = c(e_1, e_2),
                     fc = rep(fc,2),
                     pairs = rep(seq(1:length(x1)),2),
                     motif = rep(colnames(x1),2)) %>% 
        mutate(sig = ifelse(fc < -0.40 | fc > 0.4, TRUE,FALSE))
    gg <- ggplot(df,aes(x = type, y = entropy, fill = type))+
    geom_boxplot(outlier.alpha = 0)+
    geom_point(alpha = 0.5)+
    labs(x = "",
         y = "Shannon entropy (H)")+
    geom_line(aes(group=pairs,
                  color = sig)) +
    scale_fill_manual(values = c("#6a00fc","#ffa10c"),
                      name = "Sample")+
    scale_color_manual(values = c("TRUE"="black", "FALSE"="grey"))+
    stat_compare_means(aes(x = type, y = entropy),
                       method = "t.test",
                       label = "p.format")+
    geom_label_repel(
        aes(label=ifelse(sig == TRUE, as.character(motif),"")),
        segment.color="black",
        color="black",
        nudge_x = 0,
        label.size = NA,
        size = 2.5,
        fill = alpha(c("white"),0),
        parse = F,
        max.overlaps = 100)+
    theme_bw()+
    th
    return(gg)
}

x #list of motifs to separate motifs into
y #list of data.frames to to get information from
z #character vector with names of samples
stack_bar <- function(x,y,z,m){
    df <- data.frame(samples = z)
    for(i in 1:length(x)){
        seqs <- x[[i]]
        ress <- c()
        for (j in 1:length(y)){
            ddf <- y[[j]] %>% filter(motif %in% seqs) %>%
                select(-motif)
            frac <- colSums(ddf)
            res <- mean(frac)
            ress[j] <- res
            
        }
        df[ , ncol(df) + 1] <- ress
    }
    df <- df %>% pivot_longer(contains("V")) %>% 
        mutate(name = rep(c(paste(m[1]), paste(m[2]), paste(m[3])), length(z))) %>% 
        mutate(samples = as.factor(samples))
    gg <- ggplot(df, 
                 aes(fill = name, y = value, x = as.factor(samples),
                     label = round(value,3)))+
        geom_bar(position = "stack", stat = "identity")+
        labs(x = "Sample",
             y = "Motif fraction",
             caption = paste0(m[1], ": ", 
                              paste(x[[1]], sep = "", collapse = ", "),
                              "\n",
                             m[2], ": ", 
                             paste(x[[2]], sep = "", collapse = ", "),
                             "\n",
                             m[3], ": ", 
                             paste(x[[3]], sep = "", collapse = ", ")))+
        theme_bw()+
        scale_fill_manual(name = "Group",
                          values = c("#6a00fc","grey", "#ffa10c"))+
        geom_text(size = 3, position = position_stack(vjust = 0.5))+ 
        theme(plot.caption.position = "plot",
              plot.caption = element_text(hjust = 0))+   
        th
    return(gg)
}

x #list of BAMfiles returned by bamfile() for group 1
y #list of BAMfiles returned by bamfile() for group 2
z #vector of names for group 1 and 2
p #Title of plot
fragment_length <- function(x,y,z,p){
    library(ggplot2)
    library(tidyr)
    len1 <- c()
    len2 <- c()
    for (i in 1:length(x)){
        len <- mcols(x[[i]])$isize
        len <- len[len >0]
        len <- len[!is.na(len)]
        len1 <- c(len1,len)
    }
    for (i in 1:length(y)){
        len <- mcols(y[[i]])$isize
        len <- len[len >0]
        len <- len[!is.na(len)]
        len2 <- c(len2,len)
    }
    df <- data.frame(len = c(len1,len2),
                     Sample = factor(c(rep(z[1],length(len1)),
                                       rep(z[2],length(len2))),
                                     levels = c(z[1], z[2])))
    mu <- ddply(df, "Sample", summarise, grp.median=median(len))
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

x #BAM file returned by bamfile()
y #cutoff for fragment length
proportion_sub <- function(x,y){
    library(Rsamtools)
    len <- mcols(x)$isize
    len <- len[len >0]
    len <- len[!is.na(len)]
    bam_lengths_sub150 <- len[len<y]
    res <- length(bam_lengths_sub150)/length(len)
    return(res)
}
proportion_sub(NAC.1_cfChIP, 150)

x #list of bam files returned by bamfile() for group 1
y #cutoff for fragment length
z #name of samples
proportion_sub_df <- function(x, y, z){
    res <- c()
    for (i in 1:length(x)){
        res[i] <- proportion_sub(x[[i]], y)
    }
    df <- data.frame(sample = z,
                     sub_fraction = res)
    return(df)
}

x #list of proportion sub data.frame returned by proportion_sub_df() for samples
z #Name of samples
p #If samples are paired p = "Paired (default), otherwise p = "Unpaired"
proportion_sub_boxplot <- function(x, z, p = "Paired"){
    library(ggpubr)
    if(p == "Unpaired"){
        if (length(z) < 3){
            df <- data.frame(sub_fraction = c(x[[1]]$sub_fraction,
                                              x[[2]]$sub_fraction),
                             sample = factor(c(rep(z[1], nrow(x[[1]])),
                                        rep(z[2], nrow(x[[2]]))),
                                        levels = c(z[1], z[2])))
            gg <- df %>% ggplot(aes(x=sample, y=sub_fraction, fill=sample)) + 
                geom_boxplot(outlier.alpha = 0)+
                geom_jitter(alpha = 0.5, width = 0.1)+
                labs(x = "",
                     y = "Fraction under 150 bp")+
                scale_fill_manual("Sample",
                                  values = c("#ffa10c","#6a00fc"))+
                stat_compare_means(method = "t.test",
                                   comparisons = list(c(z[1], z[2])))+
                theme_bw()+
                th
        }
        else{
            df <- data.frame(sub_fraction = c(x[[1]]$sub_fraction,
                                              x[[2]]$sub_fraction,
                                              x[[3]]$sub_fraction),
                             sample = factor(c(rep(z[1], nrow(x[[1]])),
                                        rep(z[2], nrow(x[[2]])),
                                        rep(z[3], nrow(x[[3]]))),
                                        levels = c(z[1], z[2], z[3])))
            gg <- df %>% ggplot(aes(x=sample, y=sub_fraction, fill=sample)) + 
                geom_boxplot(outlier.alpha = 0)+
                geom_jitter(alpha = 0.5, width = 0.1)+
                labs(x = "",
                     y = "Fraction under 150 bp")+
                scale_fill_manual("Sample",
                                  values = c("#ffa10c","firebrick","#6a00fc"))+
                stat_compare_means(method = "t.test",
                                   comparisons = list(c(z[1],z[2]),
                                                      c(z[1],z[3]),
                                                      c(z[2],z[3])))+
                theme_bw()+
                th
        }
    }
    else{
        df <- data.frame(sub_fraction = c(x[[1]]$sub_fraction,
                                          x[[2]]$sub_fraction),
                         sample = factor(c(rep(z[1], nrow(x[[1]])),
                                    rep(z[2], nrow(x[[2]]))), 
                                    levels = c(z[1], z[2])),
                         pairs = rep(seq(1,nrow(x[[1]])),2))
        gg <- df %>% ggplot(aes(x=sample, y=sub_fraction, fill=sample)) + 
            geom_boxplot(outlier.alpha = 0)+
            geom_line(aes(group=pairs))+
            geom_point(alpha = 0.5)+
            labs(x = "",
                 y = "Fraction under 150 bp")+
            scale_fill_manual("Sample",
                              values = c("#ffa10c","#6a00fc"))+
            stat_compare_means(method = "t.test",
                               paired = TRUE,
                               comparisons = list(c(z[1], z[2])))+
            theme_bw()+
            th
    }
    return(gg)
}
x #Table consisting of BAM file names (name used when bamfile() is used),
  #genes mutated and the position of the mutation
y #Named list of samples used
fragment_length_wt_ctdna_df <- function(x,y){
    library(Rsamtools)
    library(tidyr)
    library(GenomicAlignments)
    or_wd <- getwd()
    `%ni%` <- Negate(`%in%`)
    df <- read.table(x,header = T)
    df <- df[!is.na(df$genes),]
    df6 <- data.frame(len = c(),
                      Sample = c())
    n_samples <- length(df$file)
    for (i in 1:n_samples){
        genes <- strsplit(df$genes[i],split = ",")[[1]]
        n_genes <- length(genes)
        positions <- strsplit(df$positions[i],split = ",")[[1]]
        types <- strsplit(df$Type[i],split = ",")[[1]]
        string <- strsplit(positions,split = ":")[1:n_genes]
        reads <- y[df$name[i]][[1]]
        chr <- unlist(lapply(string, FUN = function(x){return(x[1])}))
        pos <- as.numeric(unlist(lapply(string, FUN = function(x){return(x[2])})))
        for(j in 1:n_genes){
            mut_pos <- pos[j]
            mut_chr <- chr[j]
            g <- GRanges(mut_chr,IRanges(start = mut_pos, end = mut_pos),
                         strand = "*")
            sub <- subsetByOverlaps(reads,g)
            type <- types[j]
            len <- mcols(sub)$isize
            cig <- mcols(sub)$cigar
            cig <- cig[len >0]
            len <- len[len >0]
            cig <- cig[!is.na(len)]
            len <- len[!is.na(len)]
            if(grepl("indel",type)){
                indel <- strsplit(type,split = "-")[[1]][2]
                ct_len <- len[grepl(indel,cig)]
                wt_len <- len[!grepl(indel,cig)]
                df1 <- data.frame(len = c(ct_len,wt_len),
                                 Sample = c(rep("ctDNA",length(ct_len)),
                                            rep("WT", length(wt_len))))
                df6 <- rbind(df6,df1)
            }
            else{
                param = ScanBamParam(which = g, what = c("mapq", "isize","pos", "cigar"),
                                     tag = "MD")
                setwd("D:/Lung cancer input/PosDeduped")
                bam <- scanBam(df$file[i],
                               param = param)
                setwd(or_wd)
                bam_pos <- as.numeric(unlist(bam[[1]][1]))
                bam_size <- as.numeric(unlist(bam[[1]][4]))
                bam_cigar <- unlist(bam[[1]][3])
                bam_MD <- unlist(bam[[1]][5])
                df1 <- data.frame(position = bam_pos,
                                  size = bam_size,
                                  cigar = bam_cigar,
                                  MD = bam_MD)
                df1 <- df1[df1$size >0,]
                df1 <- df1[!is.na(df1$position),]
                base <- strsplit(type,split = "-")[[1]][2]
                ct_len <- df1$size[grepl(base,df1$MD)]
                wt_len <- df1$size[!grepl(base,df1$MD)]
                df4 <- data.frame(len = c(ct_len,wt_len),
                                  Sample = c(rep("ctDNA",length(ct_len)),
                                             rep("WT", length(wt_len))))
                df6 <- rbind(df6,df4)
                print(i)
            }
        }    
    }
    return(df6)
}

x #data.frame returned by fragment_length_wt_ctdna_df()
fragment_length_wt_ctdna_plot <- function(x){
    library(ggplot2)
    library(tidyr)
    library(plyr)
    colnames(x) <- c("len", "Sample")
    x$Sample <- factor(x$Sample, levels = c("ctDNA", "WT"))
    mu <- ddply(x, "Sample", summarise, grp.median=median(len))
    gg <- ggplot(data = x, aes(x = len, color = Sample))+
        geom_density(size = 1)+
        geom_vline(data=mu, aes(xintercept=grp.median, color=Sample),
                   linetype="dashed", size = 1)+
        theme_bw(base_size = 15)+
        scale_color_manual(values = c("#6a00fc","#ffa10c"))+
        labs(title = "ctDNA and WT", y = "Density", x = "Fragment length (bp)")+
        scale_x_continuous(breaks = seq(0,400,50), limits = c(0,400))+
        th
    return(gg)
}

MAF_in_bins <- function(x){
    intervals <- seq(50,400,10)
    colnames(x) <- c("len", "Sample")
    x$Sample <- factor(x$Sample, levels = c("ctDNA", "WT"))
    x_mark <- c()
    df_wt <- x %>% filter(Sample == "WT")
    df_ct <- x %>% filter(Sample == "ctDNA")
    or_MAF <- nrow(df_ct)/(nrow(df_ct)+nrow(df_wt))*100
    MAF <- c()
    for (i in 1:(length(intervals)-1)){
        lower <- intervals[i]
        upper <- intervals[i+1]
        x_mark[i] <- lower + 5
        df <- x %>% filter(lower < len & len <= upper)
        df_wt <- df %>% filter(Sample == "WT")
        df_ct <- df %>% filter(Sample == "ctDNA")
        MAF[i] <- nrow(df_ct)/(nrow(df_ct)+nrow(df_wt))*100
    }
    rects <- data.frame(xstart = c(50,150,230,340), xend = c(150,230,340,400), 
                        col = factor(c("Mutant","WT","Mutant","WT"),
                                        levels = c("Mutant", "WT")))
    df1 <- data.frame(len = x_mark,
                      MAF = MAF) %>% 
    mutate(enrichment = MAF/or_MAF)
    gg <- ggplot()+
        geom_rect(data = rects, aes(xmin = xstart, xmax = xend, 
                                    ymin = -Inf, ymax = Inf, fill = col),
                  alpha = 0.4)+
        scale_fill_manual("Enrichment",
                          values = c("#6a00fc","#ffa10c"))+ 
        geom_point(data = df1, aes(x = len, y = enrichment),
                   color = "black",
                   size = 2)+
        labs(x = "Fragment (bp)",
             y = "MAF enrichment")+
        geom_hline(yintercept = 1,
                   color = "black",
                   size = 1)+
        geom_smooth(data = df1, aes(x = len, y = enrichment),
                    method = "loess",
                    span = 0.5, se = F,
                    color = "black",
                    size = 1.2)+
        geom_vline(xintercept = c(50,150,230,340,400),
                   color = "black",
                   linetype = "dashed",
                   size = 0.5)+
        annotate(geom = "text",
                 x = c(58,160,240,350,410),
                 y = c(2.4,2.4,2.4,2.4,2.4),
                 label = c("50", "150", "230", "340", "400"),
                 color = "black",fontface = "bold",
                 size = 4)+
        scale_y_continuous(limits = c(0.4,2.4))+
        scale_x_continuous(limits = c(0,450),breaks = seq(0,450,50))+
        theme_bw()+
        th
    return(gg)
}

x #Table consisting of BAM file names (name used when bamfile() is used),
#genes mutated and the position of the mutation
y #Named list of the patient sample (length = 1)
fragment_length_wt_ctdna_patient <- function(x,y){
    library(Rsamtools)
    library(tidyr)
    library(GenomicAlignments)
    or_wd <- getwd()
    `%ni%` <- Negate(`%in%`)
    df <- read.table(x,header = T)
    df <- df[!is.na(df$genes),]
    df6 <- data.frame(len = c(),
                      Sample = c())
    df <- df %>% filter(name %in% names(y))
    if(isEmpty(df)){
        return("No mutations")
    }
    genes <- strsplit(df$genes,split = ",")[[1]]
    n_genes <- length(genes)
    positions <- strsplit(df$positions,split = ",")[[1]]
    types <- strsplit(df$Type,split = ",")[[1]]
    string <- strsplit(positions,split = ":")[1:n_genes]
    reads <- y[df$name][[1]]
    chr <- unlist(lapply(string, FUN = function(x){return(x[1])}))
    pos <- as.numeric(unlist(lapply(string, FUN = function(x){return(x[2])})))
        for(j in 1:n_genes){
            mut_pos <- pos[j]
            mut_chr <- chr[j]
            g <- GRanges(mut_chr,IRanges(start = mut_pos, end = mut_pos),
                         strand = "*")
            sub <- subsetByOverlaps(reads,g)
            type <- types[j]
            len <- mcols(sub)$isize
            cig <- mcols(sub)$cigar
            cig <- cig[len >0]
            len <- len[len >0]
            cig <- cig[!is.na(len)]
            len <- len[!is.na(len)]
            if(grepl("indel",type)){
                indel <- strsplit(type,split = "-")[[1]][2]
                ct_len <- len[grepl(indel,cig)]
                wt_len <- len[!grepl(indel,cig)]
                df1 <- data.frame(len = c(ct_len,wt_len),
                                  Sample = c(rep("ctDNA",length(ct_len)),
                                             rep("WT", length(wt_len))))
                df6 <- rbind(df6,df1)
            }
            else{
                param = ScanBamParam(which = g, what = c("mapq", "isize","pos", "cigar"),
                                     tag = "MD")
                setwd("D:/Lung cancer input/PosDeduped")
                bam <- scanBam(df$file,
                               param = param)
                setwd(or_wd)
                bam_pos <- as.numeric(unlist(bam[[1]][1]))
                bam_size <- as.numeric(unlist(bam[[1]][4]))
                bam_cigar <- unlist(bam[[1]][3])
                bam_MD <- unlist(bam[[1]][5])
                df1 <- data.frame(position = bam_pos,
                                  size = bam_size,
                                  cigar = bam_cigar,
                                  MD = bam_MD)
                df1 <- df1[df1$size >0,]
                df1 <- df1[!is.na(df1$position),]
                base <- strsplit(type,split = "-")[[1]][2]
                ct_len <- df1$size[grepl(base,df1$MD)]
                wt_len <- df1$size[!grepl(base,df1$MD)]
                df4 <- data.frame(len = c(ct_len,wt_len),
                                  Sample = c(rep("ctDNA",length(ct_len)),
                                             rep("WT", length(wt_len))))
                df6 <- rbind(df6,df4)
            }
        }    
    return(df6)
}

x #List of patient data.frames with fragment lengths of WT and ctDNA
fragment_length_wt_ctdna_boxplot <- function(x){
    x <- lapply(x,FUN = function(x){if(isa(x,"data.frame")){return(x)}
        else{return(NA)}})
    x <- x[!is.na(x)]
    ct_frac <- c()
    wt_frac <- c()
    for (i in 1:length(x)){
        df_ct <- x[[i]] %>% filter(Sample == "ctDNA")
        df_wt <- x[[i]] %>% filter(Sample == "WT")
        df_ct_sub <- df_ct %>% filter(len < 150)
        df_wt_sub <- df_wt %>% filter(len < 150)
        ct_frac[i] <- nrow(df_ct_sub)/nrow(df_ct)
        wt_frac[i] <- nrow(df_wt_sub)/nrow(df_wt)
    }
    df <- data.frame(sub_fraction = c(wt_frac,ct_frac),
                     sample = factor(c(rep("WT",length(x)),
                              rep("ctDNA", length(x))),
                                     levels = c("WT", "ctDNA")),
                     pairs = rep(1:length(x),2))
    gg <- df %>% ggplot(aes(x=sample, y=sub_fraction, fill=sample)) + 
        geom_boxplot(outlier.alpha = 0)+
        geom_line(aes(group=pairs))+
        geom_point(alpha = 0.5)+
        labs(x = "",
             y = "Fraction under 150 bp")+
        scale_fill_manual("Fragments",
                          values = c("#ffa10c","#6a00fc"))+
        stat_compare_means(method = "t.test",
                           paired = TRUE,
                           comparisons = list(c("WT", "ctDNA")))+
        theme_bw()+
        th
    return(gg)
}
x #enrichment data.frame 
y #Granges object returned by gr() for the 197 AVENIO target genes
z #name of sample to isolate from enrichment data.frame
i #name of input BAM file
fragment_length_active_inactive_df <- function(x,y,z,i){
    library(Rsamtools)
    df_sample <- x[c("genes",z)]
    df_sample <- df_sample[order(df_sample[,2]),]
    e <- df_sample
    top15 <- e$genes[(length(e$genes)-14):length(e$genes)]
    bottom15 <- e$genes[1:15]
    which <- y[elementMetadata(y)[,1] %in% c(top15)]
    index1 <- match(top15,elementMetadata(which)[,1])
    which <- which[index1]
    what <- c("mapq", "isize")
    param <- ScanBamParam(which = which, what = what)
    bam <- scanBam(i, param = param)
    length_top15 <- unname(unlist(lapply(bam,function(x) x$isize)))
    length_top15 <- length_top15[length_top15 >0]
    length_top15 <- length_top15[!is.na(length_top15)]
    which <- grs[elementMetadata(y)[,1] %in% c(bottom15)]
    index1 <- match(bottom15,elementMetadata(which)[,1])
    which <- which[index1]
    param <- ScanBamParam(which = which, what = what)
    bam <- scanBam(i, param = param)
    length_bottom15 <- unname(unlist(lapply(bam,function(x) x$isize)))
    length_bottom15 <- length_bottom15[length_bottom15 >0]
    length_bottom15 <- length_bottom15[!is.na(length_bottom15)]
    df <- data.frame(len = c(length_top15,length_bottom15),
                     Sample = c(rep("Active",length(length_top15)),
                                rep("Inactive",length(length_bottom15))))
    return(df)
}
fragment_length_active_inactive_df(enrichment_df,
                                   grs,
                                   "NAC.1",
                                   "D:/Lung cancer input/PosDeduped/PosDeduped-A-1279-input.bam")
x #Data.frame with two columns. One with fragment lengths and the other determining whether the fragment originates from an active or inactive gene
p #Title of plot
fragment_length_active_inactive_plot <- function(x,p){
    library(ggplot2)
    library(tidyr)
    library(plyr)
    df <- x
    mu <- ddply(df, "Sample", summarise, grp.median=median(len))
    gg <- ggplot(data = df, aes(x = len, color = Sample))+
        geom_density(size = 1)+
        geom_vline(data=mu, aes(xintercept=grp.median, color=Sample),
                   linetype="dashed", size = 1)+
        theme_bw(base_size = 15)+
        scale_color_manual(values = c("#6a00fc","#ffa10c"),
                           name = "Gene activity")+
        labs(title = p, y = "Density", x = "Fragment length (bp)")+
        scale_x_continuous(breaks = seq(0,400,50), limits = c(0,400))+
        th
    return(gg)
}

x #list of fragment length data.frames returned by fragment_length_active_inactive_df()
y #fragment length cutoff
z #name of samples
a #Definition of whether "Active" or "Inactive" genes are analyzed
fragment_length_active_inactive_sub_df <- function(x,y,a,z){
    res <- c()
    for (i in 1:length(x)){
        ddf <- x[[i]]
        ddf <- ddf %>% filter(Sample == a)
        bam_lengths_sub150 <- ddf$len[ddf$len<y]
        ress <- length(bam_lengths_sub150)/length(ddf$len)
        res[i] <- ress
    }
    df <- data.frame(sample = z,
                     sub_fraction = res)
    return(df)
}
x #enrichment data.frame 
y #Granges object returned by gr() for the 197 AVENIO target genes
z #name of sample to isolate from enrichment data.frame
j #name of input BAM file
p #cutoff for fragment length
sub_150_in_cfChIP_quartiles <- function(x,y,z,j,p){
    library(Rsamtools)
    df_sample <- x[c("genes",z)]
    df_sample <- df_sample[order(df_sample[,2]),]
    e <- df_sample
    number_high <- 19
    number_low <- 1
    res <- c()
    for (i in 1:10){
        if(i < 5){
            genes <- e$genes[number_low:number_high]
            number_low <- number_low + 19
            number_high <- number_high + 19
            which <- y[elementMetadata(y)[,1] %in% genes]
            index1 <- match(genes,elementMetadata(which)[,1])
            which <- which[index1]
            what <- c("mapq", "isize")
            param <- ScanBamParam(which = which, what = what)
            bam <- scanBam(j, param = param)
            len <- unname(unlist(lapply(bam,function(x) x$isize)))
            len <- len[len >0]
            len <- len[!is.na(len)]
            bam_lengths_sub150 <- len[len<p]
            ress <- length(bam_lengths_sub150)/length(len)
            res[i] <- ress
            if(i == 4){
                number_high <- number_high + 1
            }
        }
        else{
            genes <- e$genes[number_low:number_high]
            number_low <- number_low + 20
            number_high <- number_high + 20
            which <- y[elementMetadata(y)[,1] %in% genes]
            index1 <- match(genes,elementMetadata(which)[,1])
            which <- which[index1]
            what <- c("mapq", "isize")
            param <- ScanBamParam(which = which, what = what)
            bam <- scanBam(j, param = param)
            len <- unname(unlist(lapply(bam,function(x) x$isize)))
            len <- len[len >0]
            len <- len[!is.na(len)]
            bam_lengths_sub150 <- len[len<p]
            ress <- length(bam_lengths_sub150)/length(len)
            res[i] <- ress
        }
    }
    df <- data.frame(quartile = c("Q1", "Q2", "Q3", "Q4",
                                  "Q5", "Q6", "Q7", "Q8",
                                  "Q9", "Q10"),
                     sub_fraction = res,
                     sample = rep(z,10))
    return(df)
}

x #data.frame with quartiles, sub150 fraction of the quartiles and the sample
sub150_fraction_quartiles_plot <- function(x){
    ddf <- data.frame(quartile = c(),
                      sub_fraction = c(),
                      sample = c(),
                      norm_fraction = c())
    for(i in 1:length(unique(x$sample))){
        sam <- unique(x$sample)[i]
        df <- x %>% filter(sample == sam)
        base_fraction <- df$sub_fraction[1]
        df <- df %>% mutate(norm_fraction = sub_fraction/base_fraction)
        ddf <- rbind(ddf,df)
    }
    ddf <- ddf %>% 
        mutate(quartile = factor(quartile,
                                 level = c("Q1", "Q2", "Q3", "Q4",
                                           "Q5", "Q6", "Q7", "Q8",
                                           "Q9", "Q10")))
    gg <- ddf %>% ggplot(aes(x = quartile, y = norm_fraction))+
        geom_boxplot(outlier.alpha = 0,
                     fill = "#6a00fc")+
        geom_jitter(alpha = 0.5, width = 0.1)+
        theme_bw()+
        stat_compare_means(method = "anova", label.y = 0.9,
                           label.x = 8)+ 
        stat_compare_means(aes(label = after_stat(p.signif)),
                           method = "t.test",
                           paired = TRUE,
                           ref.group = "Q1",
                           size = 4.5)+
        labs(x = "Quartiles",
             y = "Normalized fraction under 150 bp")+
        scale_y_continuous(limits = c(0.9,1.5))+
        th
    
    return(gg)
    
}
sub150_fraction_quartiles_plot(collected_sub150_fraction_quartiles)
?geom_boxplot()
comparisons = list(c("Q1", "Q2"),
                   c("Q1", "Q3"),
                   c("Q1", "Q4"),
                   c("Q1", "Q5"),
                   c("Q1", "Q6"),
                   c("Q1", "Q7"),
                   c("Q1", "Q8"),
                   c("Q1", "Q9"),
                   c("Q1", "Q10"))
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
#x BAM file returned by bamfile()
#m #number of bases to be evaluated in each fragment end
mer_count <- function(x,m){
    library(GenomicRanges)
    library(chromstaR)
    library(Biostrings)
    positve <- x[strand(x) == "+"]
    negative <- x[strand(x) == "-"]
    rev_strands <- reverseComplement(mcols(negative)$seq)
    positve_counts <- nucleotideFrequencyAt(
        x = mcols(positve)$seq, at = c(1,2,3), as.array = F)
    negative_counts <- nucleotideFrequencyAt(
        x = rev_strands,
        at = c(1,2,3), as.array = F)
    return(positve_counts+negative_counts)
    
}

#x table returned by mer_count()
prop_df <- function(x){
    prop <- proportions(x)
    df <- data.frame(res = names(prop),
                     Freq = unname(prop))
    return(df)
}

#x table returned by mer_count()
count_df <- function(x){
    df <- data.frame(res = names(x),
                     Freq = unname(x))
    return(df) 
}
#x #prop_mer data.frame for group 1
#y #prop_mer data.frame for group 2

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

#x #Name of group 1
#y #Name of group 2
#z #data.frame returned by dif_motif_prop() eller diff_active_inactive_end_motif()
volcano_motif <- function(x,y,z){
    library(ggplot2)
    library(ggrepel)
    p <- -log10(z$q.value)
    z$Log10p <- p
    z <- z %>% dplyr::filter(!is.na(q.value))
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
        labs(title = "")+
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

#x #Proportion of fragment end motifs for sample 1
#y #Proportion of fragment end motifs for sample 2
#p #Name of sample 1
#q #Name of sample 2
#s #seed for set.seed()
#n #Number of neighbors
#t #Type of plot. For healthy and cancer type "both" (default), for cancer = "cancer", for healthy = "healthy"
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
            labs(title = "",
                 x = "UMAP-1", y = "UMAP-2")+
            scale_color_manual(name = "Sample",
                               labels = c("Cancer input", "Cancer cfChIP",
                                          "Healthy input", "Healthy cfChIP"),
                               values = c("#ffa10c","#6a00fc","#ffa10c","#6a00fc"))+
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
            labs(title = "",
                 x = "UMAP-1", y = "UMAP-2")+
            scale_color_manual(name = "Sample",
                               labels = c("Cancer input", "Cancer cfChIP"),
                               values = c("#ffa10c","#6a00fc"))+
            theme_bw()+
            th
    }
    if(t == "healthy"){
        gg <- ggplot(data = umap_df, aes(x = V1, y = V2, color = V3))+
            geom_point(size = 5)+
            labs(title = "",
                 x = "UMAP-1", y = "UMAP-2")+
            scale_color_manual(name = "Sample",
                               labels = c("Healthy input", "Healthy cfChIP"),
                               values = c("#ffa10c","#6a00fc"))+
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
            labs(title ="",
                 x = "UMAP-1", y = "UMAP-2")+
            scale_color_manual(name = "Sample",
                               labels = c("Cancer High expressed", "Cancer Low expressed",
                                          "Healthy High expressed", "Healthy Low expressed"),
                               values = c("#6a00fc","#ffa10c","#6a00fc","#ffa10c"))+
            scale_shape_manual(name = "Sample",
                               labels = c("Cancer High expressed", "Cancer Low expressed",
                                          "Healthy High expressed","Healthy Low expressed"),
                               values = c(19,19,17,17))+
            theme_bw()+
            th
    }
    
    return(gg)
}

#x #enrichment data.frame 
#y #Granges object returned by gr() for the 197 AVENIO target genes
#z #name of sample to isolate from enrichment data.frame
#i #name of input BAM file
#m #number of bases to be evaluated in each fragment end
end_motif_active_inactive <- function(x,y,z,i,m){
    library(Rsamtools)
    df_sample <- x[c("genes",z)]
    df_sample <- df_sample[order(df_sample[,2]),]
    e <- df_sample
    top15 <- e$genes[(length(e$genes)-14):length(e$genes)]
    bottom15 <- e$genes[1:15]
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

#x #prop_mer data.frame with paired active and inactive gene fragment end motifs
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

#x #data.frame returned by diff_active_inactive_end_motif() eller dif_motif_prop()
#q #value cutoff
#f #Can be "negative" eller "positive" for om ente FC < 1 eller FC > 1 endemotiver skal undersøge
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

#x #list of character vectors with significant sequences
#y #character vector of names of samples
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

#x #bam.file returned by bamfile()
#y #genomic ranges of AVENIO genes returned by gr()
#z #character vector of motifs to get fraction of
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
                                        

#x #name on .txt file containing each gene, the number of times each gene is sequenced and the distance each target is from TSS
#y #cutoff of distance from TSS which is too long
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
#x #data.frame returned by motif_fraction()
#y #data.frame returned by e.score()
#z #Name of sample
#b #optional: if not NULL then use character vector returned by badgene()
#r #optional: if not NULL then use character vector returned by multiple_region()
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
#x #correlation matrix with samples, rhos, p and alpha values
#z #Title of plot
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
#x #BAM file returned by bamfile()
#y #fragment length cut-off, default = 150
#z #TRUE to get fragments shorter than or equal to cutoff, whereas FALSE gives fragments longer than cutoff. Default is TRUE
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


#x #BAM file returned by bamfile()
#y #fragment length cut-off, default = 150
#m #number of bases to be evaluated in each fragment end
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

#x #list of character vectors with significant sequences
#y #character vector of names of samples
#tit #title of plot
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
            axis.ticks = element_blank())+
        scale_x_continuous(expand = expansion(mult = .2))
    return(gg)
}
#e #enrichment dataframe
#n #Name of sample
#x #input BAM file returned by bamfile()
#g #Granges object returned by gr()
#y #fragment length cut-off, default = 150
#m #number of bases to be evaluated in each fragment end. Default = 3
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

#x #list of data.frames returned by e.score()
#y #vector of names for each data.frame.
#z #vector of names for groups
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

#x #data.frame returned by e.score()
#y #name of sample
transpose <- function(x,y){
    x1 <- as.data.frame(t(x))
    colnames(x1) <- x1[1,1:length(x1)]
    x1 <- x1[-1,]
    rownames(x1) <- y
    return(x1)
}

#x #list of data.frames returned by e.score()
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

#x #Average Enrichment data.frame return from average_enrichmet() for sample 1
#y #Average Enrichment data.frame return from average_enrichmet() for sample 2
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

#x #data.frame with genes of interest returned by conf()
#g #Granges object returned by gr()
#i #input BAM file to be analyzed, returned by bamfile()
diff_gene_reads <- function(x,g,i){
    gene_g <- g[mcols(g)$SYMBOL %in% x$genes]
    sub <- subsetByOverlaps(i, gene_g,ignore.strand=TRUE)
    return(sub)
}

#x #input BAM file to be analyzed, returned by bamfile()
#y #cutoff for short fragments
sub150_fraction <- function(x,y){
    mcols(x)$abs_size <- abs(mcols(x)$isize)
    res <- mcols(x)$abs_size<y
    res <- res[!is.na(res)]
    return(mean(res))
}

#x #BAM file to be analyzed, returned by bamfile()
di_nucleosome_fraction <- function(x){
    mcols(x)$abs_size <- abs(mcols(x)$isize)
    res <- mcols(x)$abs_size<340 & mcols(x)$abs_size>230
    res <- res[!is.na(res)]
    return(mean(res))
}

#x named list of BAM files to be analyzed, returned by bamfile()
#y named list of BAM files to be analyzed, returned by bamfile()
#q named list of BAM files to be analyzed, returned by bamfile(), if not needed type q = NULL
#z Names of types of samples in x and y 
di_nucleosome_fraction_df <- function(x,y,q=NULL,z){
    df <- data.frame(name = c(),
                     proportions = c(),
                     sample = c())
    for(i in 1:length(x)){
        nam <- unlist(lapply(strsplit(names(x)[i],"_"), "[[", 1))
        res <- di_nucleosome_fraction(x[[i]])
        df1 <- data.frame(name = nam,
                          proportions = res,
                         sample = z[1])
        df <- rbind(df,df1)
    }
    for(i in 1:length(y)){
        nam <- unlist(lapply(strsplit(names(y)[i],"_"), "[[", 1))
        res <- di_nucleosome_fraction(y[[i]])
        df1 <- data.frame(name = nam,
                          proportions = res,
                         sample = z[2])
        df <- rbind(df,df1)
    }
    if(!is.null(q)){
        for(i in 1:length(q)){
            nam <- unlist(lapply(strsplit(names(q)[i],"_"), "[[", 1))
            res <- di_nucleosome_fraction(q[[i]])
            df1 <- data.frame(name = nam,
                              proportions = res,
                              sample = z[3])
            df <- rbind(df,df1)
        }
    }
    return(df)
}
#x #list of fragment length data.frames returned by fragment_length_active_inactive_df()
#aka read.table("Fragment lengths active inactive.txt", header = T)
di_nucleosome_high_low_fraction_df<- function(x){
    df <- data.frame(name = c(),
                     proportions = c(),
                     sample = c())
    for(i in 1:length(unique(x$group))){
        sample <- unique(x$group)[i]
        df_sample <- x %>% filter(group == sample)
        df_high <- df_sample %>% filter(Sample == "Active")
        df_low <- df_sample %>% filter(Sample == "Inactive")
        sub_high <- df_high %>% filter(len<340 & len>230)
        sub_low <- df_low %>% filter(len<340 & len>230)
        df1 <- data.frame(name = rep(sample,2),
                          proportions = c(nrow(sub_high)/nrow(df_high),
                                          nrow(sub_low)/nrow(df_low)),
                          sample = c("Cancer High expression","Cancer Low expression"))
        df <- rbind(df,df1)
    }
    df <- df %>% arrange(sample)
    return(df)
}

di_nucleosome_fraction_boxplot <- function(x,y,z=F){
    library(ggplot2)
    library(ggpubr)
    if(length(y)==2){
        if(z){
            df <- x %>% mutate(sample = factor(sample,levels = c(y[1],y[2]))) %>% 
                mutate(pairs = rep(1:length(unique(x$name)),2))
            gg <- ggplot(df, aes(x = sample, y = proportions,
                                 fill = sample))+
                geom_boxplot(outlier.alpha = 0)+
                geom_point(alpha = 0.5)+
                geom_line(aes(group=pairs))+
                labs(title = "",
                     x = "",
                     y = "Fraction 230-340 bp")+
                scale_fill_manual("Sample",
                                  values = c("#ffa10c","#6a00fc"))+
                stat_compare_means(method = "t.test",
                                   paired = T,
                                   comparisons = list(c(y[1], y[2])))+
                theme_bw()+
                th+
                theme(legend.position = "none")
        }
        else{
            df <- x %>% mutate(sample = factor(sample,levels = c(y[1],y[2])))
            gg <- ggplot(df, aes(x = sample, y = proportions,
                                 fill = sample))+
                geom_boxplot(outlier.alpha = 0)+
                geom_jitter(alpha = 0.5, width = 0.1)+
                labs(title = "",
                     x = "",
                     y = "Fraction 230-340 bp")+
                scale_fill_manual("Sample",
                                  values = c("#ffa10c","#6a00fc"))+
                stat_compare_means(method = "t.test",
                                   paired = FALSE,
                                   comparisons = list(c(y[1], y[2])))+
                theme_bw()+
                th+
                theme(legend.position = "none")
        }
        
    }
    if(length(y)==3){
        df <- x %>% mutate(sample = factor(sample,levels = c(y[1],y[2],y[3])))
        gg <- df %>% ggplot(aes(x=sample, y=proportions, fill=sample)) + 
            geom_boxplot(outlier.alpha = 0)+
            geom_jitter(alpha = 0.5, width = 0.1)+
            labs(title = "",
                 x = "",
                 y = "Fraction 230-340 bp")+
            scale_fill_manual("Sample",
                              values = c("#ffa10c","firebrick","#6a00fc"))+
            stat_compare_means(method = "t.test",
                               comparisons = list(c(y[1],y[2]),
                                                  c(y[1],y[3]),
                                                  c(y[2],y[3])))+
            theme_bw()+
            th+
            theme(legend.position = "none")
    }
    
    return(gg)
}

#x #data.frame with genes of interest returned by conf()
#g #Granges object returned by gr()
#l #list of input BAM file to be analyzed, returned by bamfile()
#y #cutoff for short fragments
#m #number of bases to be evaluated in each fragment end
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

#x #data.frame returned by diff_gene_epigenetic_analysis()[[2]] for sample 1
#y #data.frame returned by diff_gene_epigenetic_analysis()[[2]] for sample 2
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

#x #Data.frame of sample 1
#y #data.frame of sample 2
#z #Character vecotr of sample names
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

#x #list of bamfiles returned by bamfile()
#y #list of character vectors with motifs to group fragments by
#z #Name of groups
#t #title of plot
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
#x #Numeric vector of the motif proportions
entrop <- function(x){
    res <- -sum(x * log2(x))
    return(res)
}
#x #Proportion of fragment end motifs for sample 1
#y #Proportion of fragment end motifs for sample 2
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
#x #Proportion of fragment end motifs for sample 1
#y #Proportion of fragment end motifs for sample 2
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

#x #list of motifs to separate motifs into
#y #list of data.frames to to get information from
#z #character vector with names of samples
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

#x #named list of BAMfiles returned by bamfile() for group 1
#y #named list of BAMfiles returned by bamfile() for group 2
#z #vector of names for group 1 and 2
#p #Title of plot
#q If q = TRUE (default) concatenated lines are plotted. 
#if q = FALSE separate lines in same plot is plotted. 
#If q = "individual" individual plots for individual patients is plotted
fragment_length <- function(x,y,z,p,q){
    library(ggplot2)
    library(tidyr)
    library(plyr)
    len1 <- c()
    len2 <- c()
    name1 <- c()
    name2 <- c()
    for (i in 1:length(x)){
        len <- mcols(x[[i]])$isize
        len <- len[len >0]
        len <- len[!is.na(len)]
        len1 <- c(len1,len)
        nam <- unlist(lapply(strsplit(names(x)[i],"_"), "[[", 1))
        name <- rep(nam,length(len))
        name1 <- c(name1,name)
    }
    for (i in 1:length(y)){
        len <- mcols(y[[i]])$isize
        len <- len[len >0]
        len <- len[!is.na(len)]
        len2 <- c(len2,len)
        nam <- unlist(lapply(strsplit(names(y)[i],"_"), "[[", 1))
        name <- rep(nam,length(len))
        name2 <- c(name2,name)
    }
    rects <- data.frame(xstart = 0, xend = 150, col = "col")
    df <- data.frame(len = c(len1,len2),
                     Sample = factor(c(rep(z[1],length(len1)),
                                       rep(z[2],length(len2))),
                                     levels = c(z[1], z[2])),
                     group = c(name1,name2))
    if(isTRUE(q)){
        gg <- ggplot()+
            geom_rect(data = rects, aes(xmin = xstart, xmax = xend, 
                                        ymin = -Inf, ymax = Inf, fill = col),
                      alpha = 0.4)+
            scale_fill_manual(values = "green4",guide = "none")+
            geom_density(data = df, 
                         aes(x = len, color = Sample),
                         size = 1)+
            theme_bw(base_size = 15)+
            scale_color_manual(values = c("#6a00fc","#ffa10c"))+
            labs(title = p, y = "Density", x = "Fragment length (bp)")+
            scale_x_continuous(breaks = seq(0,400,50), limits = c(0,400))+
            th
    }
    if (isFALSE(q)){
        df <- df %>% mutate(group = paste(Sample,group))
        df <- df %>% mutate(group = factor(group, levels = unique(df$group)))
        gg <- ggplot()+
            geom_rect(data = rects, aes(xmin = xstart, xmax = xend, 
                                        ymin = -Inf, ymax = Inf, fill = col),
                      alpha = 0.4)+
            scale_fill_manual(values = "green4",guide = "none")+
            geom_density(data = df, 
                         aes(x = len, color = Sample, group = group),
                         size = 1)+
            theme_bw(base_size = 15)+
            scale_color_manual(values = c("#6a00fc","#ffa10c"))+
            labs(title = p, y = "Density", x = "Fragment length (bp)")+
            scale_x_continuous(breaks = seq(0,400,50), limits = c(0,400))+
            th
    }
    if(q == "individual"){
        gg <- ggplot()+
            geom_rect(data = rects, aes(xmin = xstart, xmax = xend, 
                                        ymin = -Inf, ymax = Inf, fill = col),
                      alpha = 0.4)+
            scale_fill_manual(values = "green4",guide = "none")+
            geom_density(data = df, 
                         aes(x = len, color = Sample),
                         size = 1)+
            theme_bw(base_size = 15)+
            scale_color_manual(values = c("#6a00fc","#ffa10c"))+
            labs(title = p, y = "Density", x = "Fragment length (bp)")+
            scale_x_continuous(breaks = seq(0,400,50), limits = c(0,400))+
            th+
            facet_wrap( ~ group,scales = "free_y")+
            theme(strip.background =element_rect(fill="white"))+
            theme(strip.text = element_text(colour = 'black',face = "bold"))
    }
    return(gg)
}

#x #BAM file returned by bamfile()
#y #cutoff for fragment length
proportion_sub <- function(x,y){
    library(Rsamtools)
    len <- mcols(x)$isize
    len <- len[len >0]
    len <- len[!is.na(len)]
    bam_lengths_sub150 <- len[len<y]
    res <- length(bam_lengths_sub150)/length(len)
    return(res)
}


#x #list of bam files returned by bamfile() for group 1
#y #cutoff for fragment length
#z #name of samples
proportion_sub_df <- function(x, y, z){
    res <- c()
    for (i in 1:length(x)){
        res[i] <- proportion_sub(x[[i]], y)
    }
    df <- data.frame(sample = z,
                     sub_fraction = res)
    return(df)
}

#x #list of proportion sub data.frame returned by proportion_sub_df() for samples
#z #Name of samples
#p #If samples are paired p = "Paired (default), otherwise p = "Unpaired"
#g #If active vs. inactive is analyzed g = "activity", default = NULL.
proportion_sub_boxplot <- function(x, z, p = "Paired", g = NULL){
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
                labs(title ="",
                     x = "",
                     y = "Fraction under 150 bp")+
                scale_fill_manual("Sample",
                                  values = c("#ffa10c","#6a00fc"))+
                stat_compare_means(method = "t.test",
                                   comparisons = list(c(z[1], z[2])))+
                theme_bw()+
                th+
                theme(legend.position = "none")
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
                labs(title = "",
                     x = "",
                     y = "Fraction under 150 bp")+
                scale_fill_manual("Sample",
                                  values = c("#ffa10c","firebrick","#6a00fc"))+
                stat_compare_means(method = "t.test",
                                   comparisons = list(c(z[1],z[2]),
                                                      c(z[1],z[3]),
                                                      c(z[2],z[3])))+
                theme_bw()+
                th+
                theme(legend.position = "none")
        }
    }
    else{
        if(is.null(g)){
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
                labs(title = "",
                     x = "",
                     y = "Fraction under 150 bp")+
                scale_fill_manual("Sample",
                                  values = c("#ffa10c","#6a00fc"))+
                stat_compare_means(method = "t.test",
                                   paired = TRUE,
                                   comparisons = list(c(z[1], z[2])))+
                theme_bw()+
                th+
                theme(legend.position = "none")
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
                labs(title = "",
                     x = "",
                     y = "Fraction under 150 bp")+
                scale_fill_manual("Gene activity",
                                  values = c("#ffa10c","#6a00fc"))+
                stat_compare_means(method = "t.test",
                                   paired = TRUE,
                                   comparisons = list(c(z[1], z[2])))+
                theme_bw()+
                th+
                theme(legend.position = "none")
        }
        
    }
    return(gg)
}
#x #Table consisting of BAM file names (name used when bamfile() is used),
  #genes mutated and the position of the mutation
#y #Named list of samples used
fragment_length_wt_ctdna_df <- function(x,y){
    library(Rsamtools)
    library(tidyr)
    library(GenomicAlignments)
    or_wd <- getwd()
    `%ni%` <- Negate(`%in%`)
    df <- read.table(x,header = T)
    df <- df[!is.na(df$genes),]
    df6 <- data.frame(len = c(),
                      Sample = c(),
                      group = c())
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
                                            rep("WT", length(wt_len))),
                                 group = rep(df$name[i],length(length(ct_len)+length(wt_len))))
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
                                             rep("WT", length(wt_len))),
                                  group = rep(df$name[i],length(length(ct_len)+length(wt_len))))
                df6 <- rbind(df6,df4)
                print(i)
            }
        }    
    }
    return(df6)
}

#x #Table consisting of BAM file names (name used when bamfile() is used),
#genes mutated and the position of the mutation
#y #Named list of samples used
fragment_length_wt_ctdna_vessies_df <- function(x,y){
    library(Rsamtools)
    library(tidyr)
    library(GenomicAlignments)
    or_wd <- getwd()
    `%ni%` <- Negate(`%in%`)
    df <- read.table(x,header = T)
    df <- df[!is.na(df$genes),]
    df6 <- data.frame(len = c(),
                      Sample = c(),
                      group = c())
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
            cig <- cig[!is.na(len)]
            len <- len[!is.na(len)]
            if(grepl("indel",type)){
                indel <- strsplit(type,split = "-")[[1]][2]
                ct_len <- abs(len[grepl(indel,cig)])
                wt_len <- abs(len[!grepl(indel,cig)])
                df1 <- data.frame(len = c(ct_len,wt_len),
                                  Sample = c(rep("ctDNA",length(ct_len)),
                                             rep("WT", length(wt_len))),
                                  group = rep(df$name[i],length(length(ct_len)+length(wt_len))))
                df6 <- rbind(df6,df1)
            }
            else{
                param = ScanBamParam(which = g, what = c("seq","mapq", "isize","pos", "cigar"),
                                     tag = "MD")
                setwd("D:/Lung cancer input/PosDeduped")
                bam <- scanBam(df$file[i],
                               param = param)
                setwd(or_wd)
                bam_pos <- as.numeric(unlist(bam[[1]][1]))
                bam_size <- as.numeric(unlist(bam[[1]][4]))
                bam_cigar <- unlist(bam[[1]][3])
                bam_MD <- unlist(bam[[1]][6])
                bam_seq <- as.character(unlist(bam[[1]][5])$seq)
                df1 <- data.frame(position = bam_pos,
                                  size = bam_size,
                                  cigar = bam_cigar,
                                  MD = bam_MD,
                                  seq = bam_seq)
                df1 <- df1[!is.na(df1$position),]
                alt <- strsplit(type,split = "-")[[1]][2]
                FragmentsWithVariants <- function(Pos, Alt,Bam){
                    base <- substring(Bam$seq, Pos-Bam$position,Pos-Bam$position)
                    lengths <- abs(Bam$size[base == Alt])
                    lengths
                }
                FragmentsWithoutVariants <- function(Pos, Alt,Bam){
                    base <- substring(Bam$seq, Pos-Bam$position,Pos-Bam$position)
                    lengths <- abs(Bam$size[base != Alt])
                    lengths
                }
                ct_len <- as.integer(unlist(mapply(FragmentsWithVariants,
                                                    Pos = mut_pos,
                                                    Alt = alt,
                                                    MoreArgs = list(Bam = df1))))
                wt_len <- as.integer(unlist(mapply(FragmentsWithoutVariants,
                                                       Pos = mut_pos,
                                                       Alt = alt,
                                                       MoreArgs = list(Bam = df1))))
                df4 <- data.frame(len = c(ct_len,wt_len),
                                  Sample = c(rep("ctDNA",length(ct_len)),
                                             rep("WT", length(wt_len))),
                                  group = rep(df$name[i],length(length(ct_len)+length(wt_len))))
                df6 <- rbind(df6,df4)
                print(i)
            }
        }    
    }
    df6 <- df6 %>% filter(len != 0)
    return(df6)
}
#fragment_length_wt_ctdna_input_vessies <- fragment_length_wt_ctdna_vessies_df("BL mutations input files.txt",named_cancer_input_list)

#x #data.frame returned by fragment_length_wt_ctdna_df()
#y If y = TRUE (default) concatenated lines are plotted. Otherwise the plot is per patient. 
fragment_length_wt_ctdna_plot <- function(x,y = TRUE){
    library(ggplot2)
    library(tidyr)
    library(plyr)
    if(y){
        colnames(x) <- c("len", "Sample", "group")
        x$Sample <- factor(x$Sample, levels = c("ctDNA", "WT"))
        rects <- data.frame(xstart = 0, xend = 150, col = "col")
        gg <- ggplot()+
            geom_rect(data = rects, aes(xmin = xstart, xmax = xend, 
                                        ymin = -Inf, ymax = Inf, fill = col),
                      alpha = 0.4)+
            scale_fill_manual(values = "green4",guide = "none")+
            geom_density(data = x, 
                         aes(x = len, color = Sample),
                         size = 1)+
            theme_bw(base_size = 15)+
            scale_color_manual("Fragment", values = c("#6a00fc","#ffa10c"))+
            labs(title = "", y = "Density", x = "Fragment length (bp)")+
            scale_x_continuous(breaks = seq(0,400,50), limits = c(0,400))+
            th
        
    }
    else{
        colnames(x) <- c("len", "Sample", "group")
        x$Sample <- factor(x$Sample, levels = c("ctDNA", "WT"))
        x <- x %>% mutate(group = unlist(lapply(strsplit(x$group,"_"), "[[", 1)))
        rects <- data.frame(xstart = 0, xend = 150, col = "col")
        gg <- ggplot()+
            geom_rect(data = rects, aes(xmin = xstart, xmax = xend, 
                                        ymin = -Inf, ymax = Inf, fill = col),
                      alpha = 0.4)+
            scale_fill_manual(values = "green4",guide = "none")+
            geom_density(data = x, 
                         aes(x = len, color = Sample),
                         size = 1)+
            theme_bw(base_size = 15)+
            scale_color_manual("Fragment", values = c("#6a00fc","#ffa10c"))+
            labs(title = "", y = "Density", x = "Fragment length (bp)")+
            scale_x_continuous(breaks = seq(0,400,50), limits = c(0,400))+
            th+
            facet_wrap( ~ group,scales = "free_y")+
            theme(strip.background =element_rect(fill="white"))+
            theme(strip.text = element_text(colour = 'black',face = "bold"))
    }
    return(gg) 
    
}



MAF_in_bins <- function(x){
    library(dplyr)
    intervals <- seq(50,400,10)
    colnames(x) <- c("len", "Sample","Prøve")
    df_wt <- x %>% dplyr::filter(Sample == "WT")
    df_ct <- x %>% dplyr::filter(Sample == "ctDNA")
    x$Sample <- factor(x$Sample, levels = c("ctDNA", "WT"))
    x_mark <- c()
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
    rects <- data.frame(xstart = c(50,150,200,340), xend = c(150,200,340,400), 
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
        labs(title ="",
             x = "Fragment length (bp)",
             y = "MAF enrichment")+
        geom_hline(yintercept = 1,
                   color = "black",
                   size = 1)+
        geom_smooth(data = df1, aes(x = len, y = enrichment),
                    method = "loess",
                    span = 0.5, se = F,
                    color = "black",
                    size = 1.2)+
        geom_vline(xintercept = c(50,150,200,340,400),
                   color = "black",
                   linetype = "dashed",
                   size = 0.5)+
        annotate(geom = "text",
                 x = c(58,160,210,350,410),
                 y = c(2.4,2.4,2.4,2.4,2.4),
                 label = c("50", "150", "200", "340", "400"),
                 color = "black",fontface = "bold",
                 size = 4)+
        scale_y_continuous(limits = c(0.4,2.4))+
        scale_x_continuous(limits = c(0,450),breaks = seq(0,450,50))+
        theme_bw()+
        th
    return(gg)
}



#x #data.frame returned by fragment_length_wt_ctdna_vessies_df()
fragment_length_wt_ctdna_boxplot <- function(x){
    ct_frac <- c()
    wt_frac <- c()
    for (i in 1:length(unique(x$group))){
        ddf <- x %>% filter(group == unique(x$group)[i])
        df_ct <-  ddf %>% filter(Sample == "ctDNA")
        df_wt <-  ddf %>% filter(Sample == "WT")
        df_ct_sub <- df_ct %>% filter(len < 150)
        df_wt_sub <- df_wt %>% filter(len < 150)
        ct_frac[i] <- nrow(df_ct_sub)/nrow(df_ct)
        wt_frac[i] <- nrow(df_wt_sub)/nrow(df_wt)
    }
    df <- data.frame(sub_fraction = c(wt_frac,ct_frac),
                     sample = factor(c(rep("WT",length(unique(x$group))),
                              rep("ctDNA", length(unique(x$group)))),
                                     levels = c("WT", "ctDNA")),
                     pairs = rep(1:length(unique(x$group)),2),
                     group = rep(unique(x$group),2))
    gg <- df %>% ggplot(aes(x=sample, y=sub_fraction, fill=sample)) + 
        geom_boxplot(outlier.alpha = 0)+
        geom_line(aes(group=pairs))+
        geom_point(alpha = 0.5)+
        labs(title = "",
             x = "",
             y = "Fraction under 150 bp")+
        scale_fill_manual("Fragments",
                          values = c("#ffa10c","#6a00fc"))+
        stat_compare_means(method = "t.test",
                           paired = TRUE,
                           comparisons = list(c("WT", "ctDNA")))+
        theme_bw()+
        th+
        theme(legend.position = "none")
    return(gg)
}

#x #enrichment data.frame 
#y #Granges object returned by gr() for the 197 AVENIO target genes
#z #name of sample to isolate from enrichment data.frame
#i #name of input BAM file
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
                                rep("Inactive",length(length_bottom15))),
                     group = rep(z,length(length_top15)+length(length_bottom15)))
    return(df)
}

#x #Data.frame with two columns. One with fragment lengths and the other determining whether the fragment originates from an active or inactive gene
#p #sample
#q If q = TRUE (default) concatenated lines are plotted. 
#if q = FALSE separate lines in same plot is plotted. 
#If q = "individual" individual plots for individual patients is plotted
fragment_length_active_inactive_plot <- function(x,p, q){
    library(ggplot2)
    library(tidyr)
    library(plyr)
    df <- x
    df <- df %>% mutate(Sample = ifelse(Sample == "Active", 
                                        paste(p, "High expressed"),
                                        paste(p, "Low expressed")))
    rects <- data.frame(xstart = 0, xend = 150, col = "col")
    if(isTRUE(q)){
        gg <- ggplot()+
            geom_rect(data = rects, aes(xmin = xstart, xmax = xend, 
                                        ymin = -Inf, ymax = Inf, fill = col),
                      alpha = 0.4)+
            scale_fill_manual(values = "green4",guide = "none")+
            geom_density(data = df, 
                         aes(x = len, color = Sample),
                         size = 1)+
            theme_bw(base_size = 15)+
            scale_color_manual(values = c("#6a00fc","#ffa10c"),
                               name = "Gene activity")+
            labs(title = "", y = "Density", x = "Fragment length (bp)")+
            scale_x_continuous(breaks = seq(0,400,50), limits = c(0,400))+
            th
    }
    if(q == "individual"){
        gg <- ggplot()+
            geom_rect(data = rects, aes(xmin = xstart, xmax = xend, 
                                        ymin = -Inf, ymax = Inf, fill = col),
                      alpha = 0.4)+
            scale_fill_manual(values = "green4",guide = "none")+
            geom_density(data = df, 
                         aes(x = len, color = Sample),
                         size = 1)+
            theme_bw(base_size = 15)+
            scale_color_manual(values = c("#6a00fc","#ffa10c"),
                               name = "Gene activity")+
            labs(title = "", y = "Density", x = "Fragment length (bp)")+
            scale_x_continuous(breaks = seq(0,400,50), limits = c(0,400))+
            th+
            facet_wrap( ~ group,scales = "free_y")+
            theme(strip.background =element_rect(fill="white"))+
            theme(strip.text = element_text(colour = 'black',face = "bold"))
    }
    
    return(gg)
}

#x #list of fragment length data.frames returned by fragment_length_active_inactive_df()
#y #fragment length cutoff
#z #name of samples
#a #Definition of whether "Active" or "Inactive" genes are analyzed
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
#x #enrichment data.frame 
#y #Granges object returned by gr() for the 197 AVENIO target genes
#z #name of sample to isolate from enrichment data.frame
#j #name of input BAM file
#p #cutoff for fragment length
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

#x #data.frame with quartiles, sub150 fraction of the quartiles and the sample
sub150_fraction_quartiles_plot <- function(x){
    library(ggpubr)
    library(rstatix)
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
    ddf_Q10 <- ddf %>% filter(quartile == "Q10")
    ddf <- ddf %>% 
        mutate(quartile = factor(quartile,
                                 level = c("Q1", "Q2", "Q3", "Q4",
                                           "Q5", "Q6", "Q7", "Q8",
                                           "Q9", "Q10")))
    res.aov <- ddf %>% friedman_test(norm_fraction ~ quartile | sample)
    gg <- ddf %>% ggplot(aes(x = quartile, y = norm_fraction))+
        geom_boxplot(outlier.alpha = 0,
                     fill = "#6a00fc")+
        geom_jitter(alpha = 0.5, width = 0.1)+
        theme_bw()+
        stat_compare_means(aes(label = after_stat(p.signif)),
                           method = "t.test",
                           paired = TRUE,
                           ref.group = "Q1",
                           size = 4.5)+
        geom_text(label = get_test_label(res.aov, type = "text"),
                  x = 8.2, y = 0.9,size = 4.5)+
        labs(title ="",
             x = "Quantiles",
             y = "Normalized fraction under 150 bp")+
        scale_y_continuous(limits = c(0.9,1.5))+
        th
    
    return(gg)
}

#x #enrichment data.frame 
#y #Granges object returned by gr() for the 197 AVENIO target genes
#z #name of sample to isolate from enrichment data.frame
#j #name of input BAM file
#m #number of bases to be evaluated in each fragment end
#a #character vector of motifs defined as "Active motifs"
active_motif_quartiles_df <- function(x,y,z,j,m,a){
    library(Rsamtools)
    df_sample <- x[c("genes",z)]
    df_sample <- df_sample[order(df_sample[,2]),]
    e <- df_sample
    number_high <- 19
    number_low <- 1
    res <- c()
    for (i in 1:10){
        print(i)
        if(i < 5){
            genes <- e$genes[number_low:number_high]
            number_low <- number_low + 19
            number_high <- number_high + 19
            which <- y[elementMetadata(y)[,1] %in% genes]
            index1 <- match(genes,elementMetadata(which)[,1])
            which <- which[index1]
            what <- c("seq", "strand")
            param <- ScanBamParam(which = which, what = what)
            bam <- scanBam(j, param = param)
            seqs <- unname(unlist(lapply(bam,function(x) x$seq)))
            strand <- unname(unlist(lapply(bam,function(x) x$strand)))
            grang <- GRanges(seqnames = "chr1", ranges = IRanges(start = 1,
                                                                       end = 2),
                                   strand = strand)
            mcols(grang)$seq <- DNAStringSet(do.call(c,unlist(seqs)))
            mers <- mer_count(grang,m)
            prop_mers <- prop_df(mers)
            mers_active <- prop_mers %>% filter(res %in% a)
            ress <- sum(mers_active$Freq)
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
            what <- c("seq", "strand")
            param <- ScanBamParam(which = which, what = what)
            bam <- scanBam(j, param = param)
            seqs <- unname(unlist(lapply(bam,function(x) x$seq)))
            strand <- unname(unlist(lapply(bam,function(x) x$strand)))
            grang <- GRanges(seqnames = "chr1", ranges = IRanges(start = 1,
                                                                 end = 2),
                             strand = strand)
            mcols(grang)$seq <- DNAStringSet(do.call(c,unlist(seqs)))
            mers <- mer_count(grang,m)
            prop_mers <- prop_df(mers)
            mers_active <- prop_mers %>% filter(res %in% a)
            ress <- sum(mers_active$Freq)
            res[i] <- ress
        }
    }
    df <- data.frame(quartile = c("Q1", "Q2", "Q3", "Q4",
                                  "Q5", "Q6", "Q7", "Q8",
                                  "Q9", "Q10"),
                     active_fraction = res,
                     sample = rep(z,10))
    return(df)
}
#a #Definition of whether "Active" or "Inactive" genes are analyzed

active_motif_fraction_quartiles_plot <- function(x,a){
    library(ggpubr)
    library(dplyr)
    library(rstatix)
    ddf <- data.frame(quartile = c(),
                      active_fraction = c(),
                      sample = c(),
                      norm_fraction = c())
    for(i in 1:length(unique(x$sample))){
        sam <- unique(x$sample)[i]
        df <- x %>% filter(sample == sam)
        base_fraction <- df$active_fraction[1]
        df <- df %>% mutate(norm_fraction = active_fraction/base_fraction)
        ddf <- rbind(ddf,df)
    }
    if(a == "Active"){
        lab <- "Normalized cfChIP motif fraction"
    }
    else{
        lab <- "Normalized input motif fraction"
    }
    ddf <- ddf %>% 
        mutate(quartile = factor(quartile,
                                 level = c("Q1", "Q2", "Q3", "Q4",
                                           "Q5", "Q6", "Q7", "Q8",
                                           "Q9", "Q10")))
    res.aov <- ddf %>% friedman_test(norm_fraction ~ quartile | sample)
    aov_res <- res.aov$p
    if(aov_res < 0.05){
        gg <- ddf %>% ggplot(aes(x = quartile, y = norm_fraction))+
            geom_boxplot(outlier.alpha = 0,
                         fill = "#6a00fc")+
            geom_jitter(alpha = 0.5, width = 0.1)+
            theme_bw()+
            geom_text(label = get_test_label(res.aov, type = "text"),
                      x = 8.2, y = 0.8,size = 4.5)+
            stat_compare_means(aes(label = after_stat(p.signif)),
                               method = "t.test",
                               paired = TRUE,
                               ref.group = "Q1",
                               size = 4.5)+
            labs(title="",
                 x = "Quantiles",
                 y = lab)+
            scale_y_continuous(limits = c(0.8,1.22))+
            th
    }
    else{
        gg <- ddf %>% ggplot(aes(x = quartile, y = norm_fraction))+
            geom_boxplot(outlier.alpha = 0,
                         fill = "#6a00fc")+
            geom_jitter(alpha = 0.5, width = 0.1)+
            theme_bw()+
            geom_text(label = get_test_label(res.aov, type = "text"),
                      x = 8.2, y = 0.8,size = 4.5)+
            labs(title="",
                 x = "Quantiles",
                 y = lab)+
            scale_y_continuous(limits = c(0.8,1.05))+
            th
    }
    
    
    return(gg)
}
active_motif_fraction_quartiles_plot(collected_inactive_inactive_cancer_fraction_quartiles, "Inactive")


#a #Definition of whether "Active" or "Inactive" genes are analyzed

high_motif_fraction_quartiles_plot <- function(x,a){
    library(ggpubr)
    library(rstatix)
    ddf <- data.frame(quartile = c(),
                      active_fraction = c(),
                      sample = c(),
                      norm_fraction = c())
    for(i in 1:length(unique(x$sample))){
        sam <- unique(x$sample)[i]
        df <- x %>% filter(sample == sam)
        base_fraction <- df$active_fraction[1]
        df <- df %>% mutate(norm_fraction = active_fraction/base_fraction)
        ddf <- rbind(ddf,df)
    }
    if(a == "Active"){
        lab <- "Normalized high expression motif fraction"
        posi <- 0.6
    }
    else{
        lab <- "Normalized low expression motif fraction"
        posi <- 1.3
    }
    ddf <- ddf %>% 
        mutate(quartile = factor(quartile,
                                 level = c("Q1", "Q2", "Q3", "Q4",
                                           "Q5", "Q6", "Q7", "Q8",
                                           "Q9", "Q10")))
    res.aov <- ddf %>% friedman_test(norm_fraction ~ quartile | sample)
    gg <- ddf %>% ggplot(aes(x = quartile, y = norm_fraction))+
        geom_boxplot(outlier.alpha = 0,
                     fill = "#6a00fc")+
        geom_jitter(alpha = 0.5, width = 0.1)+
        theme_bw()+
        geom_text(label = get_test_label(res.aov, type = "text"),
                  x = 8.2, y = posi,size = 4.5)+
        stat_compare_means(aes(label = after_stat(p.signif)),
                           method = "t.test",
                           paired = TRUE,
                           ref.group = "Q1",
                           size = 4.5)+
        labs(title="",
             x = "Quantiles",
             y = lab)+
        scale_y_continuous(limits = c(0.6,1.4))+
        th
    
    return(gg)
}

#x #enrichment data.frame 
#y #Granges object returned by gr() for the 197 AVENIO target genes
#z #name of sample to isolate from enrichment data.frame
#j #name of input BAM file
#p #cutoff for fragment length
#m #number of bases to be evaluated in each fragment end
#a #character vector of motifs defined as "Active motifs"
sub150_active_motif_fraction_quartiles_df <- function(x,y,z,j,p,m,a){
    library(Rsamtools)
    df_sample <- x[c("genes",z)]
    df_sample <- df_sample[order(df_sample[,2]),]
    e <- df_sample
    number_high <- 19
    number_low <- 1
    res <- c()
    for (i in 1:10){
        print(i)
        if(i < 5){
            genes <- e$genes[number_low:number_high]
            number_low <- number_low + 19
            number_high <- number_high + 19
            which <- y[elementMetadata(y)[,1] %in% genes]
            index1 <- match(genes,elementMetadata(which)[,1])
            which <- which[index1]
            what <- c("seq", "strand", "isize")
            param <- ScanBamParam(which = which, what = what)
            bam <- scanBam(j, param = param)
            seqs <- unname(unlist(lapply(bam,function(x) x$seq)))
            strand <- unname(unlist(lapply(bam,function(x) x$strand)))
            len <- unname(unlist(lapply(bam,function(x) x$isize)))
            grang <- GRanges(seqnames = "chr1", ranges = IRanges(start = 1,
                                                                 end = 2),
                             strand = strand)
            mcols(grang)$len <- len
            mcols(grang)$seq <- DNAStringSet(do.call(c,unlist(seqs)))
            grang <- grang[!is.na(mcols(grang)$len)]
            grang <- grang[mcols(grang)$len > 0]
            grang_sub <- grang[mcols(grang)$ len < p]
            mers <- mer_count(grang_sub,m)
            mers_active <- mers[names(mers) %in% a]
            ress <- sum(mers_active)/length(strand(grang))
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
            what <- c("seq", "strand", "isize")
            param <- ScanBamParam(which = which, what = what)
            bam <- scanBam(j, param = param)
            seqs <- unname(unlist(lapply(bam,function(x) x$seq)))
            strand <- unname(unlist(lapply(bam,function(x) x$strand)))
            len <- unname(unlist(lapply(bam,function(x) x$isize)))
            grang <- GRanges(seqnames = "chr1", ranges = IRanges(start = 1,
                                                                 end = 2),
                             strand = strand)
            mcols(grang)$len <- len
            mcols(grang)$seq <- DNAStringSet(do.call(c,unlist(seqs)))
            grang <- grang[!is.na(mcols(grang)$len)]
            grang <- grang[mcols(grang)$len > 0]
            grang_sub <- grang[mcols(grang)$ len < p]
            mers <- mer_count(grang_sub,m)
            mers_active <- mers[names(mers) %in% a]
            ress <- sum(mers_active)/length(strand(grang))
            res[i] <- ress
        }
    }
    df <- data.frame(quartile = c("Q1", "Q2", "Q3", "Q4",
                                  "Q5", "Q6", "Q7", "Q8",
                                  "Q9", "Q10"),
                     fraction = res,
                     sample = rep(z,10))
    return(df)
}

#a #Definition of whether "Active" or "Inactive" genes are analyzed

sub150_active_motif_fraction_quartiles_plot <- function(x,a){
    library(ggpubr)
    library(rstatix)
    ddf <- data.frame(quartile = c(),
                      fraction = c(),
                      sample = c(),
                      norm_fraction = c())
    for(i in 1:length(unique(x$sample))){
        sam <- unique(x$sample)[i]
        df <- x %>% filter(sample == sam)
        base_fraction <- df$fraction[1]
        df <- df %>% mutate(norm_fraction = fraction/base_fraction)
        ddf <- rbind(ddf,df)
    }
    if(a == "Active"){
        lab <- "Normalized sub 150 cfChIP motif fraction"
    }
    else{
        lab <- "Normalized sub 150 inactive motif fraction"
    }
    dfdf <- ddf %>% filter(quartile == "Q10")
    ddf <- ddf %>% 
        mutate(quartile = factor(quartile,
                                 level = c("Q1", "Q2", "Q3", "Q4",
                                           "Q5", "Q6", "Q7", "Q8",
                                           "Q9", "Q10")))
    res.aov <- ddf %>% friedman_test(norm_fraction ~ quartile | sample)
    gg <- ddf %>% ggplot(aes(x = quartile, y = norm_fraction))+
        geom_boxplot(outlier.alpha = 0,
                     fill = "#6a00fc")+
        geom_jitter(alpha = 0.5, width = 0.1)+
        theme_bw()+
        geom_text(label = get_test_label(res.aov, type = "text"),
                  x = 8.2, y = 0.8,size = 4.5)+
        stat_compare_means(aes(label = after_stat(p.signif)),
                           method = "t.test",
                           paired = TRUE,
                           ref.group = "Q1",
                           size = 4.5)+
        labs(title = "",
             x = "Quantiles",
             y = lab)+
        scale_y_continuous(limits = c(0.8,1.52))+
        th
    
    return(gg)
}


#x enrichment_df
#y name of sample to be analyzed

cancer_upregulated_genes <- function(x,y){
    library(dplyr)
    healthy <- enrichment_df %>% dplyr::select(c(genes,grep("HC",colnames(enrichment_df))))
    rownames(healthy) <- healthy$genes
    healthy <- healthy %>% dplyr::select(-genes)
    avg_h <- rowMeans(healthy)
    healthy <- data.frame(genes = names(avg_h),
                          enrichment = as.numeric(avg_h))
    case <- enrichment_df %>% dplyr::select(all_of(c("genes",y)))
    colnames(case) <- c("genes", "enrichment")
    logFC <- log2(case$enrichment/healthy$enrichment)
    fc_df <- data.frame(genes = case$genes,
                        FC = logFC) %>% arrange(-FC)
    up_genes <- fc_df$genes[1:20]
    steady <- fc_df %>% arrange(abs(FC))
    steady_genes <- steady$genes[1:20]
    res_df <- data.frame(genes = c(up_genes,steady_genes),
                         type = c(rep("upregulated",length(up_genes)),
                                  rep("steady",length(steady_genes))))
    return(res_df)
}

#x enrichment_df
#y name of sample to be analyzed
top_bottom_genes <- function(x,y){
    library(dplyr) 
    case <- enrichment_df %>% dplyr::select(all_of(c("genes",y)))
    colnames(case) <- c("genes", "enrichment")
    case <- case %>% arrange(enrichment)
    case_bottom <- case$genes[1:15]
    case <- case %>% arrange(-enrichment)
    case_top <- case$genes[1:15]
    res_df <- data.frame(genes = c(case_top,case_bottom),
                         type = c(rep("upregulated",length(case_top)),
                                  rep("steady",length(case_bottom))))
    return(res_df)
}

#x data.frame of upregulated and steady genes returned by cancer_upregulated_genes()
#y Granges object returned by bamfile() of the input sample
#g Granges object of AVENIO genes returned by gr()
#z Granges object of nucleosome positions of cfDNA in healthy individual in AVENIO genes
#t Granges object of AVENIO target regions
cleavage_nuc_dist_data <- function(x,y,g,z,t){
    y <- y[subjectHits(findOverlaps(t,y, ignore.strand = T))]
    up_x <- x %>% filter(type == "upregulated")
    up_genes <- up_x$genes
    grang_up <- g[mcols(g)$SYMBOL %in% up_genes]
    up_overlap <- findOverlaps(grang_up,y,ignore.strand = T)
    up_reads <- y[subjectHits(up_overlap)]
    qur <- queryHits(up_overlap)
    mcols(up_reads)$gene <- mcols(grang_up[queryHits(up_overlap)])$SYMBOL  
    up_reads <- up_reads[!is.na(mcols(up_reads)$isize)]
    up_reads <- up_reads[mcols(up_reads)$isize > 0]
    mcols(up_reads)$end <- mcols(up_reads)$pos + mcols(up_reads)$isize
    df <- data.frame(start_pos = c(),
                     end_pos = c(),
                     gene = c(),
                     type = c())
    for(i in 1:length(up_genes)){
        print(up_genes[i])
        up_reads_g <- up_reads[mcols(up_reads)$gene == up_genes[i]]
        up_nuc_g <- z[mcols(z)$gene == up_genes[i]]
        close_nuc <- vector(length = length(mcols(up_reads_g)$pos))
        for(j in 1:length(mcols(up_nuc_g)$pos)){
            res <- mcols(up_reads_g)$pos - 
                   mcols(up_nuc_g)$pos[j] <= 0 & close_nuc == FALSE
            if(any(res)){
                close_nuc[res] <- max(c(mcols(up_reads_g)$pos[res],
                                  mcols(up_nuc_g)$pos[j]))
            }
        }
        ddf <- data.frame(start_pos = mcols(up_reads_g)$pos - close_nuc,
                          end_pos = mcols(up_reads_g)$end - close_nuc,
                          gene = mcols(up_reads_g)$gene,
                          type = rep("upregulated",length(mcols(up_reads_g)$pos)))
        ddf <- ddf %>% filter(end_pos > 0)
        df <- rbind(df,ddf)
    }
    steady_x <- x %>% filter(type == "steady")
    steady_genes <- steady_x$genes
    grang_steady <- g[mcols(g)$SYMBOL %in% steady_genes]
    steady_overlap <- findOverlaps(grang_steady,y,ignore.strand = T)
    steady_reads <- y[subjectHits(steady_overlap)]
    qur <- queryHits(steady_overlap)
    mcols(steady_reads)$gene <- mcols(grang_steady[queryHits(steady_overlap)])$SYMBOL  
    steady_reads <- steady_reads[!is.na(mcols(steady_reads)$isize)]
    steady_reads <- steady_reads[mcols(steady_reads)$isize > 0]
    mcols(steady_reads)$end <- mcols(steady_reads)$pos + mcols(steady_reads)$isize
    for(i in 1:length(steady_genes)){
        print(steady_genes[i])
        steady_reads_g <- steady_reads[mcols(steady_reads)$gene == steady_genes[i]]
        steady_nuc_g <- z[mcols(z)$gene == steady_genes[i]]
        close_nuc <- vector(length = length(mcols(steady_reads_g)$pos))
        for(j in 1:length(mcols(steady_nuc_g)$pos)){
            res <- mcols(steady_reads_g)$pos - 
                mcols(steady_nuc_g)$pos[j] <= 0 & close_nuc == FALSE
            if(any(res)){
                close_nuc[res] <- max(c(mcols(steady_reads_g)$pos[res],
                                        mcols(steady_nuc_g)$pos[j]))
            }
        }
        ddf <- data.frame(start_pos = mcols(steady_reads_g)$pos - close_nuc,
                          end_pos = mcols(steady_reads_g)$end - close_nuc,
                          gene = mcols(steady_reads_g)$gene,
                          type = rep("steady",length(mcols(steady_reads_g)$pos)))
        ddf <- ddf %>% filter(end_pos > 0)
        df <- rbind(df,ddf)
    }
    f_df <- data.frame(breaks = c(df$start_pos, df$end_pos),
                     gene = rep(df$gene,2),
                     type = rep(df$type,2))
    return(f_df)
    
}

#x data.frame of upregulated and steady genes returned by cancer_upregulated_genes()
#y Granges object returned by bamfile() of the input sample
#g Granges object of AVENIO genes returned by gr()
#z Granges object of nucleosome positions of cfDNA in healthy individual in AVENIO genes
#t Granges object of AVENIO target regions
bed_based_cleavage_nuc_dist_data <- function(x,y,g,z,t){
    y <- y[subjectHits(findOverlaps(t,y, ignore.strand = T))]
    up_x <- x %>% filter(type == "upregulated")
    up_genes <- up_x$genes
    grang_up <- g[mcols(g)$SYMBOL %in% up_genes]
    up_overlap <- findOverlaps(grang_up,y,ignore.strand = T)
    up_reads <- y[subjectHits(up_overlap)]
    qur <- queryHits(up_overlap)
    mcols(up_reads)$gene <- mcols(grang_up[queryHits(up_overlap)])$SYMBOL  
    up_reads <- up_reads[!is.na(mcols(up_reads)$isize)]
    up_reads <- up_reads[mcols(up_reads)$isize > 0]
    mcols(up_reads)$end <- mcols(up_reads)$pos + mcols(up_reads)$isize
    up_fragments <- GRanges(seqnames = seqnames(up_reads),
                            ranges = IRanges(start = mcols(up_reads)$pos,
                                             end = mcols(up_reads)$end),
                            strand = strand(up_reads))
    mcols(up_fragments)$gene <- mcols(up_reads)$gene
    nuc_grange <- GRanges(seqnames = seqnames(z),
                          ranges = IRanges(start = (mcols(z)$pos),
                                           end = mcols(z)$pos),
                          strand = strand(z))
    mcols(nuc_grange)$gene <- z$gene
    res <- findOverlaps(nuc_grange,up_fragments, 
                 ignore.strand = TRUE,
                 type = "any",
                 select = "all")
    start_dist <- start(up_fragments)[subjectHits(res)] - end(nuc_grange)[queryHits(res)]
    end_dist <- end(up_fragments)[subjectHits(res)] - end(nuc_grange)[queryHits(res)]
    df1 <- data.frame(start_dist = start_dist,
                      end_dist = end_dist,
                      frag_id = as.character(ranges(up_fragments)[subjectHits(res)]),
                     gene = mcols(up_fragments)$gene[subjectHits(res)],
                     type= rep("upregulated",length(mcols(up_fragments)$gene[subjectHits(res)])))
    print("before up slice")
    df_start <- df1 %>% 
        group_by(frag_id) %>% 
        dplyr::slice(which.max(start_dist)) %>% 
        dplyr::select(-end_dist)
    print("up start slice done")
    df_end <- df1 %>% 
        group_by(frag_id) %>% 
        dplyr::slice(which.min(end_dist)) %>% 
        dplyr::select(-start_dist)
    print("up end slice done")
    df2 <- data.frame(breaks = c(df_start$start_dist,
                                 df_end$end_dist),
                      breakage = c(rep("start",nrow(df_start)),
                                 rep("end",nrow(df_end))),
                      frag_id = df_start$frag_id,
                      gene = df_start$gene,
                      type = df_start$type)
    steady_x <- x %>% filter(type == "steady")
    steady_genes <- steady_x$genes
    grang_steady <- g[mcols(g)$SYMBOL %in% steady_genes]
    steady_overlap <- findOverlaps(grang_steady,y,ignore.strand = T)
    steady_reads <- y[subjectHits(steady_overlap)]
    qur <- queryHits(steady_overlap)
    mcols(steady_reads)$gene <- mcols(grang_steady[queryHits(steady_overlap)])$SYMBOL  
    steady_reads <- steady_reads[!is.na(mcols(steady_reads)$isize)]
    steady_reads <- steady_reads[mcols(steady_reads)$isize > 0]
    mcols(steady_reads)$end <- mcols(steady_reads)$pos + mcols(steady_reads)$isize
    steady_fragments <- GRanges(seqnames = seqnames(steady_reads),
                            ranges = IRanges(start = mcols(steady_reads)$pos,
                                             end = mcols(steady_reads)$end),
                            strand = strand(steady_reads))
    mcols(steady_fragments)$gene <- mcols(steady_reads)$gene
    nuc_grange <- GRanges(seqnames = seqnames(z),
                          ranges = IRanges(start = (mcols(z)$pos),
                                           end = mcols(z)$pos),
                          strand = strand(z))
    mcols(nuc_grange)$gene <- z$gene
    res <- findOverlaps(nuc_grange,steady_fragments, 
                        ignore.strand = TRUE,
                        type = "any",
                        select = "all")
    start_dist <- start(steady_fragments)[subjectHits(res)] - end(nuc_grange)[queryHits(res)]
    end_dist <- end(steady_fragments)[subjectHits(res)] - end(nuc_grange)[queryHits(res)]
    df1 <- data.frame(start_dist = start_dist,
                      end_dist = end_dist,
                      frag_id = as.character(ranges(steady_fragments)[subjectHits(res)]),
                      gene = mcols(steady_fragments)$gene[subjectHits(res)],
                      type= rep("steady",length(mcols(steady_fragments)$gene[subjectHits(res)])))
    print("before steady slice")
    df_start <- df1 %>% 
        group_by(frag_id) %>% 
        dplyr::slice(which.max(start_dist)) %>% 
        dplyr::select(-end_dist)
    print("steady start slice done")
    df_end <- df1 %>% 
        group_by(frag_id) %>% 
        dplyr::slice(which.min(end_dist)) %>% 
        dplyr::select(-start_dist)
    print("steady end slice done")
    df3 <- data.frame(breaks = c(df_start$start_dist,
                                 df_end$end_dist),
                      breakage = c(rep("start",nrow(df_start)),
                                   rep("end",nrow(df_end))),
                      frag_id = df_start$frag_id,
                      gene = df_start$gene,
                      type = df_start$type)
    df <- rbind(df2,df3)
    return(df)
    
}


#x data.frame of distances from nucleosome in steady genes and upregulated genes
#y if TRUE (default) plottes alle distancer som et samlet plot, hvis FALSE plottes distancer for hver patient
distance_nucleosome_plot <- function(x, y = TRUE){
    library(ggplot2)
    x <- x %>% mutate(type2 = ifelse(type == "upregulated", "Cancer High expressed", "Cancer Low expressed")) %>%
        mutate(type2 = factor(type2, level = c("Cancer High expressed","Cancer Low expressed")))
    rects <- data.frame(xstart = c(-95,-75,75), xend = c(-75,75,95), 
                        col = factor(c("Linker","Core","Linker"),
                                     levels = c("Linker", "Core")))
    if(y){
        gg <- ggplot()+
            geom_rect(data = rects, aes(xmin = xstart, xmax = xend, 
                                        ymin = -Inf, ymax = Inf, fill = col),
                      alpha = 0.5)+
            scale_fill_manual("Region",values = c("firebrick","Green4"))+ 
            geom_density(data = x, aes(x = breaks, color = type2),
                         size = 1)+
            scale_color_manual(name = "Gene activity",
                               values = c("Cancer High expressed" = "#6a00fc",
                                          "Cancer Low expressed" = "#ffa10c"))+ 
            scale_x_continuous(breaks = seq(-150,150,50), limits = c(-150,150))+
            theme_bw()+
            geom_vline(xintercept = c(-75,75),color="black",
                       linetype="dashed", size = 1)+
            coord_cartesian(ylim = c(0.0025, 0.005))+
            labs(title = "",
                 x = "Cleavage from dyad (bp)",
                 y = "Density")+
            th+
            guides(fill = guide_legend(override.aes = list(alpha = 0.5)))
    }
    else{
        gg <- ggplot()+
            geom_rect(data = rects, aes(xmin = xstart, xmax = xend, 
                                        ymin = -Inf, ymax = Inf, fill = col), 
                      alpha = 0.5)+
            scale_fill_manual("Region",values = c("firebrick","Green4"))+
            geom_density(data = x, aes(x = breaks, color = type2),
                         size = 1)+
            scale_color_manual(name = "Gene activity",
                            values = c("Cancer High expressed" = "#6a00fc",
                                       "Cancer Low expressed" = "#ffa10c"))+ 
            scale_x_continuous(breaks = seq(-150,150,50), limits = c(-150,150))+
            theme_bw()+
            geom_vline(xintercept = c(-75,75),color="black",
                       linetype="dashed", size = 1)+
            labs(title = "",
                 x = "Cleavage from dyad (bp)",
                 y = "Density")+
            facet_wrap(~ sample, nrow = 3)+
            coord_cartesian(ylim = c(0.0025, 0.005))+
            th+
            guides(fill = guide_legend(override.aes = list(alpha = 0.5)))+
            theme(strip.background =element_rect(fill="white"))+
            theme(strip.text = element_text(colour = 'black',face = "bold"))
    }
   return(gg) 
}

#x data.frame of distances from nucleosome in steady genes and upregulated genes

distance_nucleosome_boxplot <- function(x){
    box_df <- data.frame(upregulated = c(),
                         steady = c())
    for (i in 1:length(unique(x$sample))){
        sample_df <- x %>% filter(sample == unique(sample)[i])
        sample_up <- sample_df %>% filter(type == "upregulated")
        sample_steady <- sample_df %>% filter(type == "steady")
        sample_up_nuc <- sample_up %>% filter(breaks < 75 & -75 < breaks)
        sample_steady_nuc <- sample_steady %>% filter(breaks < 75 & -75 < breaks)
        sample_up_link <- sample_up %>% filter(breaks <= 95 & 75 <= breaks | breaks <= -75 & -95 <= breaks)
        sample_steady_link <- sample_steady %>% filter(breaks <= 95 & 75 <= breaks | breaks <= -75 & -95 <= breaks)
        df <- data.frame(upregulated = (nrow(sample_up_nuc)/nrow(sample_up_link)),
                         steady = (nrow(sample_steady_nuc)/nrow(sample_steady_link)))
        box_df <- rbind(box_df, df)
        
    }
    ddf <- data.frame(fraction = c(box_df$upregulated,box_df$steady),
                      gene_activity = factor(c(rep("Cancer High expressed",nrow(box_df)),
                                        rep("Cancer Low expressed",nrow(box_df))),
                                        level = c("Cancer Low expressed", "Cancer High expressed")),
                      pairs = c(rep(seq(1:nrow(box_df)),2)))
    gg <- ggplot(ddf,aes(x = gene_activity, y = fraction, fill = gene_activity))+
        geom_boxplot(outlier.alpha = 0)+
        geom_line(aes(group=pairs))+
        geom_point(alpha = 0.5)+
        labs(title = "",
             x = "",
             y = "Core/linker break fraction")+
        scale_fill_manual("Gene activity",
                          values = c("#ffa10c","#6a00fc"))+
        stat_compare_means(method = "t.test",
                           paired = TRUE,
                           comparisons = list(c("Cancer High expressed", "Cancer Low expressed")))+
        theme_bw()+
        th+
        theme(legend.position = "none")
    return(gg)
}

#KPRP and DLGAP2 have been put in manually. DLGAP2 is based on information given by Roche
#KPRP is made by manually looking on the HCC827 input file and determine where the coverage > 100
#x Name of bedgraph file used to study the coverages. Created using samtools. I use HCC827 input file
#y cutoff for minimum coverage value.
#z Name of bed file with nucleosome positions

coverages_nuc <- function(x,y = 100,z) {
    or_wd <- getwd()
    library(hwglabr2)
    library(GenomicRanges)
    library(dplyr)
    `%ni%` <- Negate(`%in%`)
    bg <- import_bedGraph(x)
    subsets <- bg[elementMetadata(bg)[,1]>y]
    reduced_subsets <- reduce(subsets)
    setwd("C:/Users/chris/OneDrive/1PhD/Fragmentering/endemotiver")
    CH01_nuc <- read.table(z,header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="")
    CH01_nuc_range <- GRanges(seqnames = CH01_nuc$V1,strand = CH01_nuc$V6,
                              ranges = IRanges(start = CH01_nuc$V2,
                                               end = CH01_nuc$V3))
    mcols(CH01_nuc_range)$pos <- CH01_nuc$V8
    genes <- vector(length = length(mcols(CH01_nuc_range)$pos)) 
    for (i in 1:length(grs$SYMBOL)){
        rang <- grs[i]
        res <- findOverlaps(rang,CH01_nuc_range,ignore.strand = T)
        genes[subjectHits(res)] <- rep(mcols(rang)$SYMBOL[1], length(subjectHits(res)))
    }
    mcols(CH01_nuc_range)$gene <- genes
    res <- subsetByOverlaps(CH01_nuc_range, reduced_subsets, ignore.strand = TRUE)
    setwd(or_wd)
    return(res)
}
#x #name of input BAM file
#a #character vector of motifs defined as "Active motifs"
#i #character vector of motifs defined as "Input motifs"
length_fragment_with_motif_df <- function(x,a,i){
    library(Biostrings)
    library(tidyr)
    library(plyr)
    positve <- x[strand(x) == "+"]
    str_set <- BStringSet(mcols(positve)$seq)
    mer <- extractAt(x = str_set, at = IRanges(start = 1,width = 3))
    active_mers <- mer %in% a
    active_pos_fragments <- positve[unlist(active_mers)]
    input_mers <- mer %in% i
    input_pos_fragments <- positve[unlist(input_mers)]
    negative <- x[strand(x) == "-"]
    rev_strands <- reverseComplement(mcols(negative)$seq)
    str_set <- BStringSet(rev_strands)
    mer <- extractAt(x = str_set, at = IRanges(start = 1,width = 3))
    active_mers <- mer %in% a
    active_neg_fragments <- negative[unlist(active_mers)]
    input_mers <- mer %in% i
    input_neg_fragments <- negative[unlist(input_mers)]
    active_neg_fragments <- active_neg_fragments[!is.na(mcols(active_neg_fragments)$isize)]
    active_neg_fragments <- active_neg_fragments[!mcols(active_neg_fragments)$isize == 0]
    active_pos_fragments <- active_pos_fragments[!is.na(mcols(active_pos_fragments)$isize)]
    active_pos_fragments <- active_pos_fragments[!mcols(active_pos_fragments)$isize == 0]
    input_neg_fragments <- input_neg_fragments[!is.na(mcols(input_neg_fragments)$isize)] 
    input_neg_fragments <- input_neg_fragments[!mcols(input_neg_fragments)$isize == 0]
    input_pos_fragments <- input_pos_fragments[!is.na(mcols(input_pos_fragments)$isize)]
    input_pos_fragments <- input_pos_fragments[!mcols(input_pos_fragments)$isize == 0]
    df <- data.frame(len = c(mcols(active_pos_fragments)$isize,
                                abs(mcols(active_neg_fragments)$isize),
                                mcols(input_pos_fragments)$isize,
                                abs(mcols(input_neg_fragments)$isize)),
                     Sample = factor(c(rep("cfChIP",
                                  length(active_pos_fragments)+length(active_neg_fragments)),
                              rep("Input",
                                  length(input_pos_fragments)+length(input_neg_fragments))),
                              levels = c("cfChIP", "Input")))
    return(df)
}

#data.frame of reads either with input or cfChIP end motif

length_fragment_with_motif_plot <- function(x){
    library(Biostrings)
    library(tidyr)
    library(plyr)
    df <- x
    df <- df %>% mutate(Sample = paste(Sample,"motif"))
    rects <- data.frame(xstart = 0, xend = 150, col = "col")
    gg <- ggplot()+
        geom_rect(data = rects, aes(xmin = xstart, xmax = xend, 
                                    ymin = -Inf, ymax = Inf, fill = col),
                  alpha = 0.4)+
        scale_fill_manual(values = "green4",guide = "none")+
        geom_density(data = df, 
                     aes(x = len, color = Sample),
                     size = 1)+
        theme_bw(base_size = 15)+
        scale_color_manual("Fragment", values = c("#6a00fc","#ffa10c"))+
        labs(title = "", y = "Density", x = "Fragment length (bp)")+
        scale_x_continuous(breaks = seq(0,400,50), limits = c(0,400))+
        th
    return(gg)
}


#data.frame of reads either with input or cfChIP end motif

length_fragment_with_motif_boxplot <- function(x){
    library(ggpubr)
    cfChIP_res <- x %>% filter(Sample == "cfChIP")
    input_res <- x %>% filter(Sample == "Input")
    df <- data.frame(proportions = c())
    for (i in 1:length(unique(x$patient_ID))){
        sample <- unique(x$patient_ID)[i]
        cfChIP_res_sample <- cfChIP_res %>%  filter(patient_ID == sample)
        input_res_sample <- input_res %>% filter(patient_ID == sample)
        cfChIP_sub <- cfChIP_res_sample %>% filter(len < 150)
        input_sub <- input_res_sample %>% filter(len < 150)
        cfChIP_prop <- nrow(cfChIP_sub)/nrow(cfChIP_res_sample)
        input_prop <- nrow(input_sub)/nrow(input_res_sample)
        df1 <- data.frame(proportions = c(cfChIP_prop,
                                          input_prop),
                          motif = factor(c("cfChIP motifs",
                                           "Input motifs"),
                                         levels = c("Input motifs",
                                                    "cfChIP motifs")),
                          relative_prop = c(cfChIP_prop/input_prop,
                                            input_prop/input_prop),
                          pairs = rep(i,2))
        df <- rbind(df,df1)
    }
    gg <- ggplot(df, aes(x = motif, y = relative_prop,
                         fill = motif))+
        geom_boxplot(outlier.alpha = 0)+
        geom_line(aes(group=pairs))+
        geom_point(alpha = 0.5)+
        labs(title = "",
             x = "",
             y = "Normalized fraction under 150 bp")+
        scale_fill_manual("Fragment motif",
                          values = c("#ffa10c","#6a00fc"))+
        stat_compare_means(method = "t.test",
                           paired = TRUE,
                           comparisons = list(c("cfChIP motifs", "Input motifs")))+
        theme_bw()+
        th+
        theme(legend.position = "none")
        
    return(gg)
}

length_fragment_with_motif_plot_zoom <- function(x){
    library(Biostrings)
    library(tidyr)
    library(plyr)
    df <- x
    rects <- data.frame(xstart = 50, xend = 150, col = "col")
    df <- df %>% mutate(Sample = paste(Sample,"motif")) %>% 
        mutate(group = paste(Sample,patient_ID)) %>% 
        filter(len < 156)
    gg <- ggplot()+
        geom_rect(data = rects, aes(xmin = xstart, xmax = xend, 
                                    ymin = -Inf, ymax = Inf, fill = col),
                  alpha = 0.4)+
        scale_fill_manual(values = "green4",guide = "none")+
        geom_density(data = df, 
                     aes(x = len, color = Sample),
                     size = .6)+
        theme_bw(base_size = 15)+
        scale_color_manual("Fragment", values = c("#6a00fc","#ffa10c"))+
        labs(title = "", y = "Density", x = "Fragment length (bp)")+
        scale_x_continuous(breaks = seq(50,155,20), limits = c(50,155))+
        th+
        facet_wrap( ~ patient_ID,scales = "free_y")+
        theme(strip.background = element_rect(fill="white"))+
        theme(strip.text = element_text(colour = 'black',face = "bold"))
    return(gg)
}


#x #Table consisting of BAM file names (name used when bamfile() is used),
#genes mutated and the position of the mutation
#y #Named list of input samples used
#z Named list of cfChIP samples used
#g Granges object of AVENIO genes

length_of_mutated_genes_df <- function(x,y,z,g){
    library(Rsamtools)
    library(tidyr)
    library(GenomicAlignments)
    or_wd <- getwd()
    `%ni%` <- Negate(`%in%`)
    df <- read.table(x,header = T)
    df$cfChIP <- names(z)
    df <- df[!is.na(df$genes),]
    df2 <- data.frame(len = c(),
                      Sample = c(),
                      group = c())
    n_samples <- length(df$file)
    for (i in 1:n_samples){
        nam <- unlist(lapply(strsplit(df$name[i],"_"),"[[",1))
        genes <- strsplit(df$genes[i],split = ",")[[1]]
        genes_grang <- g[mcols(g)$SYMBOL %in% genes]
        reads_input <- y[df$name[i]][[1]]
        reads_cfChIP <- z[df$cfChIP[i]][[1]]
        sub_input <- subsetByOverlaps(reads_input,genes_grang,ignore.strand = T)
        len_input <- mcols(sub_input)$isize
        len_input <- len_input[len_input >0]
        len_input <- len_input[!is.na(len_input)]
        sub_cfChIP <- subsetByOverlaps(reads_cfChIP,genes_grang,ignore.strand = T)
        len_cfChIP <- mcols(sub_cfChIP)$isize
        len_cfChIP <- len_cfChIP[len_cfChIP >0]
        len_cfChIP <- len_cfChIP[!is.na(len_cfChIP)]
        
        df1 <- data.frame(len = c(len_input,len_cfChIP),
                          Sample = c(rep("Input",length(len_input)),rep("cfChIP",length(len_cfChIP))),
                          group = rep(nam,length(len_input)+length(len_cfChIP)))
        df2 <- rbind(df2,df1)
    }
    return(df2)
}

#x data.frame of fragment lengths in mutated genes for input af cfChIP samples
#y if y = TRUE concatenated plot is plotted. If = "individual" individual plots per patients is made
length_of_mutated_genes_plot <- function(x,y){
    df <- x
    rects <- data.frame(xstart = 0, xend = 150, col = "col")
    if(isTRUE(y)){
        gg <- ggplot()+
            geom_rect(data = rects, aes(xmin = xstart, xmax = xend, 
                                        ymin = -Inf, ymax = Inf, fill = col),
                      alpha = 0.4)+
            scale_fill_manual(values = "green4",guide = "none")+
            geom_density(data = df, 
                         aes(x = len, color = Sample),
                         size = 1)+
            theme_bw(base_size = 15)+
            scale_color_manual(values = c("#6a00fc","#ffa10c"))+
            labs(title = "", y = "Density", x = "Fragment length (bp)")+
            scale_x_continuous(breaks = seq(0,400,50), limits = c(0,400))+
            th
    }
    if(y == "individual"){
        gg <- ggplot()+
            geom_rect(data = rects, aes(xmin = xstart, xmax = xend, 
                                        ymin = -Inf, ymax = Inf, fill = col),
                      alpha = 0.4)+
            scale_fill_manual(values = "green4",guide = "none")+
            geom_density(data = df, 
                         aes(x = len, color = Sample),
                         size = 1)+
            theme_bw(base_size = 15)+
            scale_color_manual(values = c("#6a00fc","#ffa10c"))+
            labs(title = "", y = "Density", x = "Fragment length (bp)")+
            scale_x_continuous(breaks = seq(0,400,50), limits = c(0,400))+
            th+
            facet_wrap( ~ group,scales = "free_y")+
            theme(strip.background =element_rect(fill="white"))+
            theme(strip.text = element_text(colour = 'black',face = "bold"))
    }
    return(gg)
}

length_of_mutated_genes_boxplot <- function(x){
    library(ggpubr)
    cfChIP_res <- x %>% filter(Sample == "cfChIP")
    input_res <- x %>% filter(Sample == "Input")
    df <- data.frame(proportions = c())
    for (i in 1:length(unique(x$group))){
        sample <- unique(x$group)[i]
        cfChIP_res_sample <- cfChIP_res %>%  filter(group == sample)
        input_res_sample <- input_res %>% filter(group == sample)
        cfChIP_sub <- cfChIP_res_sample %>% filter(len < 150)
        input_sub <- input_res_sample %>% filter(len < 150)
        cfChIP_prop <- nrow(cfChIP_sub)/nrow(cfChIP_res_sample)
        input_prop <- nrow(input_sub)/nrow(input_res_sample)
        df1 <- data.frame(proportions = c(cfChIP_prop,
                                          input_prop),
                          sample = factor(c("cfChIP",
                                           "Input"),
                                         levels = c("Input",
                                                    "cfChIP")),
                          relative_prop = c(cfChIP_prop/input_prop,
                                            input_prop/input_prop),
                          pairs = rep(i,2))
        df <- rbind(df,df1)
    }
    gg <- ggplot(df, aes(x = sample, y = proportions,
                         fill = sample))+
        geom_boxplot(outlier.alpha = 0)+
        geom_line(aes(group=pairs))+
        geom_point(alpha = 0.5)+
        labs(title = "",
             x = "",
             y = "Fraction under 150 bp")+
        scale_fill_manual("Fragment motif",
                          values = c("#ffa10c","#6a00fc"))+
        stat_compare_means(method = "t.test",
                           paired = TRUE,
                           comparisons = list(c("cfChIP", "Input")))+
        theme_bw()+
        th+
        theme(legend.position = "none")
    
    return(gg)
}

#x enrichment data.frame of all samples
#y pippin enrichment data.frame of all samples
pippin_vs_cfChIP <- function(x,y){
    library(dplyr)
    x <- x %>% dplyr::select(colnames(x)[colnames(x) %in% colnames(y)])
    df <- data.frame(quartile = c(),
                     pippin_enrichment = c(),
                     sample = c())
    for(j in 2:length(x)){
        x_sample <- x %>% dplyr::select(genes,j)
        x_sample <- x_sample[order(x_sample[,2]),]
        y_sample <- y %>% dplyr::select(genes,j)
        number_high <- 19
        number_low <- 1
        res <- c()
        for (i in 1:10){
            if(i < 5){
                gene <- x_sample$genes[number_low:number_high]
                x_sample_genes <- x_sample %>% filter(genes %in% gene)
                number_low <- number_low + 19
                number_high <- number_high + 19
                y_sample_genes <- y_sample %>% filter(genes %in% gene)
                res[i] <- mean(y_sample_genes[,2])
                if(i == 4){
                    number_high <- number_high + 1
                }
            }
            else{
                gene <- x_sample$genes[number_low:number_high]
                number_low <- number_low + 20
                number_high <- number_high + 20
                y_sample_genes <- y_sample %>% filter(genes %in% gene)
                res[i] <- mean(y_sample_genes[,2])
            }
        }
        df1 <- data.frame(quartile = c("Q1", "Q2", "Q3", "Q4",
                                       "Q5", "Q6", "Q7", "Q8",
                                       "Q9", "Q10"),
                          pippin_enrichment = res,
                          sample = rep(colnames(x_sample)[2],10))
        df <- rbind(df,df1)
    }
    return(df)
}

pippin_vs_cfChIP_plot <- function(x){
    library(ggpubr)
    library(rstatix)
    normalized <- c()
    for (i in 1:length(unique(x$sample))){
        sample_ID <- unique(x$sample)[i]
        x_sample <- x %>% filter(sample == sample_ID)
        res <- x_sample$pippin_enrichment/x_sample$pippin_enrichment[1]
        normalized <- c(normalized,res)
    }
    df <- x %>% mutate(normalized = normalized) %>% 
        mutate(quartile = factor(quartile,levels = c("Q1", "Q2", "Q3", "Q4",
                                                     "Q5", "Q6", "Q7", "Q8",
                                                     "Q9", "Q10")))
    ddf <- df %>% filter(quartile == "Q10")
    res.aov <- df %>% friedman_test(normalized ~ quartile | sample)
    gg <- ggplot(df, aes(x = quartile, y = normalized))+
        geom_boxplot(outlier.alpha = 0,
                     fill = "#6a00fc")+
        geom_jitter(alpha = 0.5, width = 0.1)+
        theme_bw()+
        stat_compare_means(aes(label = after_stat(p.signif)),
                           method = "t.test",
                           paired = TRUE,
                           ref.group = "Q1",
                           size = 4.5)+
        geom_text(label = get_test_label(res.aov, type = "text"),
                  x = 8.2, y = 1,size = 4.5)+
        labs(title ="",
             x = "Quantiles",
             y = "Normalized size selection enrichment")+
        th
    return(gg)
}


#x enrichment data.frame of all samples
#y pippin enrichment data.frame of all samples
pippin_vs_cfChIP_correlation <- function(x,y){
    df <- data.frame(genes = c(),
                     cfChIP = c(),
                     pippin = c(),
                     sample = c())
    x <- x %>% dplyr::select(colnames(x)[colnames(x) %in% colnames(y)])
    valuess <- c()
    for(j in 2:length(x)){
        x_sample <- x %>% dplyr::select(genes,j)
        y_sample <- y %>% dplyr::select(genes,j)
        joint <- x_sample %>% left_join(y_sample,by = "genes")
        df1 <- data.frame(genes = joint$genes,
                          cfChIP = joint[,2],
                          pippin = joint[,3])
        cor_res <- cor.test(df1$pippin, df1$cfChIP, method = "spearman")
        rho <- cor_res$estimate
        valuess[j] <- cor_res$p.value
        sample <- rep(colnames(x_sample)[2],nrow(joint))
        sample <- paste0(sample,"\n","Spearman r = ",round(as.numeric(rho),3)) 
        df1$sample <- sample
        df <- rbind(df,df1)
    }
    gg <- ggplot(df,aes(x = pippin, y = cfChIP))+
        geom_point(size = 2,color = "#6a00fc")+
        geom_smooth(method = "lm", se = F, color = "black")+
        theme_bw()+
        labs(title = "",x = "Size selection enrichment", y = "cfChIP enrichment")+
        th+
        facet_wrap(~sample, nrow = 3,scales = "free")+
        theme(strip.background =element_rect(fill="white"))+
        theme(strip.text = element_text(colour = 'black',face = "bold"))
    
    return(gg)
}
pippin_vs_cfChIP_correlation(enrichment_df,pippin_enrichment_df)

pippin_vs_cfChIP_correlation_matrix <- function(x){
    df <- x %>% mutate(place = rep(1,nrow(x))) %>% 
        mutate(sample = factor(sample,levels=c("NAC.1","NAC.2","NAC.3","NAC.4",
                                               "NSC.1","NSC.2","NSC.3",
                                               "SSC.1","SSC.2","SSC.3","SSC.4")))
    gg <- ggplot(df, aes(x = place,y = sample,fill = rho,
                         label = format(rho, nsmall = 3)))+
        geom_tile(color = "white")+
        scale_fill_gradient("Rho",low = "white", high = "#6a00fc",
                            limits = c(0, 1))+
        theme_minimal()+
        labs(x = "",y="")+
        geom_text(col = "black",size = 3)+
        coord_equal()+
        theme(axis.text.x = element_blank(),
              legend.title = element_text(size = 12,face = "bold"),
              axis.text.y = element_text(face = "bold"),
              axis.line.x = element_blank(),
              axis.ticks = element_blank(),
              axis.line = element_blank())
    return(gg)
}


#x name of file with mutated genes, MAF in input and pippin as well as sample
pippin_effect_MAF <- function(x){
    library(dplyr)
    library(ggplot2)
    library(ggpubr)
    df <- read.table(x, header = T)
    ddf <- data.frame(gene = c(),
                      MAF = c(),
                      type = c(),
                      sample = c())
    for(i in 1:length(unique(df$sample))){
        nam <- unique(df$sample)[i]
        df_s <- df %>% filter(sample == nam)
        df1 <- data.frame(gene = rep(df_s$gene,2),
                          MAF = c(df_s$MAF_input,df_s$MAF_pippin),
                          type = c(rep("Input",nrow(df_s)),
                                   rep("Pippin",nrow(df_s))),
                          sample = rep(nam,nrow(df_s)*2))
        ddf <- rbind(ddf,df1)
    }
    ddf <- ddf %>% mutate(type = factor(type, levels = c("Input","Pippin")))
    ddf_input <- ddf %>% filter(type == "Input") %>% mutate(pairs = 1:(nrow(df))/2)
    ddf_pippin <- ddf %>% filter(type == "Pippin") %>% mutate(pairs = 1:(nrow(df))/2) %>% 
        mutate(type = rep("Size selected",nrow(ddf_input)))
    ddf <- rbind(ddf_input,ddf_pippin)
    gg <- ggplot(ddf, aes(x = type, y = MAF,
                         fill = type))+
        geom_boxplot(outlier.alpha = 0)+
        geom_line(aes(group=pairs))+
        geom_point(alpha = 0.5)+
        labs(title = "",
             x = "",
             y = "Molecular allele fraction (%)")+
        scale_fill_manual("Sample",
                          values = c("#ffa10c","#6a00fc"))+
        stat_compare_means(method = "t.test",
                           paired = TRUE,
                           comparisons = list(c("Input", "Size selected")))+
        theme_bw()+
        th+
        theme(legend.position = "none")
    ddf_input <- ddf %>% filter(type == "Input") %>% mutate(pairs = 1:(nrow(df))/2) %>% 
        summarise(mean = mean(MAF),
                  CI = confint(lm(MAF ~ 1),level = 0.95))
    print(paste("Input:",ddf_input))
    ddf_pippin <- ddf %>% filter(type == "Size selected") %>% mutate(pairs = 1:(nrow(df))/2) %>% 
        summarise(mean = mean(MAF),
                  CI = confint(lm(MAF ~ 1),level = 0.95))
    print(paste("Pippin:",ddf_pippin))
    return(gg)
}

#x input BAM file returned by bamfile()
#y lower limit for fragment size selection
#z upper limit for fragment size selection
in_silico_size_selection <- function(x,y,z){
    bins <- seq(y,z,by = 10)
    mcols(x)$size <- abs(mcols(x)$isize)
    x <- x[!is.na(mcols(x)$size)]
    x <- x[mcols(x)$size != 0]
    res_list <- list()
    for(i in 1:length(bins)){
        short_reads <- x[mcols(x)$size<bins[i]]
        print(bins[i])
        res_list <- append(res_list,list(short_reads))
    }
    names(res_list) <- as.character(bins)
    return(res_list)
}
#x Gene count table returned by gene_count()
#y .txt file of the AVENIO genes with the number of bases sequenced in each gene. based on data.frame returned by coverages()
#z Granges object returned by gr() for the 197 AVENIO target genes
e.score_mod <- function(x,y,z){
    library(GenomicAlignments)
    library(GenomicRanges)
    library(BiocParallel)
    library(dplyr)
    y <- read.table(y, header=T)
    y$coverage <- as.numeric(y$coverage)
    ChIP_reads <- x
    ChIP_reads <- ChIP_reads %>% filter(readcounts > 1)
    y <- y[match(ChIP_reads$genes, y$SYMBOL),]
    len <- sum(ChIP_reads$readcounts)
    RPKM <- c()
    for (i in 1:length(ChIP_reads$genes)){
        RPKM[i] <- (ChIP_reads$readcounts[i]*1000*1000000)/(y$coverage[i])
    }
    ChIP_reads$e <- RPKM
    res <- data.frame(genes = ChIP_reads$genes, enrichment=ChIP_reads$e)
    res <- res[order(res$enrichment),]
    return(res)
}

#x list of Granges returned by in_silico_size_selection()
#y Name of sample
#z gene count from input returned by gene_count()
#g Granges object returned by gr()

in_silico_enrichment <- function(x,y,z,g){
    df <- data.frame(gene = mcols(g)$SYMBOL)
    for(i in 1:length(x)){
        reads <- x[[i]]
        reads_counts <- c()
        overlap_res <- findOverlaps(reads,g)
        df1 <- data.frame(genes = mcols(g)$SYMBOL[subjectHits(overlap_res)]) %>% 
            dplyr::count(genes) %>% 
            left_join(z,by = "genes") %>% 
            mutate(fraction = n/(readcounts)) %>% 
            dplyr::select(genes,fraction) %>% 
            arrange(genes)
        print(names(x)[i])
        colnames(df1) <- c("gene", paste0(y,"_",names(x)[i],"bp"))
        df <- left_join(df,df1,by = "gene")
    }
    return(df)
}
#in_silico_enrichment(NAC.1_in_silico,"NAC.1",NAC.1_input_count,grs)

#x collected in silico size selction data.frame of all patients
#y collected cfChIP eller pippin enrichment data.frame
#z the sample number to plot
in_silico_correlation_plot <- function(x,y,z){
    library(tidyr)
    y <- y[!grepl("HC",colnames(y))]
    for(i in 2:length(y)){
        samp <- colnames(y[i])
        enrichment_sample <- y[c(1,i)]
        colnames(enrichment_sample) <- c("gene",samp)
        x_sample <- x %>% dplyr::select(c("gene",colnames(x)[grepl(samp,colnames(x))]))
        x_sample <- na.omit(x_sample)
        df <- x_sample %>% left_join(enrichment_sample, by = "gene") 
        colnames(df)[length(df)] <- "control"
        df <- df %>% pivot_longer(colnames(df)[grepl(samp,colnames(df))])
        ddf <- data.frame(gene = c(),
                          control = c(),
                          name = c(),
                          value = c()) 
        for(j in 1:length(unique(df$name))){
            df1 <- df %>% filter(name == unique(df$name)[j])
            cor_res <- cor.test(df1$control, df1$value, method = "spearman")
            rho <- cor_res$estimate
            samples <- paste0(df1$name,"\n","Spearman r = ",round(as.numeric(rho),3))
            df1$name <- samples
            ddf <- rbind(ddf,df1)
        }
        ddf <- ddf %>% mutate(name = factor(name, levels = unique(ddf$name)))
        gg <- ggplot(ddf, aes(x = value, y = control))+
            geom_point(size = 2,color = "#6a00fc")+
            geom_smooth(method = "lm", se = F, color = "black")+
            theme_bw()+
            labs(title = "",x = "In silico size selection fraction", y = "cfChIP enrichment")+
            th+
            facet_wrap(~name, nrow = 3,scales = "free")+
            theme(strip.background =element_rect(fill="white"))+
            theme(strip.text = element_text(colour = 'black',face = "bold"))
        if(i == z){
            return(gg)
        }
    }
    
}
#in_silico_correlation_plot(collected_in_silico_fraction,enrichment_df,6)
#x collected in silico size selction data.frame of all patients
#y collected cfChIP eller pippin enrichment data.frame
#z character either "cfChIP" eller "pippin" indicating what type of data to be compared
in_silico_matrix_cor_plot <- function(x,y,z){
    library(tidyr)
    library(ggpubr)
    y <- y[!grepl("HC",colnames(y))]
    ddf <- data.frame(sample = c(),
                      rho = c(),
                      cutoff = c())
    for(i in 2:length(y)){
        samp <- colnames(y[i])
        enrichment_sample <- y[c(1,i)]
        colnames(enrichment_sample) <- c("gene",samp)
        x_sample <- x %>% dplyr::select(c("gene",colnames(x)[grepl(samp,colnames(x))]))
        x_sample <- na.omit(x_sample)
        df <- x_sample %>% left_join(enrichment_sample, by = "gene") 
        colnames(df)[length(df)] <- "control"
        df <- df %>% pivot_longer(colnames(df)[grepl(samp,colnames(df))])
        rhos <- c()
        bp <- c()
        for(j in 1:length(unique(df$name))){
            df1 <- df %>% filter(name == unique(df$name)[j])
            cor_res <- cor.test(df1$control, df1$value, method = "spearman")
            rho <- cor_res$estimate
            rhos[j] <- round(as.numeric(rho),3)
            bp[j] <- strsplit(df1$name,"_")[[1]][2]
        }
        ddf1 <- data.frame(sample = rep(samp,length(rhos)),
                          rho = rhos,
                          cutoff = bp)
        ddf <- rbind(ddf,ddf1)
    }    
    ddf <- ddf %>% mutate(cutoff = factor(cutoff, levels = unique(ddf$cutoff)))
    df2 <- data.frame(sample = c(rep("",length(unique(ddf$cutoff))),
                                 rep("Relative Rho",length(unique(ddf$cutoff)))),
                      rho = rep(0,length(unique(ddf$cutoff)*2)),
                      cutoff = rep(unique(ddf$cutoff),2))
    ddf_group <- ddf %>% group_by(sample) %>% 
        summarise(.groups = "keep",
                  avg = median(rho),
                  relative = rho/avg,
                  cutoff = cutoff) %>% 
        group_by(cutoff) %>% 
        summarise(.groups = "keep",
                  median_relative = median(relative))
    ddf <- rbind(ddf,df2)
    ddf <- ddf %>% mutate(sample = factor(sample, levels = rev(unique(ddf$sample))))
    gg <- ggplot()+
        geom_tile(data = ddf, aes(x = cutoff,y = sample,fill = rho),
                  color = "white")+
        scale_fill_gradient("Rho",low = "white", high = "#6a00fc",
                            limits = c(0, 1))+
        theme_minimal()+
        labs(x = "",y="")+
        geom_text(data = ddf, 
                  aes(x = cutoff,y = sample,
                  label = format(rho, nsmall = 3),
                  color = ifelse(rho == 0, "no", "yes")),
                  size = 3)+
        scale_color_manual(values=c("yes"="black", "no"="white"))+
        coord_equal()+
        theme(legend.title = element_text(size = 12,face = "bold"),
              axis.text.y = element_text(face = "bold"),
              axis.line.x = element_blank(),
              axis.ticks.y.left = element_blank(),
              axis.line.y.left = element_blank())+
        geom_bar(data = ddf_group, aes(x = cutoff, y = median_relative),
                 stat = "identity", fill = "#ffa10c")+
        guides(color = "none")
    
    return(gg)
}


#x #Name of BAM file
#y #Granges object returned by gr() for the 197 AVENIO target genes
gene_count <- function(x,y){
    library(Rsamtools)
    bf <- BamFile(x)
    counts <- scanBam(bf, param = ScanBamParam(what = "pos", which = y))
    reads <- c()
    name <- c()
    for (i in 1:length(counts)){
        reads[i] <- length(counts[[i]]$pos)
        string <- names(counts)[i]
        splits <- strsplit(string, ":")
        chr <- splits[[1]][1]
        position <- as.numeric(strsplit(splits[[1]][2],"-")[[1]][1])
        g <- GRanges(chr, IRanges(start = position+10, end = position+110))
        gene <- y[y %over% g]
        name[i] <- gene$SYMBOL
    }
    df <- data.frame(genes = name, readcounts = reads)
    df <- df[order(df$genes),]
    return(df)
}
#N The number of possible motifs. default = 64
#A The number of motifs not overlapping for group1
#B The number of motifs not overlapping for group2
#k the numer of overlapping motifs
enrich_pvalue <- function(N = 64, A, B, k){
    library(gmp)
    m <- A + k
    n <- B + k
    i <- k:min(m,n)
    res <- as.numeric(sum(chooseZ(m,i)*chooseZ(N-m,n-i))/chooseZ(N,n))
    return(res)
}
enrich_pvalue(40352, 5049-261,329-261,261)

#x Granges returned by gr()
#y names list of transcript ID's, where names are the gene names
exon_extract <- function(x,y){
    library(GenomicFeatures)
    library(TxDb.Hsapiens.UCSC.hg38.knownGene)
    genome <- TxDb.Hsapiens.UCSC.hg38.knownGene
    coding_gene_names <- c("BRCA1", "BRCA2", 
                           "EGFR", "ERBB2",
                           "KRAS", "MET")
    coding_gene_transcript_ID <- y
    coding_genes <- x[mcols(x)$SYMBOL %in% coding_gene_names]
    all_exons <- exonsBy(genome,by = "tx",use.names = T)
    coding_exons <- suppressWarnings(
        unlist(
            all_exons[names(all_exons) %in% coding_gene_transcript_ID]))
    sub <- findOverlaps(coding_exons, coding_genes)
    gene_names <- mcols(coding_genes)$SYMBOL[subjectHits(sub)]
    coding_exons <- coding_exons[queryHits(sub)]
    mcols(coding_exons)$SYMBOL <- gene_names
    coding_exons <- coding_exons[order(mcols(coding_exons)$SYMBOL)]
    return(coding_exons)

}

#x BAM file of input reads returned by bamfile()
#y Granges object of genes returned by gr()
#z exons to be considered returned by exon_extract()
sub150_in_exons <- function(x,y,z,name){
    library(GenomicRanges)
    genes <- y[mcols(y)$SYMBOL %in% mcols(z)$SYMBOL]
    bam_overlap <- subsetByOverlaps(x,genes,ignore.strand = T)
    bam_overlap <- bam_overlap[!is.na(mcols(bam_overlap)$isize)]
    mcols(bam_overlap)$isize <- abs(mcols(bam_overlap)$isize)
    bam_overlap <- bam_overlap[mcols(bam_overlap)$isize < 5000]
    bam_overlap <- bam_overlap[!mcols(bam_overlap)$isize == 0]
    fragments <- GRanges(seqnames = seqnames(bam_overlap),
                         ranges = IRanges(start = start(bam_overlap),
                                          end = mcols(bam_overlap)$isize+start(bam_overlap)),
                         strand = strand(bam_overlap))
    mcols(fragments)$isize <- mcols(bam_overlap)$isize
    overlapping <- findOverlaps(fragments,genes)
    gene_names <- mcols(genes)$SYMBOL[subjectHits(overlapping)]
    fragments <- fragments[queryHits(overlapping)]
    mcols(fragments)$SYMBOL <- gene_names
    df <- data.frame(gene = c(),
                      exon = c(),
                      fraction = c(),
                      norm_fraction = c())
    for (i in 1:length(mcols(genes)$SYMBOL)){
        gen <- mcols(genes)$SYMBOL[i]
        fragment_gene <- fragments[mcols(fragments)$SYMBOL == gen]
        gene_sub150 <- mean(mcols(fragment_gene)$isize < 150)
        exons_gene <- z[mcols(z)$SYMBOL == gen]
        exon_sub <- findOverlaps(fragment_gene,exons_gene, 
                     type = "any", ignore.strand = T)
        fragment_gene <- fragment_gene[queryHits(exon_sub)]
        mcols(fragment_gene)$exon <- subjectHits(exon_sub)
        options("showHeadLines"=50)
        exon = c()
        gene = c()
        fraction = c()
        global_fraction = c()
        norm_fraction = c()
        n_fragments = c()
        for (j in 1:length(unique(mcols(fragment_gene)$exon))){
            ex <- unique(mcols(fragment_gene)$exon)[j]
            frag_exon <- fragment_gene[mcols(fragment_gene)$exon == ex]
            exon[j] <- ex
            gene[j] <- gen
            if(length(mcols(frag_exon)$isize) < 200){
                fraction[j] <- NA
                norm_fraction[j] <- NA
                n_fragments[j] <- NA
            }
            else{
                fraction[j] <- mean(mcols(frag_exon)$isize < 150)
                norm_fraction[j] <- mean(mcols(frag_exon)$isize < 150)/gene_sub150
                n_fragments[j] <- length(mcols(frag_exon)$isize)
            }
            global_fraction[j] <- gene_sub150
        }
        df1 <- data.frame(gene = gene,
                          exon = exon,
                          fraction = fraction,
                          norm_fraction = norm_fraction,
                          global_fraction = global_fraction,
                          n_fragments = n_fragments)
        df1 <- na.omit(df1)
        df <- rbind(df,df1)
    }
    df$sample <- rep(name, nrow(df))
    return(df)
}

#x Collected data.frame of normalized sub150 fraction in each exon
# Enrichment data.frame of cfChIP data
sub150_exon_plot <- function(x,y){
    library(ggplot2)
    library(dplyr)
    y <- y %>% filter(genes %in% x$gene) %>% dplyr::select(genes,
                                                    unique(x$sample))
    y_red <- y %>% dplyr::select(-genes) 
    y_red <- y_red %>% rowwise() %>% 
        mutate(med = median(c(NAC.1, NAC.2, NAC.3, NAC.4,
                              NSC.1, NSC.2, NSC.3, NSC.4,
                              SSC.1, SSC.2, SSC.3, SSC.4))) %>%
        ungroup() %>% 
        mutate(NAC.1_norm = NAC.1/med,
               NAC.2_norm = NAC.2/med,
               NAC.3_norm = NAC.3/med,
               NAC.4_norm = NAC.4/med,
               NSC.1_norm = NSC.1/med,
               NSC.2_norm = NSC.2/med,
               NSC.3_norm = NSC.3/med,
               NSC.4_norm = NSC.4/med,
               SSC.1_norm = SSC.1/med,
               SSC.2_norm = SSC.2/med,
               SSC.3_norm = SSC.3/med,
               SSC.4_norm = SSC.4/med) %>% 
        mutate(genes = y$genes)
    all_values <- c()
    for (i in 1:length(unique(x$sample))){
        samp <- x %>% filter(sample == unique(x$sample)[i])
        for (j in 1:length(unique(samp$gene))){
            gene_samp <- samp %>% filter(gene == unique(samp$gene)[j])
            ssamp <- paste0(gene_samp$sample[1],"_norm")
            y_red_samp <- y_red %>% dplyr::select(genes,ssamp) %>% 
                filter(genes == unique(samp$gene)[j])
            gene_res <- as.numeric(y_red_samp[1,2])
            new_values <- rep(gene_res,nrow(gene_samp))
            all_values <- c(all_values, new_values)
        }
    }
    x$expression <- log2(all_values)
    gg <- ggplot(x, aes(x = exon, y = norm_fraction, group = expression,
                  color = expression))+
        geom_line(linewidth = 1)+
        theme_bw()+
        th+
        scale_color_gradient("Log2 relative expression",
                             low = "#ffa10c",
                             high = "#6a00fc")+
        facet_wrap(~gene, nrow = 2,scales = "free")
        
   return(gg) 
}

fragment_endpoint_wt_ctdna_df <- function(x,y){
    library(Rsamtools)
    library(tidyr)
    library(GenomicAlignments)
    or_wd <- getwd()
    `%ni%` <- Negate(`%in%`)
    df <- read.table(x,header = T)
    df <- df[!is.na(df$genes),]
    df6 <- data.frame(len = c(),
                      Sample = c(),
                      group = c())
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
            print(paste("j =",j))
            mut_pos <- pos[j]
            mut_chr <- chr[j]
            g <- GRanges(mut_chr,IRanges(start = mut_pos, end = mut_pos),
                         strand = "*")
            sub <- subsetByOverlaps(reads,g,ignore.strand = T)
            type <- types[j]
            len <- mcols(sub)$isize
            cig <- mcols(sub)$cigar
            cig <- cig[!is.na(len)]
            sub <- sub[!is.na(len)]
            len <- len[!is.na(len)]
            cig <- cig[len !=0]
            sub <- sub[len !=0]
            len <- len[len !=0]
            len <- abs(len)
            df3 <- data.frame(len = len,
                       pos_start = start(sub),
                       pos_end = end(sub),
                       strand = as.character(strand(sub)),
                       mut_pos = rep(mut_pos,length(len)))
            df3 <- df3 %>% mutate(start_rel = ifelse(strand == "+",
                                                     pos_start-mut_pos,
                                                     pos_end-len-mut_pos)) %>% 
                mutate(end_rel = ifelse(strand == "+",
                                        pos_start+len-mut_pos,
                                        pos_end-mut_pos))
            if(grepl("indel",type)){
                indel <- strsplit(type,split = "-")[[1]][2]
                ct_fragment_start_rel <- df3$start_rel[grepl(indel,cig)]
                ct_fragment_end_rel <- df3$end_rel[grepl(indel,cig)]
                wt_fragment_start_rel <- df3$start_rel[!grepl(indel,cig)]
                wt_fragment_end_rel <- df3$end_rel[!grepl(indel,cig)]
                df1 <- data.frame(relative = c(ct_fragment_start_rel,
                                           ct_fragment_end_rel,
                                           wt_fragment_start_rel,
                                           wt_fragment_end_rel),
                                  Sample = c(rep("ctDNA",length(ct_fragment_start_rel)+
                                                     length(ct_fragment_end_rel)),
                                             rep("WT", length(wt_fragment_start_rel)+
                                                     length(wt_fragment_end_rel))),
                                  ends = c(rep("start",length(ct_fragment_start_rel)),
                                           rep("end",length(ct_fragment_start_rel)),
                                           rep("start",length(wt_fragment_start_rel)),
                                           rep("end",length(wt_fragment_start_rel))),
                                  group = rep(df$name[i],length(length(ct_fragment_start_rel)+
                                                                    length(ct_fragment_end_rel)+
                                                                    length(wt_fragment_start_rel)+
                                                      length(wt_fragment_end_rel))))
                df6 <- rbind(df6,df1)
            }
            else{
                param = ScanBamParam(which = g, what = c("strand", "qwidth","isize","pos", "cigar"),
                                     tag = "MD")
                setwd("D:/Lung cancer input/PosDeduped")
                bam <- scanBam(df$file[i],
                               param = param)
                setwd(or_wd)
                bam_pos <- as.numeric(unlist(bam[[1]][2]))
                bam_size <- as.numeric(unlist(bam[[1]][5]))
                bam_cigar <- unlist(bam[[1]][4])
                bam_MD <- unlist(bam[[1]][6])
                bam_width <- as.numeric(unlist(bam[[1]][3]))
                bam_strand <- unlist(bam[[1]][1])
                df3 <- data.frame(len = abs(bam_size),
                                  pos_start = bam_pos,
                                  pos_end = bam_pos+bam_width,
                                  strand = bam_strand,
                                  mut_pos = rep(mut_pos,length(bam_pos)),
                                  MD = bam_MD)
                df3 <- df3[!is.na(df3$pos_start),]
                df3 <- df3[df3$len !=0,]
                df3 <- df3 %>% mutate(start_rel = ifelse(strand == "+",
                                                         pos_start-mut_pos,
                                                         pos_end-len-mut_pos)) %>% 
                    mutate(end_rel = ifelse(strand == "+",
                                            pos_start+len-mut_pos,
                                            pos_end-mut_pos))
                base <- strsplit(type,split = "-")[[1]][2]
                ct_fragment_start_rel <- df3$start_rel[grepl(base,df3$MD)]
                ct_fragment_end_rel <- df3$end_rel[grepl(base,df3$MD)]
                wt_fragment_start_rel <- df3$start_rel[!grepl(base,df3$MD)]
                wt_fragment_end_rel <- df3$end_rel[!grepl(base,df3$MD)]
                df4 <- data.frame(relative = c(ct_fragment_start_rel,
                                               ct_fragment_end_rel,
                                               wt_fragment_start_rel,
                                               wt_fragment_end_rel),
                                  Sample = c(rep("ctDNA",length(ct_fragment_start_rel)+
                                                     length(ct_fragment_end_rel)),
                                             rep("WT", length(wt_fragment_start_rel)+
                                                     length(wt_fragment_end_rel))),
                                  ends = c(rep("start",length(ct_fragment_start_rel)),
                                           rep("end",length(ct_fragment_start_rel)),
                                           rep("start",length(wt_fragment_start_rel)),
                                           rep("end",length(wt_fragment_start_rel))),
                                  group = rep(df$name[i],length(length(ct_fragment_start_rel)+
                                                                    length(ct_fragment_end_rel)+
                                                                    length(wt_fragment_start_rel)+
                                                                    length(wt_fragment_end_rel))))
                df6 <- rbind(df6,df4)
                print(i)
            }
        }    
    }
    return(df6)
}


#x #data.frame returned by fragment_length_wt_ctdna_df()
#y If y = TRUE (default) concatenated lines are plotted. Otherwise the plot is per patient. 
fragment_endpoint_wt_ctdna_plot <- function(x,y = TRUE){
    library(ggplot2)
    library(tidyr)
    library(plyr)
    if(y){
        x$Sample <- factor(x$Sample, levels = c("ctDNA", "WT"))
        gg <- ggplot()+
            scale_fill_manual(values = "green4",guide = "none")+
            geom_density(data = x, 
                         aes(x = relative, color = Sample),
                         size = 1)+
            theme_bw(base_size = 15)+
            scale_color_manual("Fragment", values = c("#6a00fc","#ffa10c"))+
            labs(title = "", y = "Density", x = "Relative fragment endpoint (bp)")+
            scale_x_continuous(breaks = seq(-200,200,50), limits = c(-200,200))+
            th
        
    }
    else{
        x$Sample <- factor(x$Sample, levels = c("ctDNA", "WT"))
        x <- x %>% mutate(group = unlist(lapply(strsplit(x$group,"_"), "[[", 1)))
        gg <- ggplot()+
            scale_fill_manual(values = "green4",guide = "none")+
            geom_density(data = x, 
                         aes(x = relative, color = Sample),
                         size = 1)+
            theme_bw(base_size = 15)+
            scale_color_manual("Fragment", values = c("#6a00fc","#ffa10c"))+
            labs(title = "", y = "Density", x = "Relative fragment endpoint (bp)")+
            scale_x_continuous(breaks = seq(-200,200,50), limits = c(-200,200))+
            th+
            facet_wrap( ~ group,scales = "free_y")+
            theme(strip.background =element_rect(fill="white"))+
            theme(strip.text = element_text(colour = 'black',face = "bold"))
    }
    return(gg) 
    
}


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
prop_df <- function(x){
    df <- as.data.frame.table(proportions(x))
    return(df)
}

count_df <- function(x){
    df <- as.data.frame.table(x)
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
        theme_logo()+
        th+
        theme(axis.text.x = element_blank(),
              strip.text = element_text(face = "bold",
                                        size = 16))
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
    df <- x2 %>% left_join(y2)
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
matrix_plot(multiple_regions_correlations,"Correlation statistics for multiple region genes")
matrix_plot(all_genes_correlations,"Correlation statistics for all genes")
library(ggnewscale)
install.packages("ggnewscale")
ggnewscale::new_scale_fill()

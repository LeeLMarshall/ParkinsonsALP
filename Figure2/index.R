require(pacman)
p_load(
  "dplyr", "data.table",
  "foreach", "circlize",
  "ggplot2", "RobustRankAggreg", 
  "glue", "cowplot", "broom",
  "gridGraphics", "magick",
  "sitools", "ggrepel", "ggdendro")


theme_publication <- function(font = 7) {
  theme(
    text = element_text(family = "Helvetica", size = font, color = "black"),
    plot.title = element_text(size = font, color = "black"),
    
    panel.background = element_rect(fill = "transparent", colour = NA),
    panel.border = element_blank(),
    plot.background = element_rect(fill = "transparent", colour = NA),
    
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    
    axis.line = element_line(color = "black"),
    axis.text = element_text(size = font, color = "black"), 
    axis.ticks = element_line(size = 0.3, color = "black"),
    axis.title = element_text(size = font, color = "black"),
    
    # strip.background = element_rect(color = "white", fill = "white"),
    strip.text.x = element_text(size = font, color = "black"),
    strip.text.y = element_text(size = font, color = "black"),
    
    legend.background = element_rect(
      fill = alpha("white", 0),
      color = alpha("white", 0)
    ),
    legend.key = element_rect(color = NA, fill = NA),
    legend.key.size = unit(0.2, "cm"),    
    legend.text = element_text(size = font)
  )
}

# Load the fits
setwd("~/tmp/PD2018")
pfcr <- "Brain_PFCRep_Padlock_withGLU/www/m6*.csv" %>% Sys.glob %>% fread
pfc  <- "Brain_PFC_Padlock_CGonly/www/m6*.csv" %>% Sys.glob %>% fread
ofb  <- "Brain_OFB_Padlock_CGonly/www/m6*.csv" %>% Sys.glob %>% fread
app  <- "Appendix_PDvsControls_Padlock/www/m6_*.csv" %>% Sys.glob %>% fread


# Chromosomes we work with
chroms <- pfcr[, unique(Chr)] %>% 
  gsub("chr", "", .) %>% 
  as.numeric %>% 
  sort %>% 
  paste0("chr", .)

# Track of CG-wise logP values for each dataset
bedPFCR <- 
  pfcr[adj.P.Val < 0.05] %>% 
  .[, list(chr = Chr, start = SNP, end = SNP, 
           value = ifelse(P.Value < 1e-12, 12, -log10(P.Value)))] %>% 
  as.data.frame
bedPFC <- 
  pfc[adj.P.Val < 0.05] %>% 
  .[, list(chr = Chr, start = SNP, end = SNP, 
           value = ifelse(P.Value < 1e-12, 12, -log10(P.Value)))] %>% 
  as.data.frame
bedOFB <- 
  ofb[adj.P.Val < 0.05] %>% 
  .[, list(chr = Chr, start = SNP, end = SNP, 
           value = ifelse(P.Value < 1e-12, 12, -log10(P.Value)))] %>% 
  as.data.frame
bedAPP <- 
  app[adj.P.Val < 0.05] %>% 
  .[, list(chr = Chr, start = SNP, end = SNP, 
           value = ifelse(P.Value < 1e-12, 12, -log10(P.Value)))] %>% 
  as.data.frame

bed <- list(
  bedAPP, bedOFB, bedPFC, bedPFCR 
)


# # Track of best p value per gene
# pd <- 
#   merge(
#     pfcr[, list(Chr, SNP, ID, Gene, P.PFCR = -log10(adj.P.Val))], 
#     pfc[, list(ID, P.PFC = -log10(adj.P.Val))], 
#     by = "ID") %>% 
#   merge(., ofb[, list(ID, P.OFB = -log10(adj.P.Val))], 
#         by = "ID") %>% 
#   merge(., app[, list(ID, P.APP = -log10(adj.P.Val))],
#         by = "ID")
# 
# pd <- pd[, list(
#             Chr = unique(Chr), 
#             Start = min(SNP), 
#             End = max(SNP), 
#             P.PFCR = max(P.PFCR), 
#             P.PFC = max(P.PFC), 
#             P.OFB = max(P.OFB), 
#             P.APP = max(P.APP)
#             ), 
#          Gene]
# 
# bed2 <- pd[, list(Chr, Start, End,
#                   P.PFCR, P.PFC,
#                   P.OFB, P.APP)]
# 
# 
# # And now label the hallmark genes
# hallmarks <- aggregateRanks(
#   list(
#     pd[order(P.PFCR, decreasing=TRUE), Gene],
#     pd[order(P.APP, decreasing=TRUE), Gene],
#     pd[order(P.OFB, decreasing=TRUE), Gene]
#     # pd[order(P.PFC, decreasing=TRUE), Gene]
#   )
# )  %>% rownames
# pd[match(hallmarks, Gene), ] %>% fwrite("hallmarks.csv")
# hallmarks <- hallmarks[1:10]

# Track of gene enrichment with significant loci and hallmarks
# Enrichment of AP target genes with significant loci
getEnrichments <- function(data, name, pT = 0.2) {
  a <- data[, list(F = sum(adj.P.Val > pT, na.rm=TRUE), T = sum(adj.P.Val < pT, na.rm=TRUE))]
  b <- data[, list(F = sum(adj.P.Val > pT, na.rm=TRUE), T = sum(adj.P.Val < pT, na.rm=TRUE)), Gene]
  dt <- 
    foreach(gene = b$Gene, .combine = rbind) %do% {
      rbind(
        a - b[Gene == gene, 2:3, with=FALSE],
        b[Gene == gene, 2:3, with=FALSE]
      ) %>%
        fisher.test(alternative = "greater") %>%
        broom::tidy() %>%
        as.data.table %>% 
        .[, Gene := gene]
    }
  dt <- dt[, list(Gene, estimate, p.value)]
  dt
}

pT <- 0.05
epfcr <- getEnrichments(pfcr, "PFCR", pT)
epfc  <- getEnrichments(pfc,  "PFC",  pT)
eofb  <- getEnrichments(ofb,  "OFB",  pT)
eapp  <- getEnrichments(app,  "APP",  pT)

hallmarks <- aggregateRanks(
  list(
    epfcr[order(p.value), Gene],
    # epfc[order(p.value), Gene],
    eofb[order(p.value), Gene],
    eapp[order(p.value), Gene]
  )
)  %>% rownames
hallmarks <- hallmarks %>% head(20)

pd <- 
epfcr[, list(Gene, P.PFCR= -log10(p.value))] %>%
  merge(.,
    epfc[ , list(Gene, P.PFC = -log10(p.value))],
    by = "Gene"
  ) %>% 
  merge(., 
        eofb[, list(Gene, P.OFB = -log10(p.value))],
        by = "Gene") %>%
  merge(.,
        eapp[, list(Gene, P.APP = -log10(p.value))],
        by = "Gene")

pd <- merge(
  pfcr[, list(Chr = unique(Chr), Start = min(SNP), End = max(SNP)), Gene],
  pd, by = "Gene")


bed2 <- pd[, list(Chr, Start, End,
                  P.PFCR, P.PFC,
                  P.OFB, P.APP)]

# # Get also the PD GWAS genes
# geneList <- c("GBA","NUCKS1","SLC41A1","SIPA1L2","TMEM163","CCNT2","STK39",
#               "CHMP2B","MCCC1","TMEM175","DGKQ","FAM200B","CD38","FAM47E",
#               "SNCA","HLA-DRB6","HLA-DQA1","KLHL7","NUPL2","GPNMB","MICU3",
#               "BAG3","DLG2","MIR4697","LRRK2","OGFOD2","GCH1","TMEM229B","VPS13C",
#               "ZNF646","KAT8","ARHGAP27","CRHR1","SPPL2C","MAPT","STH","KANSL1",
#               "SYT4","LSM7","DDRGK1","ITPKB","IL1R2","SCN3A","SATB1","NCKIPSD","CDC71",
#               "ALAS1","TLR9","DNAH1","BAP1","PHF7","NISCH","STAB1","ITIH3","ITIH4","ANK2",
#               "CAMK2D","ELOVL7","ZNF184","CTSB","SORBS3","PDLIM2","C8orf58","BIN3","SH3GL2",
#               "FAM171A1","GALC","COQ7","TOX3","ATP6V0A1","PSMC3IP","TUBG2")
# # Reduce them to genes that have enrichment in any of the tracks
# geneList <- pd[Gene %in% geneList & (P.PFCR > 1.3 | P.PFC > 1.3 |
#                                     P.OFB > 1.3 | P.APP > 1.3), Gene]

# Get regions for hallmark genes
# bed3 <- 
#   pfcr[Gene %in% setdiff(geneList, hallmarks), list(unique(Chr), min(SNP), max(SNP)), Gene] %>% 
#     .[, list(Chr = V1, Start = V2, End = V3, Label = Gene, Type = '#7570b3')] %>%
#     .[order(Chr, Start)]
bed4 <-
  pfcr[Gene %in% hallmarks, list(unique(Chr), min(SNP), max(SNP)), Gene] %>% 
    .[, list(Chr = V1, Start = V2, End = V3, Label = Gene, Type = '#d95f02')] %>%
    .[order(Chr, Start)]
# bed5 <- 
#   pfcr[Gene == "SNCA", list(unique(Chr), min(SNP), max(SNP)), Gene] %>% 
#   .[, list(Chr = V1, Start = V2, End = V3, Label = Gene, Type = 'black')] %>%
#   .[order(Chr, Start)]

# bed3 <- rbind(bed3, bed4, bed5) %>% as.data.frame
# rm(bed4, bed5)
bed3 <- bed4
rm(bed4)


#### Tracks ready, do the plot ####
pA <- function() {
  colors <- c("#7570b3", "#e7298a", "#66a61e", "#e6ab02")
  circos.par("start.degree" = 90, gap.after = c(rep(1, length(chroms)-1), 15)) 
  circos.initializeWithIdeogram(species = "hg19",
                                chromosome.index = chroms,
                                plotType = c('labels', 'axis'))
  circos.genomicTrack(
    bed,
    track.height = 0.15, ylim = c(3, 12),
    panel.fun = function(region, value, ...) {
      i = getI(...)
      circos.genomicPoints(region, value, col = colors[i], pch=16, cex=0.3, ...)
    })
  circos.yaxis(side = "left", at = seq(-3, 12, by = 3),
               sector.index = get.all.sector.index()[1], labels.cex = 0.4)


  o <- with(bed2, order(Chr, Start, End))
  col_fun = colorRamp2(-(c(1, 0.2, 0.05, 0.01, 0.001) %>% log10),
                       c("grey", "#ffffb2", "#fecc5c", "#fd8d3c", "#e31a1c"))
  circos.genomicHeatmap(bed2[o,] %>% na.omit,
                        col = col_fun,
                        side = "inside",
                        border = NA)

  
  o <- with(bed3, order(Chr %>% gsub("chr", "", .) %>% as.numeric, Start))
  circos.genomicLabels(bed3[o,], labels.column = 4, 
                       side = "inside", cex = 0.5,
                       col = "black", padding = 2, 
                       labels_height = 0.1,
                       connection_height = 0.1)
  

  circos.clear()
}
pA()
###################################################
# Enrichment bias for each dataset 
# except appendix
getEnrichment <- function(dat) {
  myt <- dat[, table(adj.P.Val < 0.05, sign(logFC))]
  pd <- myt %>% as.data.table
  pd <- cbind(pd, 
              myt %>% fisher.test %>% broom::tidy() %>% .[, c("estimate", "p.value")])
  pd <- pd %>% setnames(c("Significant", "Sign", "N", "OR", "P")) %>%
    .[, Significant := factor(Significant, levels = c(FALSE, TRUE), labels = c("p > 0.05", "p < 0.05"))] %>%
    .[, Sign := factor(Sign, levels = c(-1, 1), labels = c("Hypo-\nmethylation", "Hyper-\nmethylation"))] %>% 
    .[order(Significant)] %>%
    .[, F := N / sum(N), Significant]
  pd
}



pd <- list(
  getEnrichment(pfcr)[, Name := "PFC"],
  getEnrichment(pfc)[, Name := "PFC II"],
  getEnrichment(ofb)[, Name := "OFB"]
) %>% rbindlist
pd[, OR := format(OR, digits=3)]
pd[, P := format(P, digits=3, scientific=TRUE)]
pd[, P := gsub("e", "~x~10^", P)]
pd[, Label := paste0("\"OR=\"~", OR, "~\"p=\"~", P)]

pB <- 
  ggplot(pd, aes(Sign, F, 
               fill  = paste(Significant, Sign),
               color = paste(Significant, Sign))) + 
  geom_bar(stat = "identity", position = position_dodge(width=0.5)) +
  geom_text(
    data = pd[Sign == "Hypo-\nmethylation" & Significant == "p < 0.05"],
    x = 1.5, y = 0.77, aes(label = Label),
    size = 2.2, color = "black", parse = TRUE) + 
  ylab("Fraction, %") + 
  facet_wrap(~Name, ncol = 1) + 
  # scale_fill_manual(values = c("blue", "green", "grey", "grey")) + 
  scale_fill_manual("", values=c("#33a02c", "#1f78b4", "white", "white")) + 
  scale_color_manual("", values=c("#33a02c", "#1f78b4", "grey60", "grey60")) + 
  scale_y_continuous(labels = scales::percent_format(accuracy = 2), limits = c(0, .8)) + 
  guides(fill = FALSE, color = FALSE) + 
  theme_publication() + 
  theme(axis.title.x = element_blank())

############################################
# Zoom in specific genes and their neighbourhood
# Need coordinates of the gene
gdt <- list(
  list(Gene = "TLR9", Chr = "chr3", Start = 52255096, End = 52265247),
  list(Gene = "SNCA", Chr = "chr4", Start = 90645250, End = 90759447),
  list(Gene = "ULK1", Chr = "chr12", Start = 132379279, End = 132407707),
  list(Gene = "PINK1", Chr = "chr1", Start = 20959948, End = 20978004))
gdt <- rbindlist(gdt)
gdt[, Gene := factor(Gene, levels = c("ULK1", "TLR9", "PINK1", "SNCA"))]
gdt[, FullName := paste0(Gene, " ", Chr, ":", Start, "-", End)]

smoothSignal <- function(mygene, w=10000) {
  dtA <- 
    pfcr[Gene == mygene, -1 * sign(logFC) * log10(P.Value), SNP] %>%
    .[, ksmooth(SNP, V1, kernel = "normal", bandwidth = w, n.points = 500)] %>%
    .[, Type := "PFC"]
  dtB <- pfc[Gene == mygene, -1 * sign(logFC) * log10(P.Value), SNP] %>%
    .[, ksmooth(SNP, V1, kernel = "normal", bandwidth = w, n.points = 500)] %>%
    .[, Type := "PFC II"]
  dtC <- app[Gene == mygene, -1 * sign(logFC) * log10(P.Value), SNP] %>%
    .[, ksmooth(SNP, V1, kernel = "normal", bandwidth = w, n.points = 500)] %>%
    .[, Type := "APX"]
  dtD <- ofb[Gene == mygene, -1 * sign(logFC) * log10(P.Value), SNP] %>%
    .[, ksmooth(SNP, V1, kernel = "normal", bandwidth = w, n.points = 500)] %>%
    .[, Type := "OFB"]
  dt <- rbindlist(list(dtA, dtB, dtC, dtD))
  dt[, Gene := mygene]  
}


p_load("itertools")
pd <- 
foreach(G = iter(gdt, by = "row"), .combine = rbind) %do% {
  mygene <- G$Gene
  smoothSignal(mygene)
}
pd <- merge(pd, gdt[, list(Gene, FullName)], by = "Gene")

pC <- ggplot() + 
  geom_tile(
    data = pd, 
    aes(x, Type, fill=y, color = y),
    width=1000) + 
  geom_rect(
    data = gdt,
    aes(xmin = Start, xmax = End),
    ymin = 0.5, ymax = 4.5, 
    fill = "white", color = "red", size=0.5, alpha=0
  ) + 
  ylab("") + xlab("") + 
  scale_fill_gradient2("SLP", low = "darkgreen", mid = "white", high = "darkred",
                       limits = c(-1.3, 1.3),
                       aesthetics = c("color", "fill")) + 
  scale_x_continuous(labels = f2si) + 
  theme_publication() + 
  facet_wrap(~ Gene, ncol = 2, scales="free_x") + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank())
  




# Combine plots
setwd("Figure2")
getwd()

# Save circos plot 
png("figure2A.png", width = 7.2 * 0.8, height = 7.2 * 0.8, res = 300, units = "in")
pA()
dev.off()

# quartz("", 7.2, 7.2)
dev.set(4)
legend <- get_legend(pC)
ggdraw() +
  draw_image("figure2A.png", x = 0, y = 0.25, width = 0.8, height = 0.8) + 
  draw_text("APX", x = 0.38, y = 0.855, size = 7) + 
  draw_text("OFB", x = 0.38, y = 0.841, size = 7) +
  draw_text("PFC II", x = 0.38, y = 0.829, size = 7) +
  draw_text("PFC", x = 0.38, y = 0.815, size = 7) +
  draw_plot(pB, x = 0.73, y = 0.4, width = 0.27, height = 0.6) + 
  draw_plot(pC + theme(legend.position = "none"), x = 0,   y = 0, width = 1, height = 0.3) + 
  draw_plot(legend, x = 0.77, y = 0.25, width = 0.2, height = 0.2) + 
  draw_plot_label(
    label = c("A", "B", "C"), 
    size = 10,
    x = c(0, 0.73, 0), 
    y = c(1, 1,    0.3))


ggsave("figure2.png", width = 7.2, height = 7.2, dpi = 300, units = "in")

##########################
# Overlap of targeted genes apx vs the rest
p_load("ggrepel")
p_load("broom")


foo <- function(data, title) {
  pd <- merge(
    app[, sum(adj.P.Val < 0.05, na.rm=TRUE), list(Gene)],
    data[, sum(adj.P.Val < 0.05, na.rm=TRUE), list(Gene)],
    by = c("Gene"))
  t <- pd %>%
    .[, table(V1.x > 0, V1.y > 0)] 
  t %>% 
    fisher.test %>% 
    tidy %>%
    setDT %>% 
    select(estimate, p.value, conf.low, conf.high) %>%
    
  # %>%
  #   .[, Title := title]
}

foo(pfcr, "PFC")

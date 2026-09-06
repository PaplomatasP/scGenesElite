library(ggplot2)
library(grid)
library(gridExtra)
out <- "case-studies/ck-p25"
run <- readRDS(file.path(out,"scmarker_all_cells.rds"))
ranks <- read.csv(file.path(out,"scmarker_all_cells_ranking.csv"))
top <- head(ranks,20)
top$gene <- factor(top$gene,levels=rev(top$gene))
ink <- "#18243B"
base <- theme_minimal(base_size=8,base_family="Arial") +
  theme(text=element_text(colour=ink), panel.grid.minor=element_blank(),
        panel.grid.major.y=element_blank(),
        plot.title=element_text(face="bold",size=10),
        plot.subtitle=element_text(size=7.5,colour="#526178",margin=margin(b=9)),
        axis.text=element_text(colour=ink,size=8),
        plot.margin=margin(6,6,6,6))
p1 <- ggplot(top,aes(x=gene,y=score)) +
  geom_col(fill="#5165CE",width=.68) +
  geom_text(aes(label=score),hjust=-.25,size=2.4,colour=ink) +
  coord_flip() + scale_y_continuous(limits=c(0,340),breaks=c(0,100,200,300),expand=expansion(mult=c(0,0))) +
  labs(title="A  Selected gene ranking",subtitle=paste("Top 20 of", nrow(ranks), "SCMarker-selected genes"),
       x=NULL,y="SCMarker score") + base

genes <- as.character(top$gene)
x <- log2(as.matrix(run$input[,genes,drop=FALSE])+1)
scaled <- scale(x)
groups <- factor(run$metadata$sample,
  levels=c("CK_0w_m1","CKp25_0w_m1","CK_2w_m2","CKp25_2w_m2"))
means <- sapply(levels(groups),function(g) colMeans(scaled[groups==g,,drop=FALSE]))
heat <- as.data.frame(as.table(means),stringsAsFactors=FALSE)
names(heat) <- c("gene","sample","value")
heat$gene <- factor(heat$gene,levels=levels(top$gene))
heat$sample <- factor(heat$sample,levels=levels(groups))
write.csv(heat,file.path(out,"figure_expression_values.csv"),row.names=FALSE)
p2 <- ggplot(heat,aes(x=sample,y=gene,fill=value)) +
  geom_tile(colour="white",linewidth=.35) +
  scale_x_discrete(labels=c("CK\n0 weeks","CK-p25\n0 weeks","CK\n2 weeks","CK-p25\n2 weeks"),expand=c(0,0)) +
  scale_y_discrete(expand=expansion(add=.5)) +
  scale_fill_gradient2(low="#2977A5",mid="white",high="#BA4D68",midpoint=0,
                       limits=c(-1,1),
                       breaks=c(-1,0,1),name="Mean gene-wise z-score") +
  labs(title="B  Expression across sampled groups",subtitle="96 cells per group; log2(FPKM + 1)",x=NULL,y=NULL) +
  base + theme(panel.grid=element_blank(),legend.position="bottom",
               legend.title=element_text(size=7),legend.text=element_text(size=7),
               legend.key.height=unit(2.5,"mm"),legend.key.width=unit(10,"mm"))
# Match plotting panels vertically, including the space used by B's legend.
grDevices::cairo_pdf("case-studies/ck-p25/figures/figure_layout.pdf",width=170/25.4,height=150/25.4,family="Arial")
g1 <- ggplotGrob(p1); g2 <- ggplotGrob(p2)
g1$heights <- g2$heights <- unit.pmax(g1$heights,g2$heights)
panels <- arrangeGrob(g1,g2,ncol=2,widths=c(.45,.55))
# Journal layout: panel headings and plotting area, with methods in the caption.
dev.off()
draw_manuscript <- function() {
  grid.newpage()
  pushViewport(viewport(width=.98,height=.98))
  grid.draw(panels)
  popViewport()
}
figpath <- "case-studies/ck-p25/figures/scGenesFinder_case_study"
grDevices::cairo_pdf(paste0(figpath,".pdf"),width=170/25.4,height=125/25.4,family="Arial")
draw_manuscript(); dev.off()
ragg::agg_png(paste0(figpath,"_600dpi.png"),width=170,height=125,units="mm",res=600,background="white")
draw_manuscript(); dev.off()
ragg::agg_tiff(paste0(figpath,"_600dpi.tiff"),width=170,height=125,units="mm",res=600,compression="lzw",background="white")
draw_manuscript(); dev.off()
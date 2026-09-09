# Publication plots share numerical results between screen, individual exports
# and the composite. No dataset-specific genes, groups or metrics are embedded.
publication_theme <- function() {
  ggplot2::theme_minimal(base_size=8) + ggplot2::theme(
    text=ggplot2::element_text(colour="#18243B"), panel.grid.minor=ggplot2::element_blank(),
    plot.title=ggplot2::element_text(size=10,face="bold"),
    plot.subtitle=ggplot2::element_text(size=7,colour="#526178",margin=ggplot2::margin(b=8)),
    axis.text=ggplot2::element_text(size=7,colour="#18243B"),
    plot.margin=ggplot2::margin(8,12,8,8),legend.position="bottom")
}

publication_data <- function(result, n, transform="log2", metadata=NULL, inspect_genes=NULL) {
  ctx <- result$plot_context
  if (is.null(ctx)) stop("Run the analysis again to capture the input and mapping counts.")
  n <- as.integer(n)
  if (length(n)!=1L || is.na(n) || n<2L || n>50L)
    stop("Choose between 2 and 50 ranked genes for a readable publication figure.")
  genes <- head(rownames(result$ig), n)
  expression <- ctx$expression
  if (!length(genes) || !all(genes %in% colnames(expression)))
    stop("Some ranked genes could not be matched to the mapped expression matrix.")
  x <- as.matrix(expression[,genes,drop=FALSE])
  if (!is.numeric(x) || any(!is.finite(x))) stop("Expression values must be finite numbers.")
  if (transform=="log2") {
    if (any(x<0)) stop("Negative expression values: select 'As supplied' instead of log2(x + 1).")
    x <- log2(x+1)
  } else if (transform!="none") stop("Unknown expression transformation.")
  cells <- rownames(expression)
  labels <- as.character(expression[[ncol(expression)]])
  group <- labels
  grouping <- "label"
  if (!is.null(metadata)) {
    if (!all(c("cell","sample") %in% names(metadata))) stop("Metadata CSV needs columns named cell and sample.")
    metadata$cell <- as.character(metadata$cell)
    metadata$sample <- as.character(metadata$sample)
    if (anyNA(metadata$cell) || anyDuplicated(metadata$cell)) stop("Metadata cell identifiers must be unique and non-missing.")
    idx <- match(cells,metadata$cell)
    if (anyNA(idx)) stop("Metadata must contain every analyzed cell identifier.")
    group <- metadata$sample[idx]
    if (anyNA(group) || any(!nzchar(trimws(group)))) stop("Every analyzed cell needs a sample name.")
    grouping <- "sample"
  }
  groups <- if(is.null(metadata)) unique(group) else unique(metadata$sample[metadata$cell %in% cells])
  if (length(groups)>20L) stop("Use at most 20 display groups for a readable figure.")
  scaled <- scale(x)
  # A constant gene has no expression contrast; represent it by zeros.
  constant <- !is.finite(attr(scaled,"scaled:scale")) | attr(scaled,"scaled:scale")==0
  if (any(constant)) scaled[,constant] <- 0
  means <- vapply(groups,function(g) colMeans(scaled[group==g,,drop=FALSE]),numeric(length(genes)))
  dimnames(means) <- list(genes,groups)
  heat <- as.data.frame(as.table(means),stringsAsFactors=FALSE)
  names(heat) <- c("gene","group","value")
  heat$gene <- factor(heat$gene,levels=rev(genes)); heat$group <- factor(heat$group,levels=groups)
  if (is.null(inspect_genes) || !length(inspect_genes)) inspect_genes <- head(genes,3)
  if (length(inspect_genes)>3L || !all(inspect_genes %in% genes))
    stop("Choose one to three genes from the displayed ranking for the expression plots.")
  boxes <- do.call(rbind,lapply(inspect_genes,function(g) data.frame(
    cell=cells,gene=g,group=factor(group,levels=groups),label=labels,value=x[,g])))
  boxes$gene <- factor(boxes$gene,levels=inspect_genes)
  counts <- data.frame(stage=factor(c("Input","Mapped","Selected","Displayed"),
      levels=c("Displayed","Selected","Mapped","Input")),
    n=c(ctx$input_genes,ctx$mapped_genes,nrow(result$ig),length(genes)))
  list(heat=heat,boxes=boxes,counts=counts,genes=genes,
    groups=data.frame(group=groups,cells=as.integer(table(factor(group,levels=groups)))),
    grouping=grouping,transform=transform,constant=genes[constant],method=ctx$method,
    cells=length(cells),n=length(genes),inspect_genes=inspect_genes)
}

publication_plots <- function(d, classifier) {
  theme <- publication_theme()
  reduce <- 100*(1-d$counts$n[3]/d$counts$n[2])
  p1 <- ggplot2::ggplot(d$counts,ggplot2::aes(stage,n)) +
    ggplot2::geom_col(ggplot2::aes(fill=stage),width=.5,show.legend=FALSE) +
    ggplot2::scale_fill_manual(values=c(Input="#D6DCE7",Mapped="#A9B6D4",Selected="#5267CD",Displayed="#18243B")) +
    ggplot2::geom_text(ggplot2::aes(label=format(n,big.mark=",")),hjust=-.15,size=2.7) +
    ggplot2::coord_flip() + ggplot2::scale_y_continuous(expand=ggplot2::expansion(mult=c(0,.2))) +
    ggplot2::labs(title="A  From input to selected genes",
      subtitle=sprintf("%s | %d cells | %d %s groups",d$method,d$cells,nrow(d$groups),d$grouping),
      x=NULL,y="Number of genes",caption=sprintf("%.1f%% fewer genes after selection (%s to %s)",
        reduce,format(d$counts$n[2],big.mark=","),format(d$counts$n[3],big.mark=","))) +
    theme + ggplot2::theme(panel.grid.major.y=ggplot2::element_blank(),plot.caption=ggplot2::element_text(size=7,hjust=0))
  lim <- max(1,ceiling(max(abs(d$heat$value))*10)/10)
  p2 <- ggplot2::ggplot(d$heat,ggplot2::aes(group,gene,fill=value)) +
    ggplot2::geom_tile(colour="white",linewidth=.3) +
    ggplot2::scale_fill_gradient2(low="#3C82A5",mid="white",high="#BE5C73",limits=c(-lim,lim),
      breaks=c(-lim,0,lim),name="Mean gene-wise z-score") +
    ggplot2::scale_x_discrete(expand=c(0,0)) +
    ggplot2::labs(title=paste("B  Expression context of the top",d$n),
      subtitle=paste("Means by",d$grouping,"after gene-wise standardization"),x=NULL,y=NULL) +
    theme + ggplot2::theme(panel.grid=ggplot2::element_blank(),axis.text.y=ggplot2::element_text(face="italic"),
      axis.text.x=ggplot2::element_text(angle=30,hjust=1),legend.title=ggplot2::element_text(size=7),
      legend.key.height=grid::unit(2,"mm"),legend.key.width=grid::unit(9,"mm"))
  labs <- unique(d$boxes$label)
  colours <- if (length(labs)==2L) c("#3C82A5","#BE5C73") else grDevices::hcl.colors(length(labs),"Dark 3")
  p3 <- ggplot2::ggplot(d$boxes,ggplot2::aes(group,value,fill=label)) +
    ggplot2::geom_boxplot(width=.6,outlier.size=.35,linewidth=.3) +
    ggplot2::facet_wrap(~gene,nrow=1,scales="free_y") +
    ggplot2::scale_fill_manual(values=stats::setNames(colours,labs),name="Label") +
    ggplot2::labs(title="C  Inspect individual expression patterns",
      subtitle=if(length(unique(d$groups$cells))==1L)
        sprintf("Cell distributions | %d cells per %s group",d$groups$cells[1],d$grouping)
        else sprintf("Cell distributions | %d-%d cells per %s group",min(d$groups$cells),max(d$groups$cells),d$grouping),
      x=paste("Display",d$grouping),y=if(d$transform=="log2") "log2(expression + 1)" else "Expression as supplied") +
    theme + ggplot2::theme(panel.grid.major.x=ggplot2::element_blank(),
      axis.text.x=ggplot2::element_text(angle=45,hjust=1),strip.text=ggplot2::element_text(face="italic"))
  if (is.null(classifier$error)) {
    p4 <- ConfMatrixPlot(classifier$confusionMatrix,subtitle=sprintf("Top %d genes | k-NN | accuracy %.1f%%",
      length(classifier$genes),100*classifier$confusionMatrix$overall['Accuracy']))
  } else {
    p4 <- ggplot2::ggplot() + ggplot2::annotate("text",x=0,y=0,
      label=paste(strwrap(classifier$error, width=48),collapse="\n"),size=3) +
      ggplot2::labs(title="D  Classification unavailable") + theme +
      ggplot2::theme(axis.text=ggplot2::element_blank(),axis.title=ggplot2::element_blank(),panel.grid=ggplot2::element_blank())
  }
  list(selection=p1,heatmap=p2,expression=p3,classification=p4)
}

draw_publication <- function(plots, panel="composite") {
  if(panel=="composite") {
    grid::grid.newpage()
    grid::grid.draw(gridExtra::arrangeGrob(grobs=unname(plots),ncol=2,widths=c(1.05,1),heights=c(1.3,1)))
  } else print(plots[[panel]])
}

export_publication_plot <- function(file, plots, panel="composite", format="pdf", dpi=600, n=20) {
  width <- if(panel=="composite") 180 else 160
  height <- if(panel=="composite") max(180,120+n*3) else if(panel=="heatmap") max(125,50+n*3) else 125
  if(format=="pdf") {
    if(capabilities("cairo")) grDevices::cairo_pdf(file,width=width/25.4,height=height/25.4)
    else grDevices::pdf(file,width=width/25.4,height=height/25.4)
  } else if(format=="png") {
    if(!dpi %in% c(300,600)) stop("Choose 300 or 600 dpi.")
    if(requireNamespace("ragg",quietly=TRUE)) ragg::agg_png(file,width=width,height=height,units="mm",res=dpi,background="white")
    else grDevices::png(file,width=width,height=height,units="mm",res=dpi,bg="white")
  } else stop("Choose PDF or PNG.")
  on.exit(grDevices::dev.off(),add=TRUE)
  draw_publication(plots,panel)
}

write_publication_bundle <- function(file, bundle, dpi=600) {
  file <- file.path(normalizePath(dirname(file),mustWork=TRUE),basename(file))
  directory <- tempfile("scgenes-publication-"); dir.create(directory)
  on.exit(unlink(directory,recursive=TRUE),add=TRUE)
  for(panel in c("composite",names(bundle$plots))) for(format in c("pdf","png"))
    export_publication_plot(file.path(directory,paste0(panel,".",format)),bundle$plots,panel,format,dpi,bundle$data$n)
  for(name in c("heat","boxes","counts","groups"))
    utils::write.csv(bundle$data[[name]],file.path(directory,paste0(name,".csv")),row.names=FALSE)
  utils::write.csv(bundle$ranking,file.path(directory,"ranking.csv"),row.names=FALSE)
  utils::write.csv(bundle$filtered,file.path(directory,"FilterData.csv"))
  if(is.null(bundle$classifier$error)) {
    utils::write.csv(bundle$classifier$predictions,file.path(directory,"predictions.csv"),row.names=FALSE)
    utils::write.csv(as.data.frame(bundle$classifier$confusionMatrix$table),file.path(directory,"confusion.csv"),row.names=FALSE)
    writeLines(capture.output(bundle$classifier$confusionMatrix),file.path(directory,"classification_metrics.txt"))
    saveRDS(bundle$classifier,file.path(directory,"classification.rds"))
  }
  writeLines(c(sprintf("Method: %s; displayed genes: %d",bundle$data$method,bundle$data$n),
    paste("Expression transform:",bundle$data$transform),paste("Grouping:",bundle$data$grouping),
    "Heatmap: transform mapped input expression, standardize each gene across cells, then average by group.",
    paste("Constant genes displayed as zero:",paste(bundle$data$constant,collapse=", ")),
    "Boxplots: median, interquartile range, 1.5 IQR whiskers; cells are not independent biological replicates.",
    paste("Inspected genes:",paste(bundle$data$inspect_genes,collapse=", ")),
    paste("Classifier seed:",bundle$seed),
    "k-NN: selected predictors, 80% stratified training split, center/scale, up to five-fold training CV.",
    "Classification uses the selected method's output matrix; expression plots use the mapped input values.",
    "Selection precedes the cell split; independent samples are not held out. This is not independent validation.",
    if(!is.null(bundle$classifier$error)) paste("Classification unavailable:",bundle$classifier$error),
    capture.output(sessionInfo())),file.path(directory,"README.txt"))
  zip::zip(file,files=list.files(directory),root=directory)
}

register_publication_results <- function(input,output,session,analysis_result,run_state) {
  result <- shiny::reactive({shiny::req(run_state$ready); shiny::req(analysis_result()); analysis_result()})
  metadata <- shiny::reactive({
    if(is.null(input$publicationGrouping) || input$publicationGrouping=="labels") return(NULL)
    shiny::validate(shiny::need(!is.null(input$publicationMetadata),"Upload a cell-to-sample metadata CSV."))
    utils::read.csv(input$publicationMetadata$datapath,check.names=FALSE,stringsAsFactors=FALSE)
  })
  shiny::observeEvent(list(result(),input$genes),{
    genes <- head(rownames(result()$ig),min(50,input$genes))
    chosen <- shiny::isolate(input$publicationInspect)
    chosen <- intersect(chosen,genes)
    if(!length(chosen)) chosen <- head(genes,3)
    shiny::updateSelectizeInput(session,"publicationInspect",choices=genes,selected=head(chosen,3),server=TRUE)
  })
  classifier <- shiny::reactive({
    r <- result(); n <- input$genes; seed <- input$publicationSeed
    tryCatch({
      if(length(seed)!=1L || !is.finite(seed) || seed<0 || seed>.Machine$integer.max)
        stop("Enter a valid non-negative integer seed.")
      KnnClassifier(r$newdata,r$ig,r$newdata[[ncol(r$newdata)]],genes_count=n,seed=as.integer(seed),fit_only=TRUE)
    },error=function(e) list(error=conditionMessage(e)))
  })
  data <- shiny::reactive({
    tryCatch(publication_data(result(),input$genes,input$publicationTransform,metadata(),input$publicationInspect),
      error=function(e) shiny::validate(shiny::need(FALSE,conditionMessage(e))))
  })
  plots <- shiny::reactive(publication_plots(data(),classifier()))
  bundle <- shiny::reactive({
    r <- result()
    list(data=data(),plots=plots(),classifier=classifier(),seed=input$publicationSeed,
      filtered=r$newdata,ranking=data.frame(gene=rownames(r$ig),score=r$ig[[1]]))
  })
  output$SelectionSummary <- shiny::renderPlot({print(plots()$selection)},res=120)
  output$HeatMap <- shiny::renderPlot({print(plots()$heatmap)},res=120,
    height=function() max(520,min(50,input$genes)*22))
  output$ExpressionDistributions <- shiny::renderPlot({print(plots()$expression)},res=120)
  # The classifier does not depend on heatmap grouping or expression transforms.
  output$KnnClassifier <- shiny::renderPlot({
    fit <- classifier()
    shiny::validate(shiny::need(is.null(fit$error),fit$error))
    print(ConfMatrixPlot(fit$confusionMatrix,subtitle=sprintf("Top %d genes | k-NN | accuracy %.1f%%",
      length(fit$genes),100*fit$confusionMatrix$overall['Accuracy'])))
  },res=120)
  output$HeatmapList <- DT::renderDataTable(DT::datatable(data.frame(Gene=data()$genes),rownames=FALSE))
  output$PublicationGroups <- shiny::renderTable(data()$groups)
  snapshot <- shiny::reactiveVal(NULL)
  # Clear previews whenever a setting or the underlying run changes.
  shiny::observeEvent(list(analysis_result(),run_state$ready,input$genes,input$publicationGrouping,
    input$publicationMetadata,input$publicationTransform,input$publicationInspect,input$publicationSeed),
    {snapshot(NULL)},priority=20,ignoreNULL=FALSE)
  shiny::observeEvent(input$generatePublication,{
    snapshot(bundle())
  })
  output$publicationReady <- function(...) !is.null(snapshot()) && isTRUE(run_state$ready)
  shiny::outputOptions(output,"publicationReady",suspendWhenHidden=FALSE)
  output$PublicationComposite <- shiny::renderPlot({shiny::req(snapshot());draw_publication(snapshot()$plots)},res=120)
  output$downloadPublication <- shiny::downloadHandler(
    filename=function() paste0("scGenesFinder-",input$publicationPanel,".",input$publicationFormat),
    content=function(file) {
      b <- if(input$publicationPanel=="composite") {shiny::req(snapshot());snapshot()} else bundle()
      export_publication_plot(file,b$plots,input$publicationPanel,input$publicationFormat,
        as.integer(input$downloadDpi),b$data$n)
    })
  output$downloadData <- shiny::downloadHandler(
    filename=function() paste0("scGenesFinder-results-",Sys.Date(),".zip"),
    content=function(file) write_publication_bundle(file,bundle(),as.integer(input$downloadDpi)),
    contentType="application/zip")
  invisible(list(data=data,classifier=classifier,bundle=bundle,snapshot=snapshot))
}

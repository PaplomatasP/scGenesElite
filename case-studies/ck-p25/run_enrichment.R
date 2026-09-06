# Offline KEGG over-representation with the actual mapped input as the universe.
# This is a publication audit, not a reproduction of Enrichr's default background.
out <- "case-studies/ck-p25"
sets <- strsplit(readLines(file.path(out,"KEGG_2019_Mouse.gmt"), warn=FALSE), "\t", fixed=TRUE)
# The downloaded GMT ends with an empty tab-only record, not a pathway.
sets <- Filter(function(z) length(z) >= 3L && nzchar(trimws(z[[1]])), sets)
terms <- vapply(sets, `[[`, character(1), 1)
sets <- setNames(lapply(sets, function(z) unique(toupper(z[-c(1,2)]))), terms)
for (scope in c("all_cells","week2")) {
  run <- readRDS(file.path(out,paste0("scmarker_",scope,".rds")))
  universe <- unique(toupper(names(run$input)[-ncol(run$input)]))
  for (cutoff in c(20,50,100)) {
    selected <- toupper(head(rownames(run$result$ig), cutoff))
    tab <- do.call(rbind,lapply(names(sets),function(term) {
      members <- intersect(sets[[term]],universe)
      overlap <- intersect(selected,members)
      data.frame(term=term, selected=length(selected), universe=length(universe),
        pathway_in_universe=length(members), overlap=length(overlap),
        p_value=phyper(length(overlap)-1,length(members),length(universe)-length(members),length(selected),lower.tail=FALSE),
        genes=paste(overlap,collapse=";"))
    }))
    tab$fdr <- p.adjust(tab$p_value,method="BH")
    tab <- tab[order(tab$fdr,tab$p_value,tab$term),]
    write.csv(tab,file.path(out,paste0("kegg_",scope,"_top",cutoff,".csv")),row.names=FALSE)
    cat(scope,cutoff,"FDR < .05:",sum(tab$fdr<.05),"\n")
    print(head(tab[,c("term","overlap","p_value","fdr")],8),row.names=FALSE)
  }
}

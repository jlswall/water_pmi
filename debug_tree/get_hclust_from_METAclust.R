get.hclust.from.METAclust <- function (meta, group = 4, data.trans = NULL,
                                       dist = NULL, clust = NULL,
                                       type = NULL, main = "", file = NULL,
                                       ext = NULL, width = 8, height = 8) 
{
    save <- !is.null(file)
    if (save) {
        .get.dev(file, ext, height = height, width = width)
    }
    suppressWarnings(meta.new <- filter.META(meta))
    cl <- vector()
    fac <- vector()
    num <- vector()
    for (i in 1:ncol(meta.new)) {
        if (class(meta.new[, i]) == "factor" || class(meta.new[, 
            i]) == "character") {
            fac <- unique(c(fac, i))
        }
        if (class(meta.new[, i]) == "numeric") {
            num <- unique(c(num, i))
        }
        cl <- unique(c(cl, class(meta.new[, i])))
    }
    if (is.null(data.trans)) {
        meta.new.tr <- meta.new
    }
    else {
        if (length(num) != 0 && length(fac) != 0) {
            meta.new.num <- vegan::decostand(meta.new[, num], 
                data.trans)
            meta.new.tr <- cbind(meta.new[, fac], meta.new.num)
            names(meta.new.tr) <- c(names(meta.new)[fac], names(meta.new)[num])
        }
        else {
            warning("no numeric variables to be transformed")
            meta.new.tr <- meta.new
        }
    }
    if (any(cl %in% c("factor", "character"))) {
        dis <- FD::gowdis(meta.new.tr)
    }
    else {
        dis <- vegan::vegdist(meta.new.tr, method = dist)
    }
    if (is.null(clust)) {
        hc <- stats::hclust(dis)
    }
    else {
        hc <- stats::hclust(dis, clust)
    }

    return(hc)
    ## if (length(group) != 1L) {
    ##     warning(" group must be a digit or a metadata variable; \n            will cut tree in 4 groups as default")
    ##     mem <- stats::cutree(hc, 4)
    ## }
    ## else {
    ##     if (is.numeric(group)) {
    ##         mem <- stats::cutree(hc, group)
    ##     }
    ##     else {
    ##         if (any(group %in% names(meta.new.tr))) {
    ##             mem <- factor(as.numeric(factor(meta.new.tr[[group]])))
    ##             attr(mem, "names") <- rownames(meta.new.tr)
    ##         }
    ##         else {
    ##             warning(" group is NOT a metadata variable, \n                will cut tree in 4 groups as default")
    ##             mem <- stats::cutree(hc, 4)
    ##         }
    ##     }
    ## }
    ## .pCLUST(hc, type = type, mem = mem)
    ## if (save) {
    ##     dev.off()
    ## }
    ## invisible()
}

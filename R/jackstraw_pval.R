#' Get PC scores
#'
#' Calculate scores for PC based on jackstraw results
#' @return A dataframe of PCs and their p-values
#' @importFrom reshape2 melt
#' @param object The seurat object
#' @param PCs PCs to select from
#' @param score.thresh emperical p-value cut-off to apply
#' @author Matthew Myint, Jinmiao Chen
#' @export
#'
#' @examples
#' pc_pval <- JackStraw_pval(object = sObj)
JackStraw_pval <- function (object, PCs = 1:20, score.thresh = 1e-05) {
    pAll <- GetDimReduction(object, reduction.type = "pca", slot = "jackstraw")@emperical.p.value
    pAll <- pAll[, PCs, drop = FALSE]
    pAll <- as.data.frame(pAll)
    pAll$Contig <- rownames(x = pAll)
    pAll.l <- melt(data = pAll, id.vars = "Contig")
    colnames(x = pAll.l) <- c("Contig", "PC", "Value")
    qq.df <- NULL
    score.df <- NULL
    for (i in PCs) {
        q <- qqplot(x = pAll[, i], y = runif(n = 1000), plot.it = FALSE)
        pc.score <- suppressWarnings(
            prop.test(x = c(length(x = which(x = pAll[, i] <= score.thresh)),
                            floor(x = nrow(x = pAll) * score.thresh)),
                      n = c(nrow(pAll), nrow(pAll))
            )$p.val
        )
        if (length(x = which(x = pAll[, i] <= score.thresh)) == 0) {
            pc.score <- 1
        }
        if (is.null(x = score.df)) {
            score.df <- data.frame(PC = paste0("PC", i), Score = pc.score)
        }else{
            score.df <- rbind(score.df, data.frame(PC = paste0("PC", i), Score = pc.score))
        }
        if (is.null(x = qq.df)) {
            qq.df <- data.frame(x = q$x, y = q$y, PC = paste0("PC", i))
        }else{
            qq.df <- rbind(qq.df, data.frame(x = q$x, y = q$y, PC = paste0("PC", i)))
        }
    }
    return(score.df)
}

#' PC filtering
#'
#' Choose significant PCs based on their p-value for use in TSNE
#' @return A numeric vector of significant PCs
#' @param score.df A dataframe of PCs and their scores
#' @param p.val p-value cut-off for significant PCs
#' @author Matthew Myint, Jinmiao Chen
#' @export
#'
#' @examples
#' score.df <- data.frame(PC = paste0("PC", 1:10), Score = runif(10, 0, 0.2))
#' pc_pval <- signif_PCs(object = sObj)
signif_PCs <- function (score.df, p.val = 0.05) {
    s.PCs <- which(score.df$Score < p.val)
    return(s.PCs)
}

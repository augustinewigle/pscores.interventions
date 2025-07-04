# Augustine's functions

# Function for running many values of CIVs -------------

#' @param x a list of netmeta objects
#' @param CIVs a named list of CIVs to explore. Each name corresponds to an outcome. Order should match x
#' @param correlation the correlation matrix describing the correlation between the outcomes, or NULL to assume zero correlation
#' @param type a vector describing the type of outcomes - "H" for harmful (smaller values good), "B" for beneficial (larger values good)
#' CIVs is a NAMED LIST
pscore_civs <- function(x,CIVs,correlation,type) {

  prepare_data<-prep(x)
  outcomes<-prepare_data$outcomes
  var.outcomes<-prepare_data$var.outcomes
  comm<-prepare_data$comm

  if(length(CIVs) != dim(outcomes)[3]) {

    stop("CIVs must be a list with the same length as the number of outcomes")

  }

  CIV_mat <- expand.grid(CIVs)

  pscore_df <- matrix(nrow = nrow(CIV_mat), ncol = length(comm), dimnames = list(list(), trts = comm))

  for(i in 1:nrow(CIV_mat)) {

    pscore_df[i,] <- pscores(outcomes,
                             var.outcomes,
                             correlation,
                             beta = as.numeric(-CIV_mat[i,]),
                             type = type,
                             label = as.vector(comm))

  }

  res <- cbind(pscore_df, CIV_mat) %>%
    pivot_longer(cols = 1:length(comm), names_to = "Treatment", values_to = "Pscore") %>%
    group_by(pick(names(CIVs))) %>% mutate(ranking = rank(-Pscore), poth = poth2(Pscore))

  return(list(all = res, pscores = pscore_df, CIVs = CIV_mat, labels = comm))


}

# Summary measures for CIV ranges --------------------------------------

#' @param x an object from running pscore_civs
abpmc <- function(x) {

  df <- x$pscores
  civs <- x$CIVs

  labels <- colnames(df)

  # check
  if(nrow(df) != nrow(civs)) {

    stop("Number of rows of df and civs must be equal")

  } else if(nrow(civs) < 2) {

    stop("More than one combo of CIVs must be considered to calculate the area")

  }

  # Drop the rows which only have one CIV

  nval <- apply(civs, 2, function(x) length(unique(x)))

  civs <- as.matrix(civs[,which(nval>1)])

  hs <- apply(civs, 2, function(x) (max(x)-min(x))/length(unique(x)))

  pscoremeans <- apply(df, 1, mean)

  auc <- apply(df-pscoremeans, 2, function(x, h) sum(x*prod(h)), h = hs)

  names(auc) <- labels

  return(data.frame(ABPMC = auc, Treatment = labels)) # average probability that it is better than the mean, averaged over the CIVs

}

#' @param x an object from running pscore_civs
aupc <- function(x) {

  df <- x$pscores
  civs <- x$CIVs

  labels <- colnames(df)

  # check
  if(nrow(df) != nrow(civs)) {

    stop("Number of rows of df and civs must be equal")

  } else if(nrow(civs) < 2) {

    stop("More than one combo of CIVs must be considered to calculate the area")

  }

  # Drop the rows which only have one CIV

  nval <- apply(civs, 2, function(x) length(unique(x)))

  civs <- as.matrix(civs[,which(nval>1)])

  hs <- apply(civs, 2, function(x) (max(x)-min(x))/length(unique(x)))

  auc <- apply(df, 2, function(x, h) sum(x*prod(h)), h = hs)

  names(auc) <- labels

  return(data.frame(AUPC = auc, Treatment = labels)) # average probability that it beats all other treatments, averaged over the CIVs

}

#' @param obj an object from running either aupc() or abpmc()
#' @param ordering either "ranking" for treatments on the x axis to be ordered by the measure, or "asis" to keep the same order of treatments as in obj
#' @param title optional title for the plot
plot_civsummary <- function(obj, ordering = "ranking", title = "") {

  measure <- colnames(obj)[1]

  if(ordering == "ranking") { # reorder factor so the x axis is in order of least to most preferred

    obj$Treatment <- reorder(obj$Treatment, obj[,colnames(obj) == measure])

  }

  ggplot(obj, aes(x = Treatment, y = !!sym(measure))) +
    geom_col(col = "black", fill = "skyblue3") +
    geom_hline(yintercept = 0) +
    labs(title = title) +
    theme_bw()

}

# Plotting for 1 outcome ------------------------------

#' @param x object from running pscore_civs. only accommodates objects where only one outcome has varying CIV
#' @param room a vector of length two that can be used to expand the plot area to accommodate the legend
#' @param title optional title for the plot
plot_pscores <- function(x, room = c(0.05, 0.05), title = "") {

  pscores <- x$pscores
  CIVs <- x$CIVs[,1]
  outcome_name <- names(x$CIVs)[1]

  if(ncol(x$CIVs) > 1) {

    stop("This function is meant for single-outcome P-scores")

  }

  labels <- colnames(pscores)
  highlight_ix <- sort(pscores[1,], index.return = TRUE, decreasing = TRUE)$ix[1:5] # highlight the top 5

  plot(0,0,ylim=c(min(pscores, 0), max(pscores, 1)+room[2]),
       type="l",
       xlim=c(min(CIVs),max(CIVs)+room[1]),
       main = title,
       xlab=paste0("CIV for ", outcome_name),ylab="P-score")

  for(i in 1:ncol(pscores)) {

    col_ix <- ifelse(i %in% highlight_ix, which(highlight_ix == i)+1, 1)
    lines(CIVs,pscores[,i],col=col_ix)

  }

  lines(CIVs, rowMeans(pscores), type = "b")

  legend(x = "topright", title = paste0("Top treatments\n(at CIV = ", round(CIVs[1], 1), ")"),
         legend = c(labels[highlight_ix], "Mean"),
         col = c(2:(6), 1),
         lty = c(rep(1, 5), NA),
         pch = c(rep(NA, 5), 1),
         # bty = "n",
         horiz= F)
}

#' @param x object from running pscore_civs()
#' @param title optional title for the plot
plot_pothciv <- function(x, title = "") {

  # check if multiple outcomes, only 1 can have varying CIVs

  df <- x$all
  civs <- x$CIVs

  if(ncol(civs)>1) {

    nval <- apply(civs, 2, function(x) length(unique(x)))

    if(length(which(nval>1)) == 1) {

      outcome <- colnames(civs)[which(nval>1)]

      df <- select(ungroup(df), poth, !!sym(outcome)) %>% summarise(poth = unique(poth), .by = !!sym(outcome))

    } else {

      stop("Function only compatible for objects where only one outcome has a range of CIVs")

    }

  } else {

    outcome <- colnames(civs)[1]

  }

  ggplot(df, aes(x = !!sym(outcome), y = poth)) +
    geom_point(col = "black", shape = 21, size = 2, fill = "hotpink") +
    theme_bw() +
    geom_hline(yintercept = 0) +
    labs(x = paste0("CIV (", colnames(civs)[1], ")"), y = "POTH", title = title)



}

# Functions for multiple outcomes -------------------------------------

#' @param x object from running pscore_civs
#' @param title optional title fo rplot
pscores_heatplot <- function(x, title = "") {

  if(ncol(x$CIVs) != 2) {

    stop("Only implemented for 2 outcomes")

  }

  x$all$Treatment <- factor(x$all$Treatment, levels = x$all$Treatment[order(colMeans(x$pscores), decreasing = TRUE)])

  ggplot(x$all, aes(x = !!sym(names(x$CIVs)[1]),
                y = !!sym(names(x$CIVs)[2]),
                fill = ranking)) +
    facet_wrap(~Treatment) +
    geom_tile() +
    scale_fill_viridis_c() +
    labs(x = paste0("CIV (", names(x$CIVs)[1], ")"),
         y = paste0("CIV (", names(x$CIVs)[2], ")"),
         fill = "Rank based on P-score",
         title = title)

}

#' @param x object from running pscore_civs. Should include exactly 2 outcomes
#' @param title optional title for plot
pscores_pothplot <- function(x, title) {

  if(ncol(x$CIVs) != 2) {

    stop("FUnction only designed for 2-outcome p-scores right now")

  }

  bub <- x$all %>%summarise(POTH = unique(poth), .groups = "keep")

  ggplot(bub, aes(x = !!sym(names(x$CIVs)[1]),
                  y = !!sym(names(x$CIVs)[2]),
                  size = POTH)) +
    geom_point(alpha = 0.7, shape = 21, col = "black", fill = "hotpink4") +
    scale_size(range = c(.5, 20)) +
    labs(x = paste0("CIV (", names(x$CIVs)[1], ")"),
         y = paste0("CIV (", names(x$CIVs)[2], ")"),
         title = title) +
    theme_bw()


}

# Helpers -----------------
poth2 <- function(pscores) {

  n <- length(pscores)

  sq <- pscores - mean(pscores)
  p <- sum(sq*sq)/n*12*(n-1)/(n+1)

  return(p)

}

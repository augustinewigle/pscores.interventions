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

  pscore_df <- matrix(nrow = nrow(CIV_mat), ncol = length(comm), dimnames = list(list(), comm))

  for(i in 1:nrow(CIV_mat)) {

    pscore_df[i,] <- pscores(outcomes=outcomes,
                             var.outcomes=var.outcomes,
                             correlation = correlation,
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

  ranking <- rank(-auc)

  return(data.frame(ABPMC = auc, Treatment = labels, ranking = ranking)) # average probability that it is better than the mean, averaged over the CIVs

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

  ranking <- rank(-auc)

  return(data.frame(AUPC = auc, Treatment = labels, ranking = ranking)) # average probability that it beats all other treatments, averaged over the CIVs

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
#' @param residuals Logical - should P-score residuals or raw p-scores are plotted? Default is TRUE
#' @param room a vector of length two that can be used to expand the plot area to accommodate the legend
#' @param title optional title for the plot
plot_pscores <- function(x, residuals = TRUE, room = c(0.05, 0.05), title = "") {

  if(ncol(x$CIVs) > 1) {

    stop("This function is meant for single-outcome P-scores")

  }

  if(residuals) {

    pscores <- x$pscores-mean(x$pscores)

  } else {

    pscores <- x$pscores

  }

  CIVs <- x$CIVs[,1]
  outcome_name <- names(x$CIVs)[1]

  labels <- colnames(pscores)
  highlight_ix <- sort(pscores[1,], index.return = TRUE, decreasing = TRUE)$ix[1:5] # highlight the top 5

  plot(0,0,ylim=c(min(pscores, 0), max(pscores, 1)+room[2]),
       type="l",
       xlim=c(min(CIVs),max(CIVs)+room[1]),
       main = title,
       xlab=paste0("CIV for ", outcome_name),ylab=ifelse(residuals, "Residual", "P-score"))

  for(i in 1:ncol(pscores)) {

    col_ix <- ifelse(i %in% highlight_ix, which(highlight_ix == i)+1, 1)
    lines(CIVs,pscores[,i],col=col_ix)

  }

  # Add lines to help visualize

  if(!residuals) {

    lines(CIVs, rowMeans(pscores), type = "b")

  } else {

    abline(h = 0, lty = "dashed")

  }

  # set up legend
  leglabs <- labels[highlight_ix]
  cols <- 2:6
  ltys <- rep(1,5)
  pchs <- rep(NA, 5)

  if(!residuals) {

    leglabs <- c(leglabs, "Mean")
    cols <- c(cols, 1)
    ltys <- c(ltys, NA)
    pchs <- c(pchs, 1)

  }

  legend(x = "topright", title = paste0("Top treatments\n(at CIV = ", round(CIVs[1], 1), ")"),
         legend = leglabs,
         col = cols,
         lty = ltys,
         pch = pchs,
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
    geom_line() +
    geom_point(col = "black", shape = 21, size = 2.5, fill = "hotpink") +
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
         fill = "Rank based on Extended P-score",
         title = title)

}

#' @param x object from running pscore_civs. Should include exactly 2 outcomes
#' @param newgridsize Grid size. 0 indicates to use the grid size of the original object
#' @param hightlight logical; should the minimum and maximum POTH values in the new grid be highlighted?
#' @param title optional title for plot
pscores_pothplot <- function(x, newgridsize = 0, highlight = TRUE, title = "POTH plot") {

  if(ncol(x$CIVs) != 2) {

    stop("Function only designed for 2-outcome p-scores right now")

  }

  bub <- x$all %>%summarise(POTH = unique(poth), .groups = "keep")

  if(newgridsize > 0) {

    keepx <- quantile(x$CIVs[,1], probs = 0:(newgridsize-1)/(newgridsize-1), type = 1)
    keepy <- quantile(x$CIVs[,2], probs = 0:(newgridsize-1)/(newgridsize-1), type = 1)

    ix <- which(as.data.frame(bub[,1])[,1] %in% keepx & as.data.frame(bub[,2])[,1] %in% keepy)

    bub <- bub[ix,]

  }


  if(highlight) {

    plt <- ungroup(bub) %>%
      mutate(status = case_when(POTH == max(POTH) ~ "Largest",
                                POTH == min(POTH) ~ "Smallest",
                                TRUE ~ "other"))

  } else {

    plt <- bub %>% mutate(status = "other")

  }


  g <-ggplot(plt, aes(x = !!sym(names(x$CIVs)[1]),
                  y = !!sym(names(x$CIVs)[2]),
                  size = POTH, col = status, fill = status)) +
    geom_point(alpha = 0.7, shape = 21) +
    scale_size(range = c(1, 20)) +
    scale_fill_manual(breaks = c("Largest", "Smallest"),
                       values = c(Largest = "lightskyblue3", Smallest = "hotpink3", other = "lightyellow")) +
    scale_color_manual(breaks = c("Largest", "Smallest"),
                       values = c(Largest = "lightskyblue4", Smallest = "hotpink4", other = "black")) +
    guides(size = guide_legend(override.aes = list(fill = "lightyellow"), order = 1),
           fill = guide_legend(override.aes = list(size = 10),
                               title = "Useful POTH Values"),
           colour = guide_legend(title = "Useful POTH Values")) +
    coord_cartesian(clip = "off") +
    labs(x = paste0("CIV (", names(x$CIVs)[1], ")"),
         y = paste0("CIV (", names(x$CIVs)[2], ")"),
         title = title) +
    theme_bw()


  print(g)

  return(list(plot = g, grid = bub))


}

#' @param x1 data.frame with column names `Treatment` and `ranking`
#' @param x2 data.frame with column names `Treatment` and `ranking`
#' @param name1 Name of the first group of rankings (x axis label)
#' @param name2 Name of the second group of rankings (y axis label)
pscores_compare <- function(x1, x2, name1, name2) {

  master <- left_join(x1, x2, join_by(Treatment == Treatment))

  g <-ggplot(master, aes(x = ranking.x, y = ranking.y)) +
    # geom_point(size = 4, col = "black", shape = 22) +
    geom_abline(slope = 1, intercept = 0) +
    geom_label(aes(label = Treatment), fill = "lightyellow") +
    guides(fill = guide_legend(override.aes = aes(label = "  "))) +
    labs(x = name1, y = name2) +
    # scale_fill_viridis_d(option = "turbo", begin = 0.05)+
    theme_bw()

  print(g)


}

# Helpers -----------------
poth2 <- function(pscores) {

  n <- length(pscores)

  sq <- pscores - mean(pscores)
  p <- sum(sq*sq)/n*12*(n-1)/(n+1)

  return(p)

}

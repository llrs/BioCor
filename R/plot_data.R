
# Function to make it easier to plot the data:
#' The position of the nodes is based on the similarity between them.
#' @param x Matrix with the similarities.
#'
#' @param top  a number between 0 and 1 to select the edges relating the elements of the matrix.
#' @returns A list with two elements:
#' - nodes: The position and name of the nodes
#' - edges: The information about the selected edges
#' @export
#' @rdname plot_similarity
plot_data <- function(x, top) {

  # Pick up the edges to show
  # 1: From a similarity matrix calculate the density of the edges.
  # 2: From there calculate the quantile to pick the top x% (parameter)
  # 3: Transform it into a data.frame
  stopifnot(isSymmetric(x))
  stopifnot(top < 1, top > 0)
  limit <- quantile(x[upper.tri(x)], probs = 1 - top)
  positions <- which(x >= limit & upper.tri(x), arr.ind = TRUE)
  df <- data.frame(
    A = colnames(x)[positions[, 1, drop = FALSE]],
    B = colnames(x)[positions[, 2, drop = FALSE]])
  df[["strength"]] <- x[positions]
  df[["rank"]] <- rank(x[positions])

  # Pick the distribution of the points
  # 1: calculate the dissimilarity
  # 2: calculate the PCA of the first 2 dimensions
  MDS <- as.data.frame(cmdscale(1 - x))
  colnames(MDS) <- c("x", "y")
  MDS[["node"]] <- rownames(MDS)

  # Prepare for plotting data
  A <- merge(df, MDS, by.x = "A", by.y = "node")
  B <- merge(df, MDS, by.x = "B", by.y = "node")
  m <- merge(A, B, by = c("B", "A", "strength", "rank"), suffixes = c(".start", ".end"))

  pd <- list(nodes = MDS, edges = m)
  return(pd)
}

# Plot the entities and its edges


# Function to plot the data with ggplot2 if installed.
#' Plot how similar are the data
#' @param pd The plot data from `plot_data()` function.
#' @returns A ggplot object
#' @importFrom stats cmdscale
#' @importFrom stats quantile
#' @export
#' @examples
#' if (require("org.Hs.eg.db") & require("reactome.db")) {
#'   # Extract the paths of all genes of org.Hs.eg.db from KEGG
#'   # (last update in data of June 31st 2011)
#'   genes.kegg <- as.list(org.Hs.egPATH)
#'   # Extracts the paths of all genes of org.Hs.eg.db from reactome
#'   genes.react <- as.list(reactomeEXTID2PATHID)
#'
#'   sim <- mgeneSim(c("81", "18", "10"), genes.react)
#'   pd <- plot_data(sim, top = 0.25)
#'   if (requireNamespace("ggplot2", quietly = TRUE)){
#'     plot_similarity(pd)
#'   }
#' }
plot_similarity <- function(pd) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Please install the ggplot2 package.")
  }
  .data <- NULL # Trick to avoid check notes
  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(data = pd[["edges"]],
                          ggplot2::aes(x = .data$x.start, y = .data$y.start,
                                       xend = .data$x.end, yend = .data$y.end,
                                       linewidth = .data$strength)) +
    ggplot2::geom_label(data = pd$nodes,
                        ggplot2::aes(.data$x, .data$y, label = .data$node), fill = "white") +
    ggplot2::theme_void()
  return(p)
}

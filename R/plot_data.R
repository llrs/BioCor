
# Function to make it easier to plot the data:
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
  df$strength <- x[positions]
  df$rank <- rank(df$strength)
  
  # Pick the distribution of the points
  # 1: calculate the dissimilarity
  # 2: calculate the PCA of the first 2 dimensions
  MDS <- as.data.frame(cmdscale(1 - x))
  colnames(MDS) <- c("x", "y")
  MDS$node <- rownames(MDS)
  pd <- list(node = MDS, edges = df)
  
  
  if (!requireNamespace("ggraph", quietly = TRUE)) {
    warning("tidygraph is not present: you'll need to work out how to represent the output")
    return(pd)
  }
  
  edges <- pd$edges
  colnames(edges)[1:2] <- c("from", "to")
  node <- pd$node
  colnames(node)[3] <- "name"
  node <- node[, c(1, 2:3)]
  
  tbl_graph(nodes = node, edges = edges, directed = FALSE)
}

# Plot the entities and its edges


# Function to plot the data with ggplot2 if installed.
plot_similarity <- function(pd) {
  if (!requireNamespace("tidygraph", quietly = TRUE)) {
    stop("Please install the ggraph package.")
  }

  ggraph(graph, layout = "manual", x = x, y = y) +
    geom_edge_link() + 
    geom_node_label(aes(label = name)) +
    theme_graph()
    
}
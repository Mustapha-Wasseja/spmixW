#' Row-Normalise a Spatial Weight Matrix
#'
#' Divides each row of a spatial weight matrix by its row sum so that
#' all rows sum to one. Zero-sum rows (isolates) are left unchanged.
#'
#' @param W A square matrix or sparse matrix (class \code{dgCMatrix} or similar).
#'
#' @return A sparse matrix of the same dimension with rows summing to 1
#'   (except for isolate rows that remain zero).
#'
#' @examples
#' W <- matrix(c(0, 1, 1, 0,
#'               1, 0, 0, 1,
#'               1, 0, 0, 1,
#'               0, 1, 1, 0), 4, 4)
#' Wn <- normw(W)
#' Matrix::rowSums(Wn)  # all ones
#'
#' @export
normw <- function(W) {

  W <- Matrix::Matrix(W, sparse = TRUE)
  rs <- Matrix::rowSums(W)

  # Avoid division by zero for isolates
  rs[rs == 0] <- 1

  # Efficient row-scaling: D^{-1} W where D = diag(rs)
  Dinv <- Matrix::Diagonal(x = 1 / rs)
  Dinv %*% W
}


#' Create k-Nearest-Neighbours Spatial Weight Matrix
#'
#' Constructs a row-normalised binary spatial weight matrix based on
#' k-nearest neighbours from a coordinate matrix.
#'
#' @param coords An \eqn{N \times 2} matrix of coordinates (e.g., longitude/latitude
#'   or projected x/y).
#' @param k Integer: number of nearest neighbours.
#' @param row_normalise Logical: if \code{TRUE} (default), row-normalise the
#'   resulting matrix.
#'
#' @return A sparse \eqn{N \times N} weight matrix. If \code{row_normalise = TRUE},
#'   each row sums to 1.
#'
#' @details
#' Euclidean distances are computed between all pairs of coordinates. For each
#' unit, the \code{k} closest neighbours receive weight 1 (before normalisation).
#' The matrix is generally asymmetric (i may be j's neighbour without j being
#' i's neighbour).
#'
#' @examples
#' set.seed(1)
#' coords <- cbind(runif(20), runif(20))
#' W <- make_knw(coords, k = 5)
#' dim(W)
#' range(Matrix::rowSums(W))  # all 1 if row-normalised
#'
#' @export
make_knw <- function(coords, k, row_normalise = TRUE) {

  coords <- as.matrix(coords)
  N <- nrow(coords)
  stopifnot(k < N, k >= 1L, ncol(coords) >= 2)

  # Compute pairwise Euclidean distance matrix
  d <- as.matrix(dist(coords))

  # For each row, find the indices of the k smallest distances (excluding self)
  diag(d) <- Inf  # exclude self-distance

  # Build sparse triplets
  i_idx <- integer(N * k)
  j_idx <- integer(N * k)

  for (row in seq_len(N)) {
    nn <- order(d[row, ])[seq_len(k)]
    idx <- ((row - 1L) * k + 1L):(row * k)
    i_idx[idx] <- row
    j_idx[idx] <- nn
  }

  W <- Matrix::sparseMatrix(
    i = i_idx, j = j_idx, x = rep(1, N * k),
    dims = c(N, N)
  )

  if (row_normalise) {
    W <- normw(W)
  }

  W
}

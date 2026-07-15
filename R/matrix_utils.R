#' Matrix power
#'
#' Computes \code{A^k} for a square matrix \code{A} and a non-negative integer
#' \code{k} using exponentiation by squaring (O(log k) matrix products instead
#' of O(k) in the original implementation).
#'
#' @param A A square numeric matrix.
#' @param k A non-negative integer.
#'
#' @return The matrix \code{A} raised to the power \code{k}.
#' @keywords internal
matpow <- function(A, k) {
  stopifnot(is.matrix(A), nrow(A) == ncol(A), length(k) == 1L,
            is.finite(k), k >= 0, k == round(k))
  n <- nrow(A)
  R <- diag(n)
  if (k == 0) return(R)
  B <- A
  while (k > 0) {
    if (k %% 2 == 1) R <- R %*% B
    k <- k %/% 2
    if (k > 0) B <- B %*% B
  }
  R
}

#' Stationary covariance of a stable linear system
#'
#' Solves the discrete-time Lyapunov equation \code{P = F P F' + Q} for the
#' stationary state covariance \code{P} of the system
#' \code{s_t = F s_{t-1} + eta_t}, \code{Var(eta_t) = Q}.
#'
#' @details
#' Two methods are available:
#' \describe{
#'   \item{\code{"kron"} (default)}{Exact solution of the linear system
#'     \code{(I - F \%x\% F) vec(P) = vec(Q)}. For state dimension d this is a
#'     d^2 x d^2 solve, which is fast for the moderate dimensions used here
#'     (d = 2K) and, unlike the fixed-point iteration, does not slow down as
#'     the spectral radius of \code{F} approaches 1.}
#'   \item{\code{"iterate"}}{The fixed-point iteration
#'     \code{P <- F P F' + Q} used in the original script. Convergence is
#'     linear with rate equal to the squared spectral radius of \code{F}, so
#'     this can require thousands of iterations for highly persistent systems.
#'     Kept for reference and cross-checking.}
#' }
#' The returned matrix is symmetrized. A residual check is performed for the
#' \code{"kron"} method; \code{NULL} is returned if the equation could not be
#' solved reliably (e.g., because \code{F} is not stable).
#'
#' @param Fmat Square transition matrix of the system.
#' @param Qmat Innovation covariance matrix (same dimension as \code{Fmat}).
#' @param method Either \code{"kron"} (exact, default) or \code{"iterate"}.
#' @param tol Convergence tolerance for \code{method = "iterate"}.
#' @param maxit Maximum number of iterations for \code{method = "iterate"}.
#'
#' @return The stationary covariance matrix \code{P}, or \code{NULL} if no
#'   reliable solution was found.
#'
#' @examples
#' m <- make_mats_varK(ax = c(0.5, 0.1), ay = c(0.4, 0.05),
#'                     b21 = c(0.2, 0.05), b12 = c(0.1, 0.02))
#' P <- stationary_cov(m$F, m$Q)
#'
#' @export
stationary_cov <- function(Fmat, Qmat, method = c("kron", "iterate"),
                           tol = 1e-12, maxit = 400000L) {
  method <- match.arg(method)
  stopifnot(is.matrix(Fmat), is.matrix(Qmat),
            nrow(Fmat) == ncol(Fmat),
            all(dim(Fmat) == dim(Qmat)),
            all(is.finite(Fmat)), all(is.finite(Qmat)))
  d <- nrow(Fmat)

  # A finite, positive semi-definite stationary covariance exists only if the
  # system is stable, i.e. every eigenvalue of Fmat has modulus < 1. Without
  # this guard the Kronecker solve can return a spurious algebraic solution:
  # for an unstable Fmat the equation P = F P F' + Q still has a solution
  # (whenever no eigenvalue product equals 1), but that solution is not a valid
  # covariance and the residual check does not catch it. The check costs one
  # d x d eigendecomposition, negligible next to the d^2 x d^2 Kronecker solve;
  # it also short-circuits the "iterate" method for unstable systems instead of
  # letting it diverge to Inf.
  eig <- eigen(Fmat, only.values = TRUE)$values
  if (!all(is.finite(eig)) || max(Mod(eig)) >= 1) return(NULL)

  if (method == "kron") {
    A <- diag(d * d) - kronecker(Fmat, Fmat)
    vecP <- tryCatch(solve(A, as.vector(Qmat)), error = function(e) NULL)
    if (is.null(vecP) || !all(is.finite(vecP))) return(NULL)
    P <- matrix(vecP, d, d)
    P <- (P + t(P)) / 2
    resid <- max(abs(P - Fmat %*% P %*% t(Fmat) - Qmat))
    ref <- max(1, max(abs(Qmat)), max(abs(P)))
    if (!is.finite(resid) || resid > 1e-6 * ref) return(NULL)
    return(P)
  }

  # method == "iterate": original fixed-point iteration, kept for reference
  P <- Qmat
  for (i in seq_len(maxit)) {
    P_new <- Fmat %*% P %*% t(Fmat) + Qmat
    if (!all(is.finite(P_new))) return(NULL)
    if (max(abs(P_new - P)) < tol) return((P_new + t(P_new)) / 2)
    P <- P_new
  }
  NULL
}

#' Cholesky factor for (possibly borderline) positive semi-definite matrices
#'
#' Internal helper. Symmetrizes the input and, if a plain Cholesky
#' factorization fails, retries with a small diagonal jitter. Used for drawing
#' multivariate normal variates without depending on MASS.
#'
#' @param M A symmetric positive (semi-)definite matrix.
#' @param jitter Initial jitter added to the diagonal on failure.
#'
#' @return Upper-triangular matrix \code{R} with \code{t(R) \%*\% R}
#'   (approximately) equal to \code{M}.
#' @keywords internal
chol_psd <- function(M, jitter = 1e-10) {
  M <- (M + t(M)) / 2
  out <- tryCatch(chol(M), error = function(e) NULL)
  if (!is.null(out)) return(out)
  d <- nrow(M)
  scale <- max(abs(diag(M)), 1e-12)
  for (j in c(jitter, 1e-8, 1e-6)) {
    out <- tryCatch(chol(M + diag(j * scale, d)), error = function(e) NULL)
    if (!is.null(out)) return(out)
  }
  stop("Cholesky factorization failed; covariance matrix is not positive semi-definite.")
}

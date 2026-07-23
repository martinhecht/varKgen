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
#' \code{s_t = F s_{t-1} + eta_t}, where \code{Q = Var(eta_t)} is the process
#' error covariance.
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
#' @param Qmat Process error covariance matrix (same dimension as \code{Fmat}).
#' @param method Either \code{"kron"} (exact, default) or \code{"iterate"}.
#' @param tol Convergence tolerance for \code{method = "iterate"}.
#' @param maxit Maximum number of iterations for \code{method = "iterate"}.
#' @param check_stability If \code{TRUE} (default), the spectral radius of
#'   \code{Fmat} is verified before solving. Set to \code{FALSE} only when
#'   the caller has already established stability, to avoid a second
#'   eigendecomposition per evaluation.
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
                           tol = 1e-12, maxit = 400000L,
                           check_stability = TRUE) {
  method <- match.arg(method)
  stopifnot(is.matrix(Fmat), is.matrix(Qmat),
            nrow(Fmat) == ncol(Fmat),
            all(dim(Fmat) == dim(Qmat)),
            all(is.finite(Fmat)), all(is.finite(Qmat)),
            length(tol) == 1L, is.finite(tol), tol > 0,
            length(maxit) == 1L, is.finite(maxit), maxit >= 1L,
            maxit == round(maxit),
            length(check_stability) == 1L, is.logical(check_stability),
            !is.na(check_stability))
  # Qmat must be a covariance matrix; symmetry is checked here, positive
  # semi-definiteness remains the caller's responsibility (an eigen-
  # decomposition per call would be too costly inside the loss function)
  q_scale <- max(1e-300, max(abs(Qmat)))
  if (max(abs(Qmat - t(Qmat))) > 1e-10 * q_scale) {
    stop("'Qmat' is not symmetric.", call. = FALSE)
  }
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
  if (check_stability) {
    eig <- eigen(Fmat, only.values = TRUE)$values
    if (!all(is.finite(eig)) || max(Mod(eig)) >= 1) return(NULL)
  }

  if (method == "kron") {
    A <- diag(d * d) - kronecker(Fmat, Fmat)
    vecP <- tryCatch(solve(A, as.vector(Qmat)), error = function(e) NULL)
    if (is.null(vecP) || !all(is.finite(vecP))) return(NULL)
    P <- matrix(vecP, d, d)
    P <- (P + t(P)) / 2
    resid <- max(abs(P - Fmat %*% P %*% t(Fmat) - Qmat))
    # relative reference with a machine-level absolute floor, so the check
    # stays meaningful for covariance scales far below 1
    ref <- max(.Machine$double.eps, max(abs(Qmat)), max(abs(P)))
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

#' Process error covariance for a target stationary process covariance
#'
#' Internal. Given a stable companion transition matrix \code{Fmat}, finds the
#' 2 x 2 process error covariance \code{Q2} such that the stationary state
#' covariance \code{P} of the system \code{s_t = F s_{t-1} + eta_t} (with
#' \code{Var(eta_t)} equal to \code{Q2} embedded in the top-left block and
#' zeros elsewhere) has \code{P[1, 1] = var1}, \code{P[2, 2] = var2} and
#' \code{P[1, 2] = cov12}.
#'
#' @details
#' The stationary covariance is linear in the process error covariance, so the
#' three unknowns (both error variances and the error covariance) solve a
#' 3 x 3 linear system built from the Lyapunov responses to the three basis
#' components of \code{Q2}. One LU factorization of
#' \code{I - F \%x\% F} with three right-hand sides yields the responses; the
#' full stationary covariance follows by linearity without a further solve.
#' The achieved target block is verified to within a small relative tolerance.
#'
#' The derived \code{Q2} is not guaranteed to be positive semi-definite: not
#' every target process covariance is admissible for given dynamics. The
#' smallest eigenvalue of \code{Q2} (closed form for the symmetric 2 x 2 case)
#' is returned so that callers can penalize or reject inadmissible points.
#'
#' @param Fmat Stable companion transition matrix (2K x 2K, K >= 1).
#' @param var1 Target stationary variance of process 1 (\code{P[1, 1]}).
#' @param var2 Target stationary variance of process 2 (\code{P[2, 2]}).
#' @param cov12 Target stationary covariance (\code{P[1, 2]}).
#' @param check_stability If \code{TRUE} (default), the spectral radius of
#'   \code{Fmat} is verified before solving. Set to \code{FALSE} only when
#'   the caller has already established stability, to avoid a second
#'   eigendecomposition per evaluation.
#'
#' @return \code{NULL} if \code{Fmat} is unstable or the linear solves fail;
#'   otherwise a list with \code{Q2} (2 x 2 derived process error
#'   covariance), \code{P} (full stationary state covariance, symmetrized)
#'   and \code{min_eig} (smallest eigenvalue of \code{Q2}; negative values
#'   indicate an inadmissible target).
#' @keywords internal
solve_process_error <- function(Fmat, var1, var2, cov12,
                                check_stability = TRUE) {
  stopifnot(is.matrix(Fmat), nrow(Fmat) == ncol(Fmat), nrow(Fmat) >= 2L,
            all(is.finite(Fmat)),
            length(var1) == 1L, is.finite(var1),
            length(var2) == 1L, is.finite(var2),
            length(cov12) == 1L, is.finite(cov12),
            length(check_stability) == 1L, is.logical(check_stability),
            !is.na(check_stability))
  d <- nrow(Fmat)

  if (check_stability) {
    eig <- eigen(Fmat, only.values = TRUE)$values
    if (!all(is.finite(eig)) || max(Mod(eig)) >= 1) return(NULL)
  }

  A <- diag(d * d) - kronecker(Fmat, Fmat)

  # basis right-hand sides: vec of the three symmetric basis components of the
  # top-left 2 x 2 block (column-major vec, as in as.vector())
  B <- matrix(0, d * d, 3L)
  E <- matrix(0, d, d); E[1L, 1L] <- 1
  B[, 1L] <- as.vector(E)
  E <- matrix(0, d, d); E[2L, 2L] <- 1
  B[, 2L] <- as.vector(E)
  E <- matrix(0, d, d); E[1L, 2L] <- 1; E[2L, 1L] <- 1
  B[, 3L] <- as.vector(E)

  V <- tryCatch(solve(A, B), error = function(e) NULL)
  if (is.null(V) || !all(is.finite(V))) return(NULL)

  # rows of vec(P) corresponding to P[1,1], P[2,2], P[1,2] (column-major)
  i11 <- 1L
  i22 <- d + 2L
  i12 <- d + 1L
  M <- V[c(i11, i22, i12), , drop = FALSE]

  q <- tryCatch(solve(M, c(var1, var2, cov12)), error = function(e) NULL)
  if (is.null(q) || !all(is.finite(q))) return(NULL)

  vecP <- as.vector(V %*% q)
  P <- matrix(vecP, d, d)
  P <- (P + t(P)) / 2

  Q2 <- matrix(c(q[1L], q[3L], q[3L], q[2L]), 2L, 2L)

  # verify the achieved target block (relative, with a machine-level floor)
  achieved <- c(P[1L, 1L], P[2L, 2L], P[1L, 2L])
  ref <- max(.Machine$double.eps, abs(var1), abs(var2), abs(cov12),
             max(abs(P)))
  if (max(abs(achieved - c(var1, var2, cov12))) > 1e-6 * ref) return(NULL)

  # the target-block check alone can miss an inadequate full solution for
  # ill-conditioned systems, so verify the complete Lyapunov equation
  Qfull <- matrix(0, d, d)
  Qfull[1:2, 1:2] <- Q2
  resid <- max(abs(P - Fmat %*% P %*% t(Fmat) - Qfull))
  ref_l <- max(.Machine$double.eps, max(abs(Qfull)), max(abs(P)))
  if (!is.finite(resid) || resid > 1e-6 * ref_l) return(NULL)

  # smallest eigenvalue of the symmetric 2 x 2 Q2; the closed form can
  # overflow for extreme entries, so use the symmetric eigensolver
  ev_q <- eigen(Q2, symmetric = TRUE, only.values = TRUE)$values
  if (!all(is.finite(ev_q))) return(NULL)
  min_eig <- min(ev_q)

  list(Q2 = Q2, P = P, min_eig = min_eig)
}

#' Symmetric factor of a positive semi-definite covariance matrix
#'
#' Internal helper. Returns a matrix \code{L} with
#' \code{t(L) \%*\% L = M_used}, so that \code{Z \%*\% L} has covariance
#' \code{M_used} for \code{Z} with independent standard normal entries.
#'
#' @details
#' A Cholesky factorization is not used because \code{chol()} fails for
#' exactly singular positive semi-definite matrices, and repairing this with
#' a diagonal jitter would silently simulate from \code{M + epsilon * I}
#' instead of \code{M}. Singular covariance matrices are reachable here:
#' \code{\link{make_mats_varK}} admits
#' \code{abs(pe_cov12) = sqrt(pe_var1 * pe_var2)}, and calibrated process
#' error covariances can be numerically near-singular.
#'
#' The symmetric eigendecomposition is used instead. \code{M_used} is the
#' symmetrized matrix \code{(M + t(M)) / 2} with negative eigenvalues
#' truncated to zero. Negative eigenvalues below \code{-tol * scale}, where
#' \code{scale} is the largest absolute eigenvalue, are rejected with an
#' error; only negative eigenvalues at the level of numerical noise are
#' truncated. The tolerance is therefore purely relative, and for a matrix
#' that is positive semi-definite up to roundoff the factorization is exact
#' to within that roundoff. The applied correction, i.e. the magnitude of
#' the most negative truncated eigenvalue, is attached as the attribute
#' \code{"psd_correction"} so that callers can disclose it; it is zero
#' whenever no truncation was needed. The returned factor is not triangular,
#' which is irrelevant here: only the covariance decomposition matters.
#'
#' @param M Symmetric covariance matrix.
#' @param tol Relative tolerance for negative eigenvalues, measured against
#'   the largest absolute eigenvalue.
#'
#' @return A matrix \code{L} with \code{t(L) \%*\% L = M_used}, carrying
#'   the attribute \code{"psd_correction"}.
#' @keywords internal
cov_factor_psd <- function(M, tol = 1e-10) {
  stopifnot(is.matrix(M), nrow(M) == ncol(M), nrow(M) >= 1L,
            all(is.finite(M)),
            length(tol) == 1L, is.numeric(tol), is.finite(tol), tol >= 0)
  M <- (M + t(M)) / 2
  ee <- eigen(M, symmetric = TRUE)
  if (!all(is.finite(ee$values))) {
    stop("Eigendecomposition of the covariance matrix failed.", call. = FALSE)
  }
  scale <- max(abs(ee$values))
  if (scale == 0) {
    L <- matrix(0, nrow(M), ncol(M))
    attr(L, "psd_correction") <- 0
    return(L)
  }
  if (min(ee$values) < -tol * scale) {
    stop("Covariance matrix is not positive semi-definite (smallest ",
         "eigenvalue ", format(min(ee$values), digits = 3),
         ", relative to scale ", format(scale, digits = 3), ").",
         call. = FALSE)
  }
  values <- pmax(ee$values, 0)
  L <- diag(sqrt(values), nrow = length(values)) %*% t(ee$vectors)
  attr(L, "psd_correction") <- max(0, -min(ee$values))
  L
}

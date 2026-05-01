# Internal: add standardised field name aliases to result objects
# Called at the end of every estimation function before returning
#
# Old name -> New alias:
#   pdraw    -> rho_draws (numeric vector)
#   bdraw    -> beta_draws (matrix)
#   sdraw    -> sigma_draws (numeric vector)
#   gdraw    -> gamma_draws (matrix)
#   sige     -> sigma2 (scalar)
#   rho      -> lambda (for SEM/SDEM types, keeps rho too)
#
# All old names remain for backward compatibility. New names are aliases.

.add_field_aliases <- function(result) {

  # Scalar aliases
  if (!is.null(result$sige)) result$sigma2 <- result$sige

  # SEM/SDEM: add lambda alias for rho
  if (!is.null(result$rho) && .is_sem_type(result$meth %||% "")) {
    result$lambda <- result$rho
  }

  # Draw aliases (as plain numeric vectors/matrices, not coda objects)
  if (!is.null(result$pdraw)) {
    result$rho_draws <- as.numeric(result$pdraw)
    if (.is_sem_type(result$meth %||% "")) {
      result$lambda_draws <- result$rho_draws
    }
  }
  if (!is.null(result$bdraw)) {
    result$beta_draws <- as.matrix(result$bdraw)
  }
  if (!is.null(result$sdraw)) {
    result$sigma_draws <- as.numeric(result$sdraw)
  }
  if (!is.null(result$gdraw)) {
    result$gamma_draws <- as.matrix(result$gdraw)
  }

  result
}

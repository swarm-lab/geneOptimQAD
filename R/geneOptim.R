# Generic genetic optimization function.
#
# Args:
#   fn: Function to minimize. The ouput of the function should be a single number.
#   ...: Parameters for fn.
#   N: Number of individual in the population (default: 100).
#   T: Number of "generations". A generation actually corresponds to a reproductive event between the best two individuals of a group of four individuals randomly chosen in the population (default: 500).
#   minVal: Vector of the minimum acceptable values for each gene.
#   maxVal: Vector of the maximum acceptable values for each gene.
#   pmut: Probability that a gene mutate (default: 0.25).
#   wmut: Maximum variation allowed through mutation (default: 0.25).
#   pcross: Probability of crossover (default: 0.1).
#   ncores: Number of available cores for parallel computation (default: 1).
#
# Returns:
#   A matrix of the final population.

#' Quick & dirty genetic algorithm
#'
#' Maximization of a function using a quick & dirty algorithm.
#'
#' @param fn Function to maximize. The ouput of the function should be a single
#'  numeric value.
#'
#' @param N Number of individual in the population (default: 100).
#'
#' @param G Number of "generations". A generation actually corresponds to a
#'  reproductive event between the best two individuals of a group of four
#'  individuals randomly chosen in the population (default: 500).
#'
#' @param min_val A vector of the minimum acceptable values for each gene.
#'
#' @param max_val A vector of the maximum acceptable values for each gene.
#'
#' @param p_mut Probability of a gene mutating (default: 0.25).
#'
#' @param w_mut Maximum variation allowed through mutation (default: 0.25).
#'
#' @param p_cross Probability of crossover (default: 0.1).
#'
#' @param n_cores Number of available cores for parallel computation (default: 1).
#'
#' @return A matrix of the final population.
#'
#' @seealso \code{\link{}}
#'
#' @examples
#' #TODO
geneOptim <- function(fn, ..., N = 100, G = 500, min_val, max_val,
                      p_mut = 0.25, w_mut = 0.2, p_cross = 0.1,
                      n_cores = 1, .packages = NULL) {
  require(foreach)

  if (n_cores > 1) {
    if (n_cores > 4)
      n_cores <- 4

    cl <- parallel::makeCluster(n_cores)
    doSNOW::registerDoSNOW(cl)
  }

  pop <- initPop(N, min_val, max_val)
  k <- length(min_val)
  min_fit <- {}
  mean_fit <- {}

  for (i in 1:G) {
    print(paste0("Iteration n°", i))

    # Choose 4 individuals at random and evaluate their fitness
    indices <- sample(N, 4)
    individuals <- matrix(0, nrow = 4, ncol = 2)

    foreach(j = 1:4, .packages = .packages) %dopar% {
      index <- indices[j]
      error <- fn(as.numeric(pop[index, 1:k]), ...)
      pop[index, k + 1] <- error
      individuals[j, 1] <- error
      individuals[j, 2] <- index
    }

    # Sort individuals
    individuals <- individuals[order(individuals[, 1]), ]
    indexParent1 <- individuals[1, 2]
    indexParent2 <- individuals[2, 2]
    indexChild1 <- individuals[3, 2]
    indexChild2 <- individuals[4, 2]

    # Generate two new individuals and replace the two worst ones
    pop[indexChild1, ] <- createOffspring(as.numeric(pop[indexParent1, ]), as.numeric(pop[indexParent2, ]), min_val, max_val)
    pop[indexChild2, ] <- createOffspring(as.numeric(pop[indexParent1, ]), as.numeric(pop[indexParent2, ]), min_val, max_val)
    min_fit[i] <- min(pop[, k + 1], na.rm = TRUE)
    mean_fit[i] <- mean(pop[, k + 1], na.rm = TRUE)
    pop <- pop[order(pop[, k + 1]), ]

    # Display the current best parameters
    progress(pop, k, min_fit, mean_fit)
    Sys.sleep(0)
  }

  if (n_cores > 1)
    parallel::stopCluster(cl)

  pop
}





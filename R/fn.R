#' Create initial population
#'
#' This function creates an initial population of randomly generated individuals.
#'
#' @param N Number of individuals in the population.
#'
#' @param min_val A vector of the minimum acceptable values for each gene.
#'
#' @param max_val A vector of the maximum acceptable values for each gene.
#'
#' @return A N x length(min_val) matrix. Each line corresponds to an individual
#'  in the initial population. Each column correspond to the value of each gene,
#'  randomly chosen between min_val and max_val. The last column corresponds to
#'  the fitness of each individual (NA after initialization).
#'
#' @seealso \code{\link{}}
#'
#' @examples
#' #TODO
initPop <- function(N, min_val, max_val) {
  if (length(min_val) != length(max_val) | !is.vector(min_val) | !is.vector(max_val) | !is.numeric(min_val) | !is.numeric(max_val))
    stop("min_val and max_val must be numeric vectors of the same length.")

  k <- length(min_val)
  pop <- matrix(0, nrow = N, ncol = k + 1)

  pop[, k + 1] <- NA
  pop[, 1:k] <- t(replicate(N, min_val + (max_val - min_val) * runif(k)))

  pop <- as.data.frame(pop)
  names(pop) <- c(paste0("Gene", 1:k), "Fitness")
  pop
}


#' Create new offspring
#'
#' Generate a new offspring from the genome of two parents.
#'
#' @param genome1 Vector containing the gene values for individual 1.
#'
#' @param genome2 Vector containing the gene values for individual 2.
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
#' @return A vector containing the gene values of the new individual.
#'
#' @seealso \code{\link{}}
#'
#' @examples
#' #TODO
createOffspring <- function(genome1, genome2, min_val, max_val, p_mut = 0.25, w_mut = 0.2, p_cross = 0.1) {
  k <- length(min_val)
  individual <- matrix(0, nrow = 1, ncol = k + 1)
  individual[1, k + 1] <- NA

  # Randomize genomes from parents
  individual[1, 1:k] <- ifelse(runif(k) < 0.5, genome1, genome2)

  # Mutate
  idx <- which(runif(k) < p_mut)
  individual[1, idx] <- individual[1, idx] * (1 - w_mut / 2) + w_mut * runif(length(idx))

  # Crossover
  idx <- which(runif(k) < p_cross)
  individual[1, idx] <- min_val[idx] + (max_val[idx] - min_val[idx]) * runif(length(idx))

  # Check consistency
  idx <- which(individual[1, 1:k] < min_val)
  individual[1, idx] <- min_val[idx]

  idx <- which(individual[1, 1:k] > max_val)
  individual[1, idx] <- max_val[idx]

  individual
}


#' Track optimization progress
#'
#' Display current state and history of the genetic optimization.
#'
#' @param pop Matrix of the gene values of the current population. Each line
#'  correspond to an individual in the initial population. Each column correspond
#'  to the value of each gene. The last column corresponds to the fitness of each
#'  individual.
#'
#' @param k Number of genes.
#'
#' @param min_fit Vector of the minimum fitness for each generation.
#'
#' @param mean_fit Vector of the average fitness for each generation.
#'
#' @return Plot the minimum and average fitness for the last 100 generations.
#'
#' @seealso \code{\link{}}
#'
#' @examples
#' #TODO
progress <- function(pop, k, min_fit, mean_fit) {
  print("Current best parameters:")
  print(paste(names(pop)[1:k], sep = "	"))
  print(paste(pop[1, 1:k], sep = "	"))
  print(paste("Min fit =", min_fit[length(min_fit)], "- Mean fit =", mean_fit[length(mean_fit)]))

  if (length(mean_fit) > 100) {
    plot(mean_fit[(length(mean_fit) - 100):length(mean_fit)] ~ c((length(mean_fit) - 100):length(mean_fit)),
         type = "l", ylim = c(0, max(mean_fit[(length(mean_fit) - 100):length(mean_fit)])))
    lines(min_fit[(length(min_fit) - 100):length(min_fit)] ~ c((length(min_fit) - 100):length(min_fit)), col = "red")
  } else {
    plot(mean_fit ~ c(1:length(mean_fit)), ylim = c(0, max(mean_fit)), type = "l")
    lines(min_fit ~ c(1:length(min_fit)), col = "red")
  }
}






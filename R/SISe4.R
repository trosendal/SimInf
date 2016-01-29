## siminf, a framework for stochastic disease spread simulations
## Copyright (C) 2015  Pavol Bauer
## Copyright (C) 2015  Stefan Engblom
## Copyright (C) 2015  Stefan Widgren
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.

##' Class \code{"SISe4"}
##'
##' Class to handle the SISe4 \code{\link{siminf_model}} model.
##' @include siminf_model.R
##' @include AllGenerics.R
##' @export
setClass("SISe4", contains = c("siminf_model"))

##' Create a SISe4 model
##'
##' Create a SISe4 model to be used by the simulation framework.
##'
##'
##' The argument init must be a \code{data.frame} with the following
##' columns:
##' \describe{
##' \item{id}{Node identifier that uniquely identifies each node. The
##' node identifiers must be zero-based, i.e. the first identifier
##' must be equal to zero.}
##' \item{S_1}{The number of sucsceptible in age category 1}
##' \item{I_1}{The number of infected in age category 1}
##' \item{S_2}{The number of sucsceptible in age category 2}
##' \item{I_2}{The number of infected in age category 2}
##' \item{S_3}{The number of sucsceptible in age category 3}
##' \item{I_3}{The number of infected in age category 3}
##' \item{S_3}{The number of sucsceptible in age category 4}
##' \item{I_3}{The number of infected in age category 4}
##' }
##'
##' @template beta-section
##' @param init A \code{data.frame} with the initial state in each
##' node, see details.
##' @param tspan An increasing sequence of points in time where the
##' state of the system is to be returned.
##' @param events a \code{data.frame} with the scheduled events, see
##' \code{\link{siminf_model}}.
##' @param phi A numeric vector with the initial environmental
##' infectious pressure in each node. Default NULL which gives 0 in
##' each node.
##' @param upsilon_1 Indirect transmission rate of the environmental
##' infectious pressure in age category 1
##' @param upsilon_2 Indirect transmission rate of the environmental
##' infectious pressure in age category 2
##' @param upsilon_3 Indirect transmission rate of the environmental
##' infectious pressure in age category 3
##' @param upsilon_4 Indirect transmission rate of the environmental
##' infectious pressure in age category 3
##' @param gamma_1 The recovery rate from infected to susceptible for
##' age category 1
##' @param gamma_2 The recovery rate from infected to susceptible for
##' age category 2
##' @param gamma_3 The recovery rate from infected to susceptible for
##' age category 3
##' @param gamma_4 The recovery rate from infected to susceptible for
##' age category 3 
##' @param alpha Shed rate from infected individuals
##' @param epsilon The background environmental infectious pressure
##' @param lambda_sd1 the introduction rate for subdomain1 
##' @param lambda_sd2 the introduction rate for subdomain2
##' @param lambda_sd3 the introduction rate for subdomain3
##' @param lambda_sd4 the introduction rate for subdomain4
##' @param lambda_sd5 the introduction rate for subdomain5
##' @param lambda_sd6 the introduction rate for subdomain6
##' @param lambda_sd7 the introduction rate for subdomain7
##' @param lambda_sd8 the introduction rate for subdomain8
##' @param lambda_sd9 the introduction rate for subdomain9
##' @param lambda_sd10 the introduction rate for subdomain10
##' @param lambda_sd11 the introduction rate for subdomain11
##' @param lambda_sd12 the introduction rate for subdomain12
##' @param lambda_sd13 the introduction rate for subdomain13
##' @param lambda_sd14 the introduction rate for subdomain14
##' @template beta-param
##' @return \code{SISe4}
##' @include check_arguments.R
##' @export
SISe4 <- function(init,
                  tspan,
                  events    = NULL,
                  phi       = NULL,
                  upsilon_1 = NULL,
                  upsilon_2 = NULL,
                  upsilon_3 = NULL,
                  upsilon_4 = NULL,
                  gamma_1   = NULL,
                  gamma_2   = NULL,
                  gamma_3   = NULL,
                  gamma_4   = NULL,
                  alpha     = NULL,
                  beta_t1   = NULL,
                  beta_t2   = NULL,
                  beta_t3   = NULL,
                  beta_t4   = NULL,
                  end_t1    = NULL,
                  end_t2    = NULL,
                  end_t3    = NULL,
                  end_t4    = NULL,
                  epsilon   = NULL,
                  lambda_sd1  = NULL,
                  lambda_sd2  = NULL,
                  lambda_sd3  = NULL,
                  lambda_sd4  = NULL,
                  lambda_sd5  = NULL,
                  lambda_sd6  = NULL,
                  lambda_sd7  = NULL,
                  lambda_sd8  = NULL,
                  lambda_sd9  = NULL,
                  lambda_sd10 = NULL,
                  lambda_sd11 = NULL,
                  lambda_sd12 = NULL,
                  lambda_sd13 = NULL,
                  lambda_sd14 = NULL)
{
    compartments <- c("S_1", "I_1", "S_2", "I_2", "S_3", "I_3", "S_4", "I_4")

    ## Check arguments.

    ## Check init
    if (!all(c("id", compartments) %in% names(init)))
        stop("Missing columns in init")

    ## Check initial infectious pressure
    if (is.null(phi))
        phi <- rep(0, nrow(init))
    check_infectious_pressure_arg(nrow(init), phi)

    ## Check 'gdata' parameters
    check_gdata_arg(upsilon_1, upsilon_2, upsilon_3, upsilon_4,
                    gamma_1, gamma_2, gamma_3, gamma_4,
                    alpha, beta_t1, beta_t2, beta_t3, beta_t4, epsilon,
                    lambda_sd1,  lambda_sd2,  lambda_sd3,  lambda_sd4,
                    lambda_sd5,  lambda_sd6,  lambda_sd7,  lambda_sd8,
                    lambda_sd9,  lambda_sd10, lambda_sd11, lambda_sd12,
                    lambda_sd13, lambda_sd14 )

    ## Check interval endpoints
    check_integer_arg(end_t1, end_t2, end_t3, end_t4)
    if (identical(length(end_t1), 1L))
        end_t1 <- rep(end_t1, nrow(init))
    if (identical(length(end_t2), 1L))
        end_t2 <- rep(end_t2, nrow(init))
    if (identical(length(end_t3), 1L))
        end_t3 <- rep(end_t3, nrow(init))
    if (identical(length(end_t4), 1L))
        end_t4 <- rep(end_t4, nrow(init))
    check_end_t_arg(nrow(init), end_t1, end_t2, end_t3, end_t4)

    ## Arguments seems ok...go on

    init <- init[,c("id", compartments)]

    E <- Matrix(c(1, 0, 0, 0, 1, 0, 0, 0,
                  1, 0, 0, 0, 0, 0, 0, 0,
                  0, 1, 0, 0, 0, 1, 0, 0,
                  0, 1, 0, 0, 0, 0, 0, 0,
                  0, 0, 1, 0, 0, 0, 1, 0,
                  0, 0, 1, 0, 0, 0, 0, 0,
                  0, 0, 0, 1, 0, 0, 0, 1,
                  0, 0, 0, 1, 0, 0, 0, 0),
                nrow   = 8,
                ncol   = 8,
                byrow  = TRUE,
                sparse = TRUE)
    E <- as(E, "dgCMatrix")
    colnames(E) <- as.character(seq_len(dim(E)[2]))
    rownames(E) <- compartments

    N <- Matrix(c(0, 0,
                  0, 0,
                  2, 0,
                  2, 0,
                  0, 2,
                  0, 2,
                  0, 0,
                  0, 0),
                nrow   = 8,
                ncol   = 2,
                byrow  = TRUE,
                sparse = TRUE)
    N <- as(N, "dgCMatrix")
    colnames(N) <- as.character(seq_len(dim(N)[2]))
    rownames(N) <- compartments

    ## These are the transition parameters on the rows and columns
    ## that must be updated when one of them is applied
    
    G <- Matrix(c(1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                  1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                  1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                  0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0,
                  0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0,
                  0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0,
                  0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0,
                  0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0,
                  0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0,
                  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1,
                  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1,
                  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1),
                nrow   = 12,
                ncol   = 12,
                byrow  = TRUE,
                sparse = TRUE)
    G <- as(G, "dgCMatrix")
    colnames(G) <- as.character(seq_len(dim(G)[2]))
    rownames(G) <- c("S_1 -> I_1", "I_1 -> S_1", "S_1* -> I_1*",
                     "S_2 -> I_2", "I_2 -> S_2", "S_2* -> I_2*",
                     "S_3 -> I_3", "I_3 -> S_3", "S_3* -> I_3*",
                     "S_4 -> I_4", "I_4 -> S_4", "S_4* -> I_4*")

    ## These are the reactions in the columns and the compartments on
    ## the rows
    
    S <- Matrix(c(-1,  1, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,
                   1, -1,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,
                   0,  0,  0, -1,  1, -1,  0,  0,  0,  0,  0,  0,
                   0,  0,  0,  1, -1,  1,  0,  0,  0,  0,  0,  0,
                   0,  0,  0,  0,  0,  0, -1,  1, -1,  0,  0,  0,
                   0,  0,  0,  0,  0,  0,  1, -1,  1,  0,  0,  0,
                   0,  0,  0,  0,  0,  0,  0,  0,  0, -1,  1, -1,
                   0,  0,  0,  0,  0,  0,  0,  0,  0,  1, -1,  1),
                nrow   = 8,
                ncol   = 12,
                byrow  = TRUE,
                sparse = TRUE)
    S <- as(S, "dgCMatrix")
    colnames(S) <- as.character(seq_len(dim(S)[2]))
    rownames(S) <- compartments

    v0 <- matrix(phi, nrow  = 1, byrow = TRUE)
    storage.mode(v0) <- "double"

    ldata <- matrix(c(end_t1, end_t2, end_t3, end_t4),
                    nrow  = 4,
                    byrow = TRUE)
    storage.mode(ldata) <- "double"

    gdata <- c(upsilon_1, upsilon_2, upsilon_3, upsilon_4,
               gamma_1, gamma_2, gamma_3, gamma_4,
               alpha,
               beta_t1, beta_t2, beta_t3, beta_t4,
## 14 subdomains               
               lambda_sd1,
               lambda_sd2,
               lambda_sd3,
               lambda_sd4,
               lambda_sd5,
               lambda_sd6,
               lambda_sd7,
               lambda_sd8,
               lambda_sd9,
               lambda_sd10,
               lambda_sd11,
               lambda_sd12,
               lambda_sd13,
               lambda_sd14)
    storage.mode(gdata) <- "double"

    model <- siminf_model(G      = G,
                          S      = S,
                          init   = init,
                          E      = E,
                          N      = N,
                          tspan  = tspan,
                          events = events,
                          sd     = init$sd,
                          ldata  = ldata,
                          gdata  = gdata,
                          v0     = v0)

    return(as(model, "SISe4"))
}

##' @rdname run-methods
##' @export
setMethod("run",
          signature(model = "SISe4"),
          function(model, threads, seed)
          {
              ## check that siminf_model contains all data structures
              ## required by the siminf solver and that they make sense
              validObject(model);

              .Call(SISe4_run, model, threads, seed)
          }
)

##' @rdname susceptible-methods
##' @export
setMethod("susceptible",
          signature("SISe4"),
          function(model, age = 1:4, i = NULL, by = 1, ...)
          {
              if (identical(dim(model@U), c(0L, 0L)))
                  stop("Please run the model first, the 'U' matrix is empty")

              age_categories <- 1:4
              stopifnot(all(age %in% age_categories))

              result <- NULL
              j <- seq(from = 1, to = dim(model@U)[2], by = by)

              for (k in age_categories) {
                  ## Are we interested in this age category?
                  if (k %in% age) {
                      ## Select rows for the specific age category
                      ii <- seq(from = 1 + (k - 1) * 2, to = dim(model@U)[1], by = 8)

                      ## Are we only interested in susceptible from
                      ## specific nodes?
                      if (!is.null(i))
                          ii <- ii[i]

                      ## Extract susceptible and add to result
                      if (is.null(result)) {
                          result <- as.matrix(model@U[ii, j, drop = FALSE])
                      } else {
                          result <- result + as.matrix(model@U[ii, j, drop = FALSE])
                      }
                  }
              }

              result
          }
)

##' @rdname infected-methods
##' @export
setMethod("infected",
          signature("SISe4"),
          function(model, age = 1:4, i = NULL, by = 1, ...)
          {
              if (identical(dim(model@U), c(0L, 0L)))
                  stop("Please run the model first, the 'U' matrix is empty")

              age_categories <- 1:4
              stopifnot(all(age %in% age_categories))

              result <- NULL
              j <- seq(from = 1, to = dim(model@U)[2], by = by)

              for (k in age_categories) {
                  ## Are we interested in this age category?
                  if (k %in% age) {
                      ## Select rows for the specific age category
                      ii <- seq(from = k * 2, to = dim(model@U)[1], by = 8)

                      ## Are we only interested in infected from
                      ## specific nodes?
                      if (!is.null(i))
                          ii <- ii[i]

                      ## Extract infected and add to result
                      if (is.null(result)) {
                          result <- as.matrix(model@U[ii, j, drop = FALSE])
                      } else {
                          result <- result + as.matrix(model@U[ii, j, drop = FALSE])
                      }
                  }
              }

              result
          }
)

##' @rdname prevalence-methods
##' @export
setMethod("prevalence",
          signature("SISe4"),
          function(model, age = 1:4, wnp = FALSE, i = NULL, by = 1, ...)
          {
              I <- infected(model, age = age, i = i, by = by)
              S <- susceptible(model, age = age, i = i, by = by)

              if (identical(wnp, FALSE)) {
                  I <- colSums(I)
                  S <- colSums(S)
              }

              I / (S + I)
          }
)

##' @name plot-methods
##' @aliases plot plot-methods plot,SISe4-method
##' @docType methods
##' @importFrom graphics plot
##' @export
setMethod("plot",
          signature(x = "SISe4"),
          function(x, t0 = NULL, ...)
      {
          callNextMethod(x,
                         t0 = t0,
                         legend = expression(S[1], I[1],
                                             S[2], I[2],
                                             S[3], I[3],
                                             S[4], I[4]),
                         ...)
      }
)

## This file is part of SimInf, a framework for stochastic
## disease spread simulations.
##
## Copyright (C) 2015 Pavol Bauer
## Copyright (C) 2017 -- 2019 Robin Eriksson
## Copyright (C) 2015 -- 2019 Stefan Engblom
## Copyright (C) 2015 -- 2020 Stefan Widgren
##
## SimInf is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## SimInf is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <https://www.gnu.org/licenses/>.

##' Class \code{"SimInf_model"}
##'
##' Class to handle data for the \code{SimInf_model}.
##' @template G-slot
##' @template S-slot
##' @template U-slot
##' @template U_sparse-slot
##' @slot V The result matrix for the real-valued continuous
##'     state. \code{V[, j]} contains the real-valued state of the
##'     system at \code{tspan[j]}. Numeric matrix
##'     (\eqn{N_n}\code{dim(ldata)[1]} \eqn{\times}
##'     \code{length(tspan)}).
##' @slot V_sparse If the model was configured to write the solution
##'     to a sparse matrix (\code{dgCMatrix}) the \code{V_sparse}
##'     contains the data and \code{V} is empty. The layout of the
##'     data in \code{V_sparse} is identical to \code{V}.
##' @template ldata-slot
##' @template gdata-slot
##' @template tspan-slot
##' @template u0-slot
##' @slot v0 The initial value for the real-valued continuous state.
##'     Numeric matrix (\code{dim(ldata)[1]} \eqn{\times N_n}).
##' @slot events Scheduled events \code{\linkS4class{SimInf_events}}
##' @template C_code-slot
##' @include SimInf_events.R
##' @export
##' @importFrom methods validObject
##' @importClassesFrom Matrix dgCMatrix
setClass("SimInf_model",
         slots = c(G        = "dgCMatrix",
                   S        = "dgCMatrix",
                   U        = "matrix",
                   U_sparse = "dgCMatrix",
                   ldata    = "matrix",
                   gdata    = "numeric",
                   tspan    = "numeric",
                   u0       = "matrix",
                   V        = "matrix",
                   V_sparse = "dgCMatrix",
                   v0       = "matrix",
                   events   = "SimInf_events",
                   C_code   = "character"))

##' Check if a SimInf_model object is valid
##'
##' @param object The SimInf_model_object to check.
##' @include valid.R
##' @noRd
valid_SimInf_model_object <- function(object) {
    ## Check events
    validObject(object@events)

    errors <- c(valid_tspan(object),
                valid_u0(object),
                valid_U(object),
                valid_v0(object),
                valid_V(object),
                valid_S(object),
                valid_G(object),
                valid_ldata(object),
                valid_gdata(object))

    if (length(errors))
        return(errors)
    TRUE
}

## Assign the function as the validity method for the class.
setValidity("SimInf_model", valid_SimInf_model_object)

##' Create a \code{SimInf_model}
##'
##' @template G-param
##' @template S-param
##' @template U-param
##' @template ldata-param
##' @template gdata-param
##' @template tspan-param
##' @param u0 The initial state vector. Either a matrix (\eqn{N_c
##'     \times N_n}) or a a \code{data.frame} with the number of
##'     individuals in each compartment in every node.
##' @param events A \code{data.frame} with the scheduled events.
##' @param V The result matrix for the real-valued continous
##'     compartment state (\eqn{N_n}\code{dim(ldata)[1]} \eqn{\times}
##'     \code{length(tspan)}).  \code{V[, j]} contains the real-valued
##'     state of the system at \code{tspan[j]}.
##' @param v0 The initial continuous state vector in every node.
##'     (\code{dim(ldata)[1]} \eqn{\times N_N}). The continuous state
##'     vector is updated by the specific model during the simulation
##'     in the post time step function.
##' @param E Sparse matrix to handle scheduled events, see
##'     \code{\linkS4class{SimInf_events}}.
##' @param N Sparse matrix to handle scheduled events, see
##'     \code{\linkS4class{SimInf_events}}.
##' @param C_code Character vector with optional model C code. If
##'     non-empty, the C code is written to a temporary C-file when
##'     the \code{run} method is called.  The temporary C-file is
##'     compiled and the resulting DLL is dynamically loaded. The DLL
##'     is unloaded and the temporary files are removed after running
##'     the model.
##' @return \linkS4class{SimInf_model}
##' @include init.R
##' @export
##' @importFrom methods as
##' @importFrom methods is
##' @importFrom methods new
SimInf_model <- function(G,
                         S,
                         tspan,
                         events = NULL,
                         ldata  = NULL,
                         gdata  = NULL,
                         U      = NULL,
                         u0     = NULL,
                         v0     = NULL,
                         V      = NULL,
                         E      = NULL,
                         N      = NULL,
                         C_code = NULL) {
    u0 <- init_x0(u0)
    G <- init_sparse_matrix(G)
    S <- init_sparse_matrix(S)
    ldata <- init_data_matrix(ldata)
    gdata <- init_data_vector(gdata)
    U <- init_output_matrix(U)
    v0 <- init_x0(v0, "double", TRUE)
    V <- init_output_matrix(V, "double")
    C_code <- init_C_code(C_code)
    tspan <- init_tspan(tspan)

    ## Check events
    if (!any(is.null(events), is.data.frame(events)))
        stop("'events' must be NULL or a data.frame.", call. = FALSE)
    events <- SimInf_events(E = E, N = N, events = events, t0 = tspan$t0)

    new("SimInf_model",
        G      = G,
        S      = S,
        U      = U,
        ldata  = ldata,
        gdata  = gdata,
        tspan  = tspan$tspan,
        u0     = u0,
        v0     = v0,
        V      = V,
        events = events,
        C_code = C_code)
}

##' Extract number of nodes in a model
##'
##' Extract number of nodes in a model.
##' @param model the \code{model} object to extract the number of
##'     nodes from.
##' @return the number of nodes in the model.
##' @export
##' @examples
##' ## Create an 'SIR' model with 100 nodes, with 99 susceptible,
##' ## 1 infected and 0 recovered in each node.
##' u0 <- data.frame(S = rep(99, 100), I = rep(1, 100), R = rep(0, 100))
##' model <- SIR(u0 = u0, tspan = 1:10, beta = 0.16, gamma = 0.077)
##'
##' ## Display the number of nodes in the model.
##' Nn(model)
Nn <- function(model) {
    check_model_argument(model)
    dim(model@u0)[2]
}

## Number of compartments
Nc <- function(model) {
    check_model_argument(model)
    dim(model@S)[1]
}

## Number of transitions
Nt <- function(model) {
    check_model_argument(model)
    dim(model@G)[1]
}

## Number of continuous state variables
Nd <- function(model) {
    check_model_argument(model)
    dim(model@v0)[1]
}

##' Extract global data from a \code{SimInf_model} object
##'
##' The global data is a numeric vector that is common to all nodes.
##' The global data vector is passed as an argument to the transition
##' rate functions and the post time step function.
##' @param model The \code{model} to get global data from.
##' @return a numeric vector
##' @export
##' @examples
##' ## Create an SIR model
##' model <- SIR(u0 = data.frame(S = 99, I = 1, R = 0),
##'              tspan = 1:5, beta = 0.16, gamma = 0.077)
##'
##' ## Set 'beta' to a new value
##' gdata(model, "beta") <- 2
##'
##' ## Extract the global data vector that is common to all nodes
##' gdata(model)
gdata <- function(model) {
    check_model_argument(model)
    model@gdata
}

##' Set a global data parameter for a \code{SimInf_model} object
##'
##' The global data is a numeric vector that is common to all nodes.
##' The global data vector is passed as an argument to the transition
##' rate functions and the post time step function.
##' @param model The \code{model} to set a global model parameter for.
##' @param parameter The name of the parameter to set.
##' @param value A numeric value.
##' @return a \code{SimInf_model} object
##' @export
##' @examples
##' ## Create an SIR model
##' model <- SIR(u0 = data.frame(S = 99, I = 1, R = 0),
##'              tspan = 1:5, beta = 0.16, gamma = 0.077)
##'
##' ## Set 'beta' to a new value
##' gdata(model, "beta") <- 2
##'
##' ## Extract the global data vector that is common to all nodes
##' gdata(model)
"gdata<-" <- function(model, parameter, value) {
    check_model_argument(model)

    ## Check paramter argument
    if (missing(parameter))
        stop("Missing 'parameter' argument.", call. = FALSE)
    if (!is.character(parameter))
        stop("'parameter' argument must be a character.", call. = FALSE)

    ## Check value argument
    if (missing(value))
        stop("Missing 'value' argument.", call. = FALSE)
    if (!is.numeric(value))
        stop("'value' argument must be a numeric.", call. = FALSE)

    model@gdata[parameter] <- value

    model
}

##' Extract local data from a node
##'
##' The local data is a numeric vector that is specific to a node.
##' The local data vector is passed as an argument to the transition
##' rate functions and the post time step function.
##' @param model The \code{model} to get local data from.
##' @param node index to node to extract local data from.
##' @return a numeric vector
##' @export
##' @examples
##' ## Create an 'SISe' model with 1600 nodes.
##' model <- SISe(u0 = u0_SISe(), tspan = 1:100, events = events_SISe(),
##'               phi = 0, upsilon = 1.8e-2, gamma = 0.1, alpha = 1,
##'               beta_t1 = 1.0e-1, beta_t2 = 1.0e-1, beta_t3 = 1.25e-1,
##'               beta_t4 = 1.25e-1, end_t1 = c(91, 101), end_t2 = c(182, 185),
##'               end_t3 = c(273, 275), end_t4 = c(365, 360), epsilon = 0)
##'
##' ## Display local data from the first two nodes.
##' ldata(model, node = 1)
##' ldata(model, node = 2)
ldata <- function(model, node) {
    check_model_argument(model)

    ## Check node argument
    if (missing(node))
        stop("Missing 'node' argument.", call. = FALSE)
    if (!is.numeric(node) || !identical(length(node), 1L) || node < 1)
        stop("Invalid 'node' argument.", call. = FALSE)

    model@ldata[, node]
}

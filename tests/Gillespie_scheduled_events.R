## SimInf, a framework for stochastic disease spread simulations
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

library(SimInf)

## For debugging
sessionInfo()

## TEST 1: model with inert events that are selected from the first
## comparment only: PASSES
##
##Model with 2 compartments
m <- SimInf::mparse(transitions = c("S->c1*S->I", "I->c2*I->S"),
                    compartments = c("S","I"),
                    c1 = 0.01, c2 = 0.01)
## Number of timesteps:
tsteps <- 10000
## The events
events <- data.frame(
    event      = rep(2, tsteps), ## A bunch of "internal transfers"
    time       = 1:tsteps,       ## 1 event per timestep in the simulation
    node       = 1,              ## This model only has 1 node
    dest       = 1,              ## Just 1 node
    n          = 1,              ## how many animals are moved
    proportion = 0,              ## This is not used when n > 0
    select     = 1,              ## 1st column in select matrix
    shift      = 1)              ## 1st column in the shift matrix
E <- as(matrix(as.numeric(c(1, 0)), nrow = 2), "dgCMatrix")    ## Select animals from compartment S
N <- matrix(as.numeric(c(0, 0)), nrow = 2)                     ## Do no shifting
model <- SimInf::init(m,
                      u0 = data.frame(S = 5000, I = 5000),
                      tspan = 1:tsteps,
                      events = events,
                      E = E,
                      N = N)
r1 <- run(model, threads = 1, seed = 1)
## Sans events
model <- SimInf::init(m,
                      u0 = data.frame(S = 5000, I = 5000),
                      tspan = 1:tsteps)
r2 <- run(model, threads = 1, seed = 1)
## We expect the outcome to be the same with or without events
stopifnot(identical(U(r1), U(r2)))
##
##
##TEST 2: shift an S to I and then an I to S on the same day everyday:
##PASSES
##
## Model with 2 compartments
m <- SimInf::mparse(transitions = c("S->c1*S->I", "I->c2*I->S"),
                    compartments = c("S","I"),
                    c1 = 0.01, c2 = 0.01)

## Number of timesteps:
tsteps <- 10000

## The events
events <- data.frame(
    event      = rep(2, tsteps),                ## A bunch of "internal transfers"
    time       = rep(1:tsteps, each = 2),       ## 2 events per timestep in the simulation
    node       = 1,                             ## This model only has 1 node
    dest       = 1,                             ## Just 1 node
    n          = 1,                             ## how many animals are moved
    proportion = 0,                             ## This is not used when n > 0
    select     = c(1, 2),                       ## 1st then 2nd column in select matrix
    shift      = c(1, 2))                       ## 1st then 2nd column in the shift matrix
E <- as(matrix(as.numeric(c(1, 0, 0, 1)), nrow = 2), "dgCMatrix")    ## Select animals from compartment S or I
N <- matrix(as.numeric(c(1, 0, 0, -1)), nrow = 2)                    ## Do 1 shift of S to I and I to S.
model <- SimInf::init(m,
                      u0 = data.frame(S = 5000, I = 5000),
                      tspan = 1:tsteps,
                      events = events,
                      E = E,
                      N = N)
r1 <- run(model, threads = 1, seed = 1)
## Sans events
model <- SimInf::init(m,
                      u0 = data.frame(S = 5000, I = 5000),
                      tspan = 1:tsteps)
r2 <- run(model, threads = 1, seed = 1)
## We expect the outcome to be the same with or without events
stopifnot(identical(U(r1), U(r2)))
## TEST 3: Model with inert event selected from both compartments:
## FAILS
##
## Model with 2 compartments
m <- SimInf::mparse(transitions = c("S->c1*S->I", "I->c2*I->S"),
                    compartments = c("S","I"),
                    c1 = 0.01, c2 = 0.01)
## Number of timesteps:
tsteps <- 10000
## The events
events <- data.frame(
    event      = rep(2, tsteps), ## A bunch of "internal transfers"
    time       = 1:tsteps,       ## 1 event per timestep in the simulation
    node       = 1,              ## This model only has 1 node
    dest       = 1,              ## Just 1 node
    n          = 1,              ## how many animals are moved
    proportion = 0,              ## This is not used when n > 0
    select     = 1,              ## 1st column in select matrix
    shift      = 1)              ## 1st column in the shift matrix
E <- as(matrix(as.numeric(c(1, 1)), nrow = 2), "dgCMatrix")    ## Select animals from compartment S and I
N <- matrix(as.numeric(c(0, 0)), nrow = 2)                     ## Do no shifting
model <- SimInf::init(m,
                      u0 = data.frame(S = 5000, I = 5000),
                      tspan = 1:tsteps,
                      events = events,
                      E = E,
                      N = N)
r1 <- run(model, threads = 1, seed = 1)
## Sans events
model <- SimInf::init(m,
                      u0 = data.frame(S = 5000, I = 5000),
                      tspan = 1:tsteps)
r2 <- run(model, threads = 1, seed = 1)
## We expect the outcome to be the same with or without events
stopifnot(identical(U(r1), U(r2)))
##

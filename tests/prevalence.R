## SimInf, a framework for stochastic disease spread simulations
## Copyright (C) 2015 - 2019  Stefan Widgren
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
## along with this program.  If not, see <https://www.gnu.org/licenses/>.

library("SimInf")

## For debugging
sessionInfo()

## Define a tolerance
tol = 1e-8

model <- SIR(u0 = data.frame(S = c(8, 5, 0), I = c(0, 1, 0), R = c(0, 0, 4)),
             tspan = 1:5, beta = 0.1, gamma = 0.1)

res <- tools::assertError(prevalence(model, I~.|R == 0))
stopifnot("Please run the model first, the trajectory is empty"
          == res[[1]]$message)

model@U <- matrix(c(8L, 8L, 8L, 8L, 8L,
                    0L, 0L, 0L, 0L, 0L,
                    0L, 0L, 0L, 0L, 0L,
                    5L, 4L, 3L, 2L, 1L,
                    1L, 2L, 3L, 3L, 3L,
                    0L, 0L, 0L, 1L, 2L,
                    0L, 0L, 0L, 0L, 0L,
                    0L, 0L, 0L, 0L, 0L,
                    4L, 4L, 4L, 4L, 4L),
                  ncol = 5,
                  byrow = TRUE,
                  dimnames = list(c("S", "I", "R", "S", "I", "R", "S", "I", "R"),
                                  c("1", "2", "3", "4", "5")))

res <- tools::assertError(prevalence(model))
stopifnot("Missing 'formula' argument"
          == res[[1]]$message)

res <- tools::assertError(prevalence(model, "I~."))
stopifnot("'formula' argument is not a 'formula'"
          == res[[1]]$message)

res <- tools::assertError(prevalence(model, ~I))
stopifnot("Invalid formula specification."
          == res[[1]]$message)

res <- tools::assertError(SimInf:::parse_formula_item("", "S"))
stopifnot("No compartments in formula specification."
          == res[[1]]$message)

p <- prevalence(model, I~.)$prevalence
stopifnot(all(abs(p - c(1/18, 2/18, 3/18, 3/18, 3/18)) < tol))

p <- prevalence(model, I~.|R == 0)$prevalence
stopifnot(all(abs(p - c(1/14, 2/14, 3/14, 0/8, 0/8)) < tol))

p <- prevalence(model, I~.|R > 0)$prevalence
stopifnot(all(abs(p - c(0/4, 0/4, 0/4, 3/10, 3/10)) < tol))

stopifnot(all(is.nan(prevalence(model, I~.|R == 5)$prevalence)))

res <- tools::assertError(prevalence(model, I~.|TRUE == 0))
stopifnot("The condition must be either 'TRUE' or 'FALSE' for every node and time step."
          == res[[1]]$message)

p <- prevalence(model, I~.| S == 0 | R == 0)$prevalence
stopifnot(all(abs(p - c(1/18, 2/18, 3/18, 0/12, 0/12)) < tol))

p <- prevalence(model, I~.| S == 0 | R == 0, node = 2)$prevalence
stopifnot(all(abs(p[1:3] - c(1/6, 2/6, 3/6)) < tol))
stopifnot(all(is.nan(p[4:5])))
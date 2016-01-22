/*
 *  SimInf, a framework for stochastic disease spread simulations
 *  Copyright (C) 2015  Pavol Bauer
 *  Copyright (C) 2015 - 2016  Stefan Engblom
 *  Copyright (C) 2015 - 2016  Stefan Widgren
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "SISe4.h"
#include "siminf_forward_euler_linear_decay.h"

/* Offset in integer compartment state vector */
enum {S_1, I_1, S_2, I_2, S_3, I_3, S_4, I_4};

/* Offset in real-valued continuous state vector */
enum {PHI};

/* Offsets in node local data (ldata) to parameters in the model */
enum {END_T1, END_T2, END_T3, END_T4};

/* Offsets in global data (gdata) to parameters in the model */
enum {UPSILON_1, UPSILON_2, UPSILON_3, UPSILON_4, 
      GAMMA_1, GAMMA_2, GAMMA_3, GAMMA_4,
      ALPHA, BETA_T1, BETA_T2, BETA_T3, BETA_T4, EPSILON, LAMBDA};
/* There are more than one LAMBDA parameter they are indexed by adding
   the subdomain to the LAMBDA index*/

/**
 * In age category 1; susceptible to infected: S -> I
 *
 * @param u The compartment state vector in node.
 * @param v The continuous state vector in node.
 * @param ldata The local data vector for the node.
 * @param gdata The global data vector.
 * @param t Current time.
 * @param sd The sub-domain of node.
 * @return propensity.
 */
double SISe4_S_1_to_I_1(
    const int *u,
    const double *v,
    const double *ldata,
    const double *gdata,
    double t,
    int sd)
{
    return gdata[UPSILON_1] * v[PHI] * u[S_1];
}

/**
 * In age category 1; susceptible to infected: S -> I
 * Not driven by PHI 
 * @param u The compartment state vector in node.
 * @param v The continuous state vector in node.
 * @param ldata The local data vector for the node.
 * @param gdata The global data vector.
 * @param t Current time.
 * @param sd The sub-domain of node.
 * @return propensity.
 */
double SISe4_S_1_to_I_1_intro(
    const int *u,
    const double *v,
    const double *ldata,
    const double *gdata,
    double t,
    int sd)
{
    return gdata[LAMBDA + sd] * u[S_1];
}


/**
 * In age category 2; susceptible to infected: S -> I
 *
 * @param u The compartment state vector in node.
 * @param v The continuous state vector in node.
 * @param ldata The local data vector for the node.
 * @param gdata The global data vector.
 * @param t Current time.
 * @param sd The sub-domain of node.
 * @return propensity.
 */
double SISe4_S_2_to_I_2(
    const int *u,
    const double *v,
    const double *ldata,
    const double *gdata,
    double t,
    int sd)
{
    return gdata[UPSILON_2] * v[PHI] * u[S_2];
}

/**
 * In age category 2; susceptible to infected: S -> I
 * Not driven by PHI 
 * @param u The compartment state vector in node.
 * @param v The continuous state vector in node.
 * @param ldata The local data vector for the node.
 * @param gdata The global data vector.
 * @param t Current time.
 * @param sd The sub-domain of node.
 * @return propensity.
 */
double SISe4_S_2_to_I_2_intro(
    const int *u,
    const double *v,
    const double *ldata,
    const double *gdata,
    double t,
    int sd)
{
    return gdata[LAMBDA + sd] * u[S_2];
}

/**
 *  In age category 3; susceptible to infected: S -> I
 *
 * @param u The compartment state vector in node.
 * @param v The continuous state vector in node.
 * @param ldata The local data vector for the node.
 * @param gdata The global data vector.
 * @param t Current time.
 * @param sd The sub-domain of node.
 * @return propensity.
 */
double SISe4_S_3_to_I_3(
    const int *u,
    const double *v,
    const double *ldata,
    const double *gdata,
    double t,
    int sd)
{
    return gdata[UPSILON_3] * v[PHI] * u[S_3];
}

/**
 * In age category 3; susceptible to infected: S -> I
 * Not driven by PHI 
 * @param u The compartment state vector in node.
 * @param v The continuous state vector in node.
 * @param ldata The local data vector for the node.
 * @param gdata The global data vector.
 * @param t Current time.
 * @param sd The sub-domain of node.
 * @return propensity.
 */
double SISe4_S_3_to_I_3_intro(
    const int *u,
    const double *v,
    const double *ldata,
    const double *gdata,
    double t,
    int sd)
{
    return gdata[LAMBDA + sd] * u[S_3];
}

/**
 *  In age category 4; susceptible to infected: S -> I
 *
 * @param u The compartment state vector in node.
 * @param v The continuous state vector in node.
 * @param ldata The local data vector for the node.
 * @param gdata The global data vector.
 * @param t Current time.
 * @param sd The sub-domain of node.
 * @return propensity.
 */
double SISe4_S_4_to_I_4(
    const int *u,
    const double *v,
    const double *ldata,
    const double *gdata,
    double t,
    int sd)
{
    return gdata[UPSILON_4] * v[PHI] * u[S_4];
}

/**
 * In age category 4; susceptible to infected: S -> I
 * Not driven by PHI 
 * @param u The compartment state vector in node.
 * @param v The continuous state vector in node.
 * @param ldata The local data vector for the node.
 * @param gdata The global data vector.
 * @param t Current time.
 * @param sd The sub-domain of node.
 * @return propensity.
 */
double SISe4_S_4_to_I_4_intro(
    const int *u,
    const double *v,
    const double *ldata,
    const double *gdata,
    double t,
    int sd)
{
    return gdata[LAMBDA + sd] * u[S_4];
}

/**
 *  In age category 1; infected to susceptible: I -> S
 *
 * @param u The compartment state vector in node.
 * @param v The continuous state vector in node.
 * @param ldata The local data vector for the node.
 * @param gdata The global data vector.
 * @param t Current time.
 * @param sd The sub-domain of node.
 * @return propensity.
 */
double SISe4_I_1_to_S_1(
    const int *u,
    const double *v,
    const double *ldata,
    const double *gdata,
    double t,
    int sd)
{
    return gdata[GAMMA_1] * u[I_1];
}

/**
 * In age category 2; infected to susceptible: I -> S
 *
 * @param u The compartment state vector in node.
 * @param v The continuous state vector in node.
 * @param ldata The local data vector for the node.
 * @param gdata The global data vector.
 * @param t Current time.
 * @param sd The sub-domain of node.
 * @return propensity.
 */
double SISe4_I_2_to_S_2(
    const int *u,
    const double *v,
    const double *ldata,
    const double *gdata,
    double t,
    int sd)
{
    return gdata[GAMMA_2] * u[I_2];
}

/**
 * In age category 3; infected to susceptible: I -> S
 *
 * @param u The compartment state vector in node.
 * @param v The continuous state vector in node.
 * @param ldata The local data vector for the node.
 * @param gdata The global data vector.
 * @param t Current time.
 * @param sd The sub-domain of node.
 * @return propensity
 */
double SISe4_I_3_to_S_3(
    const int *u,
    const double *v,
    const double *ldata,
    const double *gdata,
    double t,
    int sd)
{
    return gdata[GAMMA_3] * u[I_3];
}

/**
 * In age category 4; infected to susceptible: I -> S
 *
 * @param u The compartment state vector in node.
 * @param v The continuous state vector in node.
 * @param ldata The local data vector for the node.
 * @param gdata The global data vector.
 * @param t Current time.
 * @param sd The sub-domain of node.
 * @return propensity
 */
double SISe4_I_4_to_S_4(
    const int *u,
    const double *v,
    const double *ldata,
    const double *gdata,
    double t,
    int sd)
{
    return gdata[GAMMA_4] * u[I_4];
}

/**
 * Update environmental infectious pressure phi
 *
 * @param v_new The continuous state vector in the node after the post
 * time step
 * @param u The compartment state vector in the node.
 * @param v The current continuous state vector in the node.
 * @param ldata The local data vector for the node.
 * @param gdata The global data vector.
 * @param node The node.
 * @param t Current time.
 * @param sd The sub-domain of the node.
 * @return error code (<0), or 1 if node needs to update the
 * transition rates, or 0 when it doesn't need to update the
 * transition rates.
 */
int SISe4_post_time_step(
    double *v_new,
    const int *u,
    const double *v,
    const double *ldata,
    const double *gdata,
    int node,
    double t,
    int sd)
{
    const int day = (int)t % 365;
    const double I_n = u[I_1] + u[I_2] + u[I_3] + u[I_4];
    const double n = I_n + u[S_1] + u[S_2] + u[S_3] + u[S_4];
    const double phi = v[PHI];

    /* Time dependent beta in each of the four intervals of the
     * year. Forward Euler step. */
    v_new[PHI] = siminf_forward_euler_linear_decay(
        phi, day,
        ldata[END_T1], ldata[END_T2], ldata[END_T3], ldata[END_T4],
        gdata[BETA_T1], gdata[BETA_T2], gdata[BETA_T3], gdata[BETA_T4]);

    if (n > 0.0)
        v_new[PHI] += gdata[ALPHA] * I_n / n + gdata[EPSILON];
    else
        v_new[PHI] += gdata[EPSILON];

    if (isfinite(v_new[PHI]))
        return phi != v_new[PHI]; /* 1 if needs update */
    return SIMINF_ERR_V_IS_NOT_FINITE;
}

/**
 * Run simulation for the SISe4 model
 *
 * This function is called from R with '.Call'
 * @param model The SISe4 model
 * @param threads Number of threads
 * @param seed Random number seed.
 * @return S4 class SISe4 with the simulated trajectory in U
 */
SEXP SISe4_run(SEXP model, SEXP threads, SEXP seed)
{
    int err = 0;
    SEXP result, class_name;
    PropensityFun t_fun[] = {&SISe4_S_1_to_I_1, &SISe4_I_1_to_S_1, &SISe4_S_1_to_I_1_intro,
                             &SISe4_S_2_to_I_2, &SISe4_I_2_to_S_2, &SISe4_S_2_to_I_2_intro,
                             &SISe4_S_3_to_I_3, &SISe4_I_3_to_S_3, &SISe4_S_3_to_I_3_intro,
			     &SISe4_S_4_to_I_4, &SISe4_I_4_to_S_4, &SISe4_S_4_to_I_4_intro};

    if (R_NilValue == model || S4SXP != TYPEOF(model))
        Rf_error("Invalid SISe4 model");

    class_name = getAttrib(model, R_ClassSymbol);
    if (strcmp(CHAR(STRING_ELT(class_name, 0)), "SISe4") != 0)
        Rf_error("Invalid SISe4 model: %s", CHAR(STRING_ELT(class_name, 0)));

    result = PROTECT(duplicate(model));

    err = siminf_run(result, threads, seed, t_fun, &SISe4_post_time_step);

    UNPROTECT(1);

    if (err)
        siminf_error(err);

    return result;
}

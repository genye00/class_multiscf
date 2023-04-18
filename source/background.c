/** @file background.c Documented background module
 *
 * * Julien Lesgourgues, 17.04.2011
 * * routines related to ncdm written by T. Tram in 2011
 * * new integration scheme written by N. Schoeneberg in 2020
 *
 * Deals with the cosmological background evolution.
 * This module has two purposes:
 *
 * - at the beginning, to initialize the background, i.e. to integrate
 *    the background equations, and store all background quantities
 *    as a function of conformal time inside an interpolation table.
 *
 * - to provide routines which allow other modules to evaluate any
 *    background quantity for a given value of the conformal time (by
 *    interpolating within the interpolation table), or to find the
 *    correspondence between redshift and conformal time.
 *
 *
 * The overall logic in this module is the following:
 *
 * 1. most background parameters that we will call {A}
 * (e.g. rho_gamma, ..) can be expressed as simple analytical
 * functions of the scale factor 'a' plus a few variables that we will
 * call {B} (e.g. (phi, phidot) for quintessence, or some temperature
 * for exotic particles, etc...). [Side note: for simplicity, all variables
 * {B} are declared redundently inside {A}.]
 *
 * 2. in turn, quantities {B} can be found as a function of the the
 * scale factor [or rather (a/a_0)] by integrating the
 * background equations. Thus {B} also includes the density of species
 * which energy conservation equation must be integrated explicitely,
 * like the density of fluids or of decaying dark matter.
 *
 * 3. some other quantities that we will call {C} (like e.g. proper
 * and conformal time, the sound horizon, the analytic scale-invariant
 * growth factor) also require an explicit integration with respect to
 * (a/a_0) [or rather log(a/a_p)], since they cannot be inferred
 * analytically from (a/a_0) and parameters {B}. The difference
 * between {B} and {C} parameters is that {C} parameters do not need
 * to be known in order to get {A}.
 *
 * So, we define the following routines:
 *
 * - background_functions() returns all background quantities {A} as a
 *    function of (a/a_0) and of quantities {B}.
 *
 * - background_solve() integrates the quantities {B} and {C} with
 *    respect to log(a/a_0); this integration requires many calls to
 *    background_functions().
 *
 * - the result is stored in the form of a big table in the background
 *    structure. There is one column for the scale factor, and one for
 *    each quantity {A} or {C} [Side note: we don;t include {B} here
 *    because the {B} variables are already decalred redundently also
 *    as {A} quantitites.]
 *
 * Later in the code:
 *
 * - If we know the variables (a/a_0) + {B} and need some quantity {A}
 *    (but not {C}), the quickest and most precise way is to call
 *    directly background_functions() (for instance, in simple models,
 *    if we want H at a given value of the scale factor).
 *
 * - If we know 'tau' and want any other quantity, we can call
 *    background_at_tau(), which interpolates in the table and returns
 *    all values.
 *
 * - If we know 'z' but not the {B} variables, or if we know 'z' and
 *    we want {C} variables, we need to call background_at_z(), which
 *    interpolates in the table and returns all values.
 *
 * - Finally, it can be useful to get 'tau' for a given redshift 'z'
 *    or vice-versa: this can be done with background_tau_of_z() or
 *    background_z_of_tau().
 *
 *
 * In order to save time, background_at_tau() ans background_at_z()
 * can be called in three modes: short_info, normal_info, long_info
 * (returning only essential quantities, or useful quantities, or
 * rarely useful quantities). Each line in the interpolation table is
 * a vector whose first few elements correspond to the short_info
 * format; a larger fraction contribute to the normal format; and the
 * full vector corresponds to the long format. The guideline is that
 * short_info returns only geometric quantities like a, H, H'; normal
 * format returns quantities strictly needed at each step in the
 * integration of perturbations; long_info returns quantities needed
 * only occasionally.
 *
 * In summary, the following functions can be called from other modules:
 *
 * -# background_init() at the beginning background_at_tau(),
 * -# background_at_z(), background_tau_of_z(), background_z_of_tau() at any later time
 * -# background_free() at the end, when no more calls to the previous functions are needed
 *
 * For units and normalisation conventions, there are two guiding principles:
 *
 * 1) All quantities are expressed in natural units in which everything is in powers of Mpc, e.g.:
 *
 * - t stands for (cosmological or proper time)*c in Mpc
 * - tau stands for (conformal time)*c in Mpc
 * - H stands for (Hubble parameter)/c in \f$ Mpc^{-1} \f$
 * - etc.
 *
 * 2) New since v3.0: all quantities that should normally scale with some power of
 * a_0^n are renormalised by a_0^{-n}, in order to be independent of a_0, e.g.
 *
 * - a in the code stands for \f$ a/a_0 \f$ in reality
 * - tau in the code stands for \f$ a_0 \tau c \f$ in Mpc
 * - any prime in the code stands for \f$ (1/a_0) d/d\tau \f$
 * - r stands for any comoving radius times a_0
 * - etc.
 */

#include "background.h"

/**
 * Background quantities at given redshift z.
 *
 * Evaluates all background quantities at a given value of
 * redshift by reading the pre-computed table and interpolating.
 *
 * @param pba           Input: pointer to background structure (containing pre-computed table)
 * @param z             Input: redshift
 * @param return_format Input: format of output vector (short_info, normal_info, long_info)
 * @param inter_mode     Input: interpolation mode (normal or closeby)
 * @param last_index    Input/Output: index of the previous/current point in the interpolation array (input only for closeby mode, output for both)
 * @param pvecback      Output: vector (assumed to be already allocated)
 * @return the error status
 */

int background_at_z(
    struct background *pba,
    double z,
    enum vecback_format return_format,
    enum interpolation_method inter_mode,
    int * last_index,
    double * pvecback /* vector with argument pvecback[index_bg] (must be already allocated with a size compatible with return_format) */
)
{

    /** Summary: */

    /** - define local variables */

    /* size of output vector, controlled by input parameter return_format */
    int pvecback_size;

    /* log(a) (in fact, given our normalisation conventions, this is log(a/a_0)) */
    double loga;

    /** - check that log(a) = log(1/(1+z)) = -log(1+z) is in the pre-computed range */
    loga = -log(1 + z);

    class_test(loga < pba->loga_table[0],
               pba->error_message,
               "out of range: a/a_0 = %e < a_min/a_0 = %e, you should decrease the precision parameter a_ini_over_a_today_default\n",
               1. / (1. + z), exp(pba->loga_table[0]));

    class_test(loga > pba->loga_table[pba->bt_size - 1],
               pba->error_message,
               "out of range: a/a_0 = %e > a_max/a_0 = %e\n", 1. / (1. + z),
               exp(pba->loga_table[pba->bt_size - 1]));

    /** - deduce length of returned vector from format mode */

    if(return_format == normal_info)
    {
        pvecback_size = pba->bg_size_normal;
    }
    else
    {
        if(return_format == short_info)
        {
            pvecback_size = pba->bg_size_short;
        }
        else
        {
            pvecback_size = pba->bg_size;
        }
    }


    /** - interpolate from pre-computed table with array_interpolate()
        or array_interpolate_growing_closeby() (depending on
        interpolation mode) */

    if(inter_mode == inter_normal)
    {
        class_call(array_interpolate_spline(
                       pba->loga_table,
                       pba->bt_size,
                       pba->background_table,
                       pba->d2background_dloga2_table,
                       pba->bg_size,
                       loga,
                       last_index,
                       pvecback,
                       pvecback_size,
                       pba->error_message),
                   pba->error_message,
                   pba->error_message);
    }
    if(inter_mode == inter_closeby)
    {
        class_call(array_interpolate_spline_growing_closeby(
                       pba->loga_table,
                       pba->bt_size,
                       pba->background_table,
                       pba->d2background_dloga2_table,
                       pba->bg_size,
                       loga,
                       last_index,
                       pvecback,
                       pvecback_size,
                       pba->error_message),
                   pba->error_message,
                   pba->error_message);
    }

    return _SUCCESS_;
}

/**
 * Background quantities at given conformal time tau.
 *
 * Evaluates all background quantities at a given value of
 * conformal time by reading the pre-computed table and interpolating.
 *
 * @param pba           Input: pointer to background structure (containing pre-computed table)
 * @param tau           Input: value of conformal time
 * @param return_format Input: format of output vector (short_info, normal_info, long_info)
 * @param inter_mode     Input: interpolation mode (normal or closeby)
 * @param last_index    Input/Output: index of the previous/current point in the interpolation array (input only for closeby mode, output for both)
 * @param pvecback      Output: vector (assumed to be already allocated)
 * @return the error status
 */

int background_at_tau(
    struct background *pba,
    double tau,
    enum vecback_format return_format,
    enum interpolation_method inter_mode,
    int * last_index,
    double * pvecback /* vector with argument pvecback[index_bg] (must be already allocated with a size compatible with return_format) */
)
{

    /** Summary: */

    /** - define local variables */
    double z;

    /** - Get current redshift */
    class_call(background_z_of_tau(pba, tau, &z),
               pba->error_message,
               pba->error_message);

    /** - Get background at corresponding redshift */
    class_call(background_at_z(pba, z, return_format, inter_mode, last_index, pvecback),
               pba->error_message,
               pba->error_message);

    return _SUCCESS_;
}

/**
 * Conformal time at given redshift.
 *
 * Returns tau(z) by interpolation from pre-computed table.
 *
 * @param pba Input: pointer to background structure
 * @param z   Input: redshift
 * @param tau Output: conformal time
 * @return the error status
 */

int background_tau_of_z(
    struct background *pba,
    double z,
    double * tau
)
{

    /** Summary: */

    /** - define local variables */

    /* necessary for calling array_interpolate(), but never used */
    int last_index;

    /** - check that \f$ z \f$ is in the pre-computed range */
    class_test(z < pba->z_table[pba->bt_size - 1],
               pba->error_message,
               "out of range: z=%e < z_min=%e\n", z, pba->z_table[pba->bt_size - 1]);

    class_test(z > pba->z_table[0],
               pba->error_message,
               "out of range: z=%e > z_max=%e\n", z, pba->z_table[0]);

    /** - interpolate from pre-computed table with array_interpolate() */
    class_call(array_interpolate_spline(
                   pba->z_table,
                   pba->bt_size,
                   pba->tau_table,
                   pba->d2tau_dz2_table,
                   1,
                   z,
                   &last_index,
                   tau,
                   1,
                   pba->error_message),
               pba->error_message,
               pba->error_message);

    return _SUCCESS_;
}
/**
 * Redshift at given conformal time.
 *
 * Returns z(tau) by interpolation from pre-computed table.
 *
 * @param pba Input: pointer to background structure
 * @param tau Input: conformal time
 * @param z   Output: redshift
 * @return the error status
 */

int background_z_of_tau(
    struct background *pba,
    double tau,
    double * z
)
{

    /** Summary: */

    /** - define local variables */

    /* necessary for calling array_interpolate(), but never used */
    int last_index;

    /** - check that \f$ tau \f$ is in the pre-computed range */
    class_test(tau < pba->tau_table[0],
               pba->error_message,
               "out of range: tau=%e < tau_min=%e\n", tau, pba->tau_table[0]);

    class_test(tau > pba->tau_table[pba->bt_size - 1],
               pba->error_message,
               "out of range: tau=%e > tau_max=%e\n", tau, pba->tau_table[pba->bt_size - 1]);

    /** - interpolate from pre-computed table with array_interpolate() */
    class_call(array_interpolate_spline(
                   pba->tau_table,
                   pba->bt_size,
                   pba->z_table,
                   pba->d2z_dtau2_table,
                   1,
                   tau,
                   &last_index,
                   z,
                   1,
                   pba->error_message),
               pba->error_message,
               pba->error_message);

    return _SUCCESS_;
}

/**
 * Function evaluating all background quantities which can be computed
 * analytically as a function of a and of {B} quantities (see
 * discussion at the beginning of this file).
 *
 * @param pba           Input: pointer to background structure
 * @param a             Input: scale factor (in fact, with our normalisation conventions, this is (a/a_0) )
 * @param pvecback_B    Input: vector containing all {B} quantities
 * @param return_format Input: format of output vector
 * @param pvecback      Output: vector of background quantities (assumed to be already allocated)
 * @return the error status
 */

int background_functions(
    struct background * pba,
    double a,
    double * pvecback_B, /* vector with argument pvecback[index_bi] */
    enum vecback_format return_format,
    double * pvecback /* vector with argument pvecback[index_bg] (must be already allocated with a size compatible with return_format) */
)
{

    /** Summary: */

    /** - define local variables */

    /* total density */
    double rho_tot;
    /* critical density */
    double rho_crit;
    /* total pressure */
    double p_tot;
    /* total relativistic density */
    double rho_r;
    /* total non-relativistic density */
    double rho_m;
    /* background ncdm quantities */
    double rho_ncdm, p_ncdm, pseudo_p_ncdm;
    /* index for n_ncdm species */
    int n_ncdm;
    /* index for ndescf components */
    int n_ndescf;
    /* fluid's time-dependent equation of state parameter */
    double w_fld, dw_over_da, integral_fld;
    /* scalar field quantities */
    double phi, phi_prime;
    /* Since we only know a_prime_over_a after we have rho_tot,
       it is not possible to simply sum up p_tot_prime directly.
       Instead we sum up dp_dloga = p_prime/a_prime_over_a. The formula is
       p_prime = a_prime_over_a * dp_dloga = a_prime_over_a * Sum [ (w_prime/a_prime_over_a -3(1+w)w)rho].
       Note: The scalar field contribution must be added in the end, as an exception!*/
    double dp_dloga;

    /** - initialize local variables */
    rho_tot = 0.;
    p_tot = 0.;
    dp_dloga = 0.;
    rho_r = 0.;
    rho_m = 0.;

    class_test(a <= 0.,
               pba->error_message,
               "a = %e instead of strictly positive", a);

    /** - pass value of \f$ a\f$ to output */
    pvecback[pba->index_bg_a] = a;

    /** - compute each component's density and pressure */

    /* photons */
    pvecback[pba->index_bg_rho_g] = pba->Omega0_g * pow(pba->H0, 2) / pow(a, 4);
    rho_tot += pvecback[pba->index_bg_rho_g];
    p_tot += (1. / 3.) * pvecback[pba->index_bg_rho_g];
    dp_dloga += -(4. / 3.) * pvecback[pba->index_bg_rho_g];
    rho_r += pvecback[pba->index_bg_rho_g];

    /* baryons */
    pvecback[pba->index_bg_rho_b] = pba->Omega0_b * pow(pba->H0, 2) / pow(a, 3);
    rho_tot += pvecback[pba->index_bg_rho_b];
    p_tot += 0;
    rho_m += pvecback[pba->index_bg_rho_b];

    /* cdm */
    if(pba->has_cdm == _TRUE_)
    {
        pvecback[pba->index_bg_rho_cdm] = pba->Omega0_cdm * pow(pba->H0, 2) / pow(a, 3);
        rho_tot += pvecback[pba->index_bg_rho_cdm];
        p_tot += 0.;
        rho_m += pvecback[pba->index_bg_rho_cdm];
    }

    /* dcdm */
    if(pba->has_dcdm == _TRUE_)
    {
        /* Pass value of rho_dcdm to output */
        pvecback[pba->index_bg_rho_dcdm] = pvecback_B[pba->index_bi_rho_dcdm];
        rho_tot += pvecback[pba->index_bg_rho_dcdm];
        p_tot += 0.;
        rho_m += pvecback[pba->index_bg_rho_dcdm];
    }

    /* dr */
    if(pba->has_dr == _TRUE_)
    {
        /* Pass value of rho_dr to output */
        pvecback[pba->index_bg_rho_dr] = pvecback_B[pba->index_bi_rho_dr];
        rho_tot += pvecback[pba->index_bg_rho_dr];
        p_tot += (1. / 3.) * pvecback[pba->index_bg_rho_dr];
        dp_dloga += -(4. / 3.) * pvecback[pba->index_bg_rho_dr];
        rho_r += pvecback[pba->index_bg_rho_dr];
    }

    /* Non dark energy scalar fields, by Gen Ye */
    if(pba->has_ndescf == _TRUE_)
    {
        for(n_ndescf = 0; n_ndescf < pba->N_ndescf; n_ndescf++)
        {
            /* This block of information is irrelevant in fld approx, keep track of it anyway */
            phi = pvecback_B[pba->index_bi_phi_ndescf1 + n_ndescf];
            phi_prime = pvecback_B[pba->index_bi_phi_prime_ndescf1 + n_ndescf];
            pvecback[pba->index_bg_phi_ndescf1 + n_ndescf] = phi;
            pvecback[pba->index_bg_phi_prime_ndescf1 + n_ndescf] = phi_prime;
            pvecback[pba->index_bg_V_ndescf1 + n_ndescf] = V_ndescf(pba, phi, n_ndescf);
            pvecback[pba->index_bg_dV_ndescf1 + n_ndescf] = dV_ndescf(pba, phi, n_ndescf);
            pvecback[pba->index_bg_ddV_ndescf1 + n_ndescf] = ddV_ndescf(pba, phi, n_ndescf);
            /* this block */

            if(pba->ndescf_use_kg[n_ndescf] == _TRUE_)
            {

                pvecback[pba->index_bg_rho_ndescf1 + n_ndescf] = (phi_prime * phi_prime / (2 * a * a) + V_ndescf(
                            pba, phi, n_ndescf)) / 3.;
                pvecback[pba->index_bg_p_ndescf1 + n_ndescf] = (phi_prime * phi_prime / (2 * a * a) - V_ndescf(pba,
                        phi, n_ndescf)) / 3.;
                /* rho_ndescf becomes the dynamic variable in fld approx, in the KG regime we update it here. Is it safe? */
                pvecback_B[pba->index_bi_rho_ndescf1 + n_ndescf] = pvecback[pba->index_bg_rho_ndescf1 + n_ndescf];
            }
            else
            {
                pvecback[pba->index_bg_rho_ndescf1 + n_ndescf] = pvecback_B[pba->index_bi_rho_ndescf1 + n_ndescf];
                pvecback[pba->index_bg_p_ndescf1 + n_ndescf] = pba->ndescf_w_fld[n_ndescf] *
                        pvecback_B[pba->index_bi_rho_ndescf1 + n_ndescf];
            }

            rho_tot += pvecback[pba->index_bg_rho_ndescf1 + n_ndescf];
            p_tot += pvecback[pba->index_bg_p_ndescf1 + n_ndescf];
            dp_dloga += 0.0; /** <-- This depends on a_prime_over_a, so we cannot add it now! */
            rho_r += 3.*pvecback[pba->index_bg_p_ndescf1 + n_ndescf];
            rho_m += pvecback[pba->index_bg_rho_ndescf1 + n_ndescf] - 3.*pvecback[pba->index_bg_p_ndescf1 +
                     n_ndescf];
        }
    }

    /* Scalar field */
    if(pba->has_scf == _TRUE_)
    {
        phi = pvecback_B[pba->index_bi_phi_scf];
        phi_prime = pvecback_B[pba->index_bi_phi_prime_scf];
        pvecback[pba->index_bg_phi_scf] = phi; // value of the scalar field phi
        pvecback[pba->index_bg_phi_prime_scf] =
            phi_prime; // value of the scalar field phi derivative wrt conformal time
        pvecback[pba->index_bg_V_scf] = V_scf(pba,
                                              phi); //V_scf(pba,phi); //write here potential as function of phi
        pvecback[pba->index_bg_dV_scf] = dV_scf(pba,
                                                phi); // dV_scf(pba,phi); //potential' as function of phi
        pvecback[pba->index_bg_ddV_scf] = ddV_scf(pba,
                                          phi); // ddV_scf(pba,phi); //potential'' as function of phi
        pvecback[pba->index_bg_rho_scf] = (phi_prime * phi_prime / (2 * a * a) + V_scf(pba,
                                           phi)) / 3.; // energy of the scalar field. The field units are set automatically by setting the initial conditions
        pvecback[pba->index_bg_p_scf] = (phi_prime * phi_prime / (2 * a * a) - V_scf(pba,
                                         phi)) / 3.; // pressure of the scalar field
        rho_tot += pvecback[pba->index_bg_rho_scf];
        p_tot += pvecback[pba->index_bg_p_scf];
        dp_dloga += 0.0; /** <-- This depends on a_prime_over_a, so we cannot add it now! */
        //divide relativistic & nonrelativistic (not very meaningful for oscillatory models)
        rho_r += 3.*pvecback[pba->index_bg_p_scf]; //field pressure contributes radiation
        rho_m += pvecback[pba->index_bg_rho_scf] - 3.*
                 pvecback[pba->index_bg_p_scf]; //the rest contributes matter
        //printf(" a= %e, Omega_scf = %f, \n ",a, pvecback[pba->index_bg_rho_scf]/rho_tot );
    }

    /* ncdm */
    if(pba->has_ncdm == _TRUE_)
    {

        /* Loop over species: */
        for(n_ncdm = 0; n_ncdm < pba->N_ncdm; n_ncdm++)
        {

            /* function returning background ncdm[n_ncdm] quantities (only
               those for which non-NULL pointers are passed) */
            class_call(background_ncdm_momenta(
                           pba->q_ncdm_bg[n_ncdm],
                           pba->w_ncdm_bg[n_ncdm],
                           pba->q_size_ncdm_bg[n_ncdm],
                           pba->M_ncdm[n_ncdm],
                           pba->factor_ncdm[n_ncdm],
                           1. / a - 1.,
                           NULL,
                           &rho_ncdm,
                           &p_ncdm,
                           NULL,
                           &pseudo_p_ncdm),
                       pba->error_message,
                       pba->error_message);

            pvecback[pba->index_bg_rho_ncdm1 + n_ncdm] = rho_ncdm;
            rho_tot += rho_ncdm;
            pvecback[pba->index_bg_p_ncdm1 + n_ncdm] = p_ncdm;
            p_tot += p_ncdm;
            pvecback[pba->index_bg_pseudo_p_ncdm1 + n_ncdm] = pseudo_p_ncdm;
            /** See e.g. Eq. A6 in 1811.00904. */
            dp_dloga += (pseudo_p_ncdm - 5 * p_ncdm);

            /* (3 p_ncdm1) is the "relativistic" contribution to rho_ncdm1 */
            rho_r += 3.* p_ncdm;

            /* (rho_ncdm1 - 3 p_ncdm1) is the "non-relativistic" contribution
               to rho_ncdm1 */
            rho_m += rho_ncdm - 3.* p_ncdm;
        }
    }

    /* Lambda */
    if(pba->has_lambda == _TRUE_)
    {
        pvecback[pba->index_bg_rho_lambda] = pba->Omega0_lambda * pow(pba->H0, 2);
        rho_tot += pvecback[pba->index_bg_rho_lambda];
        p_tot -= pvecback[pba->index_bg_rho_lambda];
    }

    /* fluid with w(a) and constant cs2 */
    if(pba->has_fld == _TRUE_)
    {

        /* get rho_fld from vector of integrated variables */
        pvecback[pba->index_bg_rho_fld] = pvecback_B[pba->index_bi_rho_fld];

        /* get w_fld from dedicated function */
        class_call(background_w_fld(pba, a, &w_fld, &dw_over_da, &integral_fld), pba->error_message,
                   pba->error_message);
        pvecback[pba->index_bg_w_fld] = w_fld;

        // Obsolete: at the beginning, we had here the analytic integral solution corresponding to the case w=w0+w1(1-a/a0):
        // pvecback[pba->index_bg_rho_fld] = pba->Omega0_fld * pow(pba->H0,2) / pow(a,3.*(1.+pba->w0_fld+pba->wa_fld)) * exp(3.*pba->wa_fld*(a-1.));
        // But now everthing is integrated numerically for a given w_fld(a) defined in the function background_w_fld.

        rho_tot += pvecback[pba->index_bg_rho_fld];
        p_tot += w_fld * pvecback[pba->index_bg_rho_fld];
        dp_dloga += (a * dw_over_da - 3 * (1 + w_fld) * w_fld) * pvecback[pba->index_bg_rho_fld];
    }

    /* relativistic neutrinos (and all relativistic relics) */
    if(pba->has_ur == _TRUE_)
    {
        pvecback[pba->index_bg_rho_ur] = pba->Omega0_ur * pow(pba->H0, 2) / pow(a, 4);
        rho_tot += pvecback[pba->index_bg_rho_ur];
        p_tot += (1. / 3.) * pvecback[pba->index_bg_rho_ur];
        dp_dloga += -(4. / 3.) * pvecback[pba->index_bg_rho_ur];
        rho_r += pvecback[pba->index_bg_rho_ur];
    }

    /* interacting dark matter */
    if(pba->has_idm_dr == _TRUE_)
    {
        pvecback[pba->index_bg_rho_idm_dr] = pba->Omega0_idm_dr * pow(pba->H0, 2) / pow(a, 3);
        rho_tot += pvecback[pba->index_bg_rho_idm_dr];
        p_tot += 0.;
        rho_m += pvecback[pba->index_bg_rho_idm_dr];
    }

    /* interacting dark radiation */
    if(pba->has_idr == _TRUE_)
    {
        pvecback[pba->index_bg_rho_idr] = pba->Omega0_idr * pow(pba->H0, 2) / pow(a, 4);
        rho_tot += pvecback[pba->index_bg_rho_idr];
        p_tot += (1. / 3.) * pvecback[pba->index_bg_rho_idr];
        rho_r += pvecback[pba->index_bg_rho_idr];
    }

    /** - compute expansion rate H from Friedmann equation: this is the
        only place where the Friedmann equation is assumed. Remember
        that densities are all expressed in units of \f$ [3c^2/8\pi G] \f$, ie
        \f$ \rho_{class} = [8 \pi G \rho_{physical} / 3 c^2]\f$ */
    pvecback[pba->index_bg_H] = sqrt(rho_tot - pba->K / a / a);

    /** - compute derivative of H with respect to conformal time */
    pvecback[pba->index_bg_H_prime] = - (3. / 2.) * (rho_tot + p_tot) * a + pba->K / a;

    /* Total energy density*/
    pvecback[pba->index_bg_rho_tot] = rho_tot;

    /* Total pressure */
    pvecback[pba->index_bg_p_tot] = p_tot;

    /* Derivative of total pressure w.r.t. conformal time */
    pvecback[pba->index_bg_p_tot_prime] = a * pvecback[pba->index_bg_H] * dp_dloga;
    if(pba->has_scf == _TRUE_)
    {
        /** The contribution of scf was not added to dp_dloga, add p_scf_prime here: */
        pvecback[pba->index_bg_p_prime_scf] = pvecback[pba->index_bg_phi_prime_scf] *
                                              (-pvecback[pba->index_bg_phi_prime_scf] * pvecback[pba->index_bg_H] / a - 2. /
                                               3.*pvecback[pba->index_bg_dV_scf]);
        pvecback[pba->index_bg_p_tot_prime] += pvecback[pba->index_bg_p_prime_scf];
    }
    /* Loop over contribution to dp_dloga for each ndescf, by Gen Ye */
    if(pba->has_ndescf == _TRUE_)
    {
        for(n_ndescf = 0; n_ndescf < pba->N_ndescf; n_ndescf++)
        {
            if(pba->ndescf_use_kg[n_ndescf] == _TRUE_)
            {
                pvecback[pba->index_bg_p_prime_ndescf1 + n_ndescf] = pvecback[pba->index_bg_phi_prime_ndescf1 +
                        n_ndescf] * (-pvecback[pba->index_bg_phi_prime_ndescf1 + n_ndescf] * pvecback[pba->index_bg_H] / a -
                                     2. / 3.*pvecback[pba->index_bg_dV_ndescf1 + n_ndescf]);
            }
            else
            {
                /* p_prime = w * rho_prime = -3*a*H*w*(1+w)*rho in fld approx */
                w_fld = pba->ndescf_w_fld[n_ndescf];
                pvecback[pba->index_bg_p_prime_ndescf1 + n_ndescf] = -3.* a * pvecback[pba->index_bg_H] * w_fld *
                        (1. + w_fld) * pvecback[pba->index_bg_rho_ndescf1 + n_ndescf];
            }
            pvecback[pba->index_bg_p_tot_prime] += pvecback[pba->index_bg_p_prime_ndescf1 + n_ndescf];
        }
    }
    /* end of modification by Gen Ye */

    /** - compute critical density */
    rho_crit = rho_tot - pba->K / a / a;
    class_test(rho_crit <= 0.,
               pba->error_message,
               "rho_crit = %e instead of strictly positive", rho_crit);

    /** - compute relativistic density to total density ratio */
    pvecback[pba->index_bg_Omega_r] = rho_r / rho_crit;

    /** - compute other quantities in the exhaustive, redundant format */
    if(return_format == long_info)
    {

        /** - store critical density */
        pvecback[pba->index_bg_rho_crit] = rho_crit;

        /** - compute Omega_m */
        pvecback[pba->index_bg_Omega_m] = rho_m / rho_crit;

        /** - cosmological time */
        pvecback[pba->index_bg_time] = pvecback_B[pba->index_bi_time];

        /** - comoving sound horizon */
        pvecback[pba->index_bg_rs] = pvecback_B[pba->index_bi_rs];

        /** - growth factor */
        pvecback[pba->index_bg_D] = pvecback_B[pba->index_bi_D];

        /** - velocity growth factor */
        pvecback[pba->index_bg_f] = pvecback_B[pba->index_bi_D_prime] / (pvecback_B[pba->index_bi_D] * a *
                                    pvecback[pba->index_bg_H]);

        /* one can put other variables here */
        /*  */
        /*  */

    }

    return _SUCCESS_;

}

/**
 * Single place where the fluid equation of state is
 * defined. Parameters of the function are passed through the
 * background structure. Generalisation to arbitrary functions should
 * be simple.
 *
 * @param pba            Input: pointer to background structure
 * @param a              Input: current value of scale factor (in fact, with our conventions, of (a/a_0))
 * @param w_fld          Output: equation of state parameter w_fld(a)
 * @param dw_over_da_fld Output: function dw_fld/da
 * @param integral_fld   Output: function \f$ \int_{a}^{a_0} da 3(1+w_{fld})/a \f$
 * @return the error status
 */

int background_w_fld(
    struct background * pba,
    double a,
    double * w_fld,
    double * dw_over_da_fld,
    double * integral_fld
)
{

    double Omega_ede = 0.;
    double dOmega_ede_over_da = 0.;
    double d2Omega_ede_over_da2 = 0.;
    double a_eq, Omega_r, Omega_m;

    /** - first, define the function w(a) */
    switch(pba->fluid_equation_of_state)
    {
    case CLP:
        *w_fld = pba->w0_fld + pba->wa_fld * (1. - a);
        break;
    case EDE:
        // Omega_ede(a) taken from eq. (10) in 1706.00730
        Omega_ede = (pba->Omega0_fld - pba->Omega_EDE * (1. - pow(a, -3.*pba->w0_fld)))
                    / (pba->Omega0_fld + (1. - pba->Omega0_fld) * pow(a, 3.*pba->w0_fld))
                    + pba->Omega_EDE * (1. - pow(a, -3.*pba->w0_fld));

        // d Omega_ede / d a taken analytically from the above
        dOmega_ede_over_da = - pba->Omega_EDE * 3.*pba->w0_fld * pow(a,
                             -3.*pba->w0_fld - 1.) / (pba->Omega0_fld + (1. - pba->Omega0_fld) * pow(a, 3.*pba->w0_fld))
                             - (pba->Omega0_fld - pba->Omega_EDE * (1. - pow(a,
                                     -3.*pba->w0_fld))) * (1. - pba->Omega0_fld) * 3.*pba->w0_fld * pow(a,
                                             3.*pba->w0_fld - 1.) / pow(pba->Omega0_fld + (1. - pba->Omega0_fld) * pow(a, 3.*pba->w0_fld), 2)
                             + pba->Omega_EDE * 3.*pba->w0_fld * pow(a, -3.*pba->w0_fld - 1.);

        // find a_equality (needed because EDE tracks first radiation, then matter)
        Omega_r = pba->Omega0_g * (1. + 3.044 * 7. / 8.*pow(4. / 11.,
                                   4. / 3.)); // assumes LambdaCDM + eventually massive neutrinos so light that they are relativistic at equality; needs to be generalised later on.
        Omega_m = pba->Omega0_b;
        if(pba->has_cdm == _TRUE_) Omega_m += pba->Omega0_cdm;
        if(pba->has_idm_dr == _TRUE_) Omega_m += pba->Omega0_idm_dr;
        if(pba->has_dcdm == _TRUE_)
            class_stop(pba->error_message,
                       "Early Dark Energy not compatible with decaying Dark Matter because we omitted to code the calculation of a_eq in that case, but it would not be difficult to add it if necessary, should be a matter of 5 minutes");
        a_eq = Omega_r / Omega_m; // assumes a flat universe with a=1 today

        // w_ede(a) taken from eq. (11) in 1706.00730
        *w_fld = - dOmega_ede_over_da * a / Omega_ede / 3. / (1. - Omega_ede) + a_eq / 3. / (a + a_eq);
        break;
    /* Acoustic Dark Energy (ADE) from arXiv:1905.12618 */
    case ADE:
        *w_fld = (1 + pba->wf_ADE) / pow(1 + pow(pba->a_ADE / a, 6 * (1 + pba->wf_ADE)), 0.5) - 1.;
        break;
    }


    /** - then, give the corresponding analytic derivative dw/da (used
        by perturbation equations; we could compute it numerically,
        but with a loss of precision; as long as there is a simple
        analytic expression of the derivative of the previous
        function, let's use it! */
    switch(pba->fluid_equation_of_state)
    {
    case CLP:
        *dw_over_da_fld = - pba->wa_fld;
        break;
    case EDE:
        d2Omega_ede_over_da2 = 0.;
        *dw_over_da_fld = - d2Omega_ede_over_da2 * a / 3. / (1. - Omega_ede) / Omega_ede
                          - dOmega_ede_over_da / 3. / (1. - Omega_ede) / Omega_ede
                          + dOmega_ede_over_da * dOmega_ede_over_da * a / 3. / (1. - Omega_ede) / (1. - Omega_ede) / Omega_ede
                          + a_eq / 3. / (a + a_eq) / (a + a_eq);
        break;
    case ADE:
        *dw_over_da_fld = (3 * pba->a_ADE * pow(pba->a_ADE / a,
                                                -1 + 6 * (1 + pba->wf_ADE)) * pow(1 + pba->wf_ADE, 2)) / (a * a * pow(1 + pow(pba->a_ADE / a,
                                                        6 * (1 + pba->wf_ADE)), 1.5));
        break;
    }

    /** - finally, give the analytic solution of the following integral:
          \f$ \int_{a}^{a0} da 3(1+w_{fld})/a \f$. This is used in only
          one place, in the initial conditions for the background, and
          with a=a_ini. If your w(a) does not lead to a simple analytic
          solution of this integral, no worry: instead of writing
          something here, the best would then be to leave it equal to
          zero, and then in background_initial_conditions() you should
          implement a numerical calculation of this integral only for
          a=a_ini, using for instance Romberg integration. It should be
          fast, simple, and accurate enough. */
    switch(pba->fluid_equation_of_state)
    {
    case CLP:
        *integral_fld = 3.*((1. + pba->w0_fld + pba->wa_fld) * log(1. / a) + pba->wa_fld * (a - 1.));
        break;
    case EDE:
        class_stop(pba->error_message,
                   "EDE implementation not finished: to finish it, read the comments in background.c just before this line\n");
        break;
    case ADE:
        *integral_fld = 0; // numerically integrated in background_initial_conditions()
        break;
    }

    /** note: of course you can generalise these formulas to anything,
        defining new parameters pba->w..._fld. Just remember that so
        far, HyRec explicitely assumes that w(a)= w0 + wa (1-a/a0); but
        Recfast does not assume anything */

    return _SUCCESS_;
}

/**
 * Initialize the background structure, and in particular the
 * background interpolation table.
 *
 * @param ppr Input: pointer to precision structure
 * @param pba Input/Output: pointer to initialized background structure
 * @return the error status
 */

int background_init(
    struct precision * ppr,
    struct background * pba
)
{

    /** Summary: */

    /** - write class version */
    if(pba->background_verbose > 0)
    {
        printf("Running CLASS version %s\n", _VERSION_);
        printf("Computing background\n");
    }

    /** - if shooting failed during input, catch the error here */
    class_test_except(pba->shooting_failed == _TRUE_,
                      pba->error_message,
                      background_free_input(pba),
                      "Shooting failed, try optimising input_get_guess(). Error message:\n\n%s",
                      pba->shooting_error);

    /** - assign values to all indices in vectors of background quantities */
    /* Modification by Gen Ye: also initialize arrays freed in background_free_noinput here */
    class_call(background_indices(ppr, pba),
               pba->error_message,
               pba->error_message);

    /** - check that input parameters make sense and write additional information about them */
    class_call(background_checks(ppr, pba),
               pba->error_message,
               pba->error_message);

    /** - integrate the background over log(a), allocate and fill the background table */
    class_call(background_solve(ppr, pba),
               pba->error_message,
               pba->error_message);

    /** - find and store a few derived parameters at radiation-matter equality */
    class_call(background_find_equality(ppr, pba),
               pba->error_message,
               pba->error_message);

    /* - write a summary of the budget of the universe */
    class_call(background_output_budget(pba),
               pba->error_message,
               pba->error_message);

    return _SUCCESS_;

}

/**
 * Free all memory space allocated by background_init() and by input_read_parameters().
 *
 *
 * @param pba Input: pointer to background structure (to be freed)
 * @return the error status
 */

int background_free(
    struct background *pba
)
{

    class_call(background_free_noinput(pba),
               pba->error_message,
               pba->error_message);

    class_call(background_free_input(pba),
               pba->error_message,
               pba->error_message);

    return _SUCCESS_;
}

/**
 * Free only the memory space NOT allocated through
 * input_read_parameters(), but through background_init()
 *
 * @param pba Input: pointer to background structure (to be freed)
 * @return the error status
 */

int background_free_noinput(
    struct background *pba
)
{

    free(pba->tau_table);
    free(pba->z_table);
    free(pba->loga_table);
    free(pba->d2tau_dz2_table);
    free(pba->d2z_dtau2_table);
    free(pba->background_table);
    free(pba->d2background_dloga2_table);

    return _SUCCESS_;
}
/**
 * Free pointers inside background structure which were
 * allocated in input_read_parameters()
 *
 * @param pba Input: pointer to background structure
 * @return the error status
 */

int background_free_input(
    struct background *pba
)
{

    int k;

    if(pba->Omega0_ncdm_tot != 0.)
    {
        for(k = 0; k < pba->N_ncdm; k++)
        {
            free(pba->q_ncdm[k]);
            free(pba->w_ncdm[k]);
            free(pba->q_ncdm_bg[k]);
            free(pba->w_ncdm_bg[k]);
            free(pba->dlnf0_dlnq_ncdm[k]);
        }
        free(pba->ncdm_quadrature_strategy);
        free(pba->ncdm_input_q_size);
        free(pba->ncdm_qmax);
        free(pba->q_ncdm);
        free(pba->w_ncdm);
        free(pba->q_ncdm_bg);
        free(pba->w_ncdm_bg);
        free(pba->dlnf0_dlnq_ncdm);
        free(pba->q_size_ncdm);
        free(pba->q_size_ncdm_bg);
        free(pba->M_ncdm);
        free(pba->T_ncdm);
        free(pba->ksi_ncdm);
        free(pba->deg_ncdm);
        free(pba->Omega0_ncdm);
        free(pba->m_ncdm_in_eV);
        free(pba->factor_ncdm);
        if(pba->got_files != NULL)
            free(pba->got_files);
        if(pba->ncdm_psd_files != NULL)
            free(pba->ncdm_psd_files);
        if(pba->ncdm_psd_parameters != NULL)
            free(pba->ncdm_psd_parameters);
    }

    /* Free ndescf pointers, by Gen Ye */
    if(pba->N_ndescf > 0)
    {
        free(pba->Omega0_ndescf);
        free(pba->m_ndescf_in_eV);
        free(pba->ndescf_parameters);
        free(pba->ndescf_w_fld);
        free(pba->ndescf_m_fld);
        free(pba->ndescf_use_kg);
        free(pba->ndescf_loga_fld);
        free(pba->ndescf_tau_fld);
        free(pba->ndescf_tau_avgfld);
        free(pba->ndescf_sorted_index);
    }

    if(pba->scf_parameters != NULL) // modified so that EDE is also freed, by Gen Ye
        free(pba->scf_parameters);
    return _SUCCESS_;
}

/**
 * Assign value to each relevant index in vectors of background quantities.
 *
 * @param pba Input: pointer to background structure
 * @return the error status
 */

int background_indices(
    struct precision* ppr,
    struct background* pba
)
{

    /** Summary: */

    /** - define local variables */

    /* a running index for the vector of background quantities */
    int index_bg;
    /* a running index for the vector of background quantities to be integrated */
    int index_bi;

    /* indexing ndescf field, by Gen Ye */
    int counter;

    /** - initialize all flags: which species are present? */

    pba->has_cdm = _FALSE_;
    pba->has_ncdm = _FALSE_;
    pba->has_dcdm = _FALSE_;
    pba->has_dr = _FALSE_;
    pba->has_ndescf = _FALSE_; // ndescf flags, by Gen Ye
    pba->has_scf = _FALSE_;
    pba->has_lambda = _FALSE_;
    pba->has_fld = _FALSE_;
    pba->has_ur = _FALSE_;
    pba->has_idr = _FALSE_;
    pba->has_idm_dr = _FALSE_;
    pba->has_curvature = _FALSE_;

    if(pba->Omega0_cdm != 0.)
        pba->has_cdm = _TRUE_;

    if(pba->Omega0_ncdm_tot != 0.)
        pba->has_ncdm = _TRUE_;

    if(pba->Omega0_dcdmdr != 0.)
    {
        pba->has_dcdm = _TRUE_;
        if(pba->Gamma_dcdm != 0.)
            pba->has_dr = _TRUE_;
    }

    /* ndescf by Gen Ye */
    if(pba->N_ndescf > 0) // or do if Omega0_ndescf_tot > 0.?
    {
        pba->has_ndescf = _TRUE_;
        for(counter = 0; counter < pba->N_ndescf; counter++)
        {
            pba->ndescf_use_kg[counter] = _TRUE_;
        }
    }

    // force scf evolution if signaled, by Gen Ye
    if((pba->Omega0_scf != 0.) || (pba->force_ede == _TRUE_))
        pba->has_scf = _TRUE_;

    if(pba->Omega0_lambda != 0.)
        pba->has_lambda = _TRUE_;

    if(pba->Omega0_fld != 0.)
        pba->has_fld = _TRUE_;

    if(pba->Omega0_ur != 0.)
        pba->has_ur = _TRUE_;

    if(pba->Omega0_idr != 0.)
        pba->has_idr = _TRUE_;

    if(pba->Omega0_idm_dr != 0.)
        pba->has_idm_dr = _TRUE_;

    if(pba->sgnK != 0)
        pba->has_curvature = _TRUE_;

    /** - initialize all indices */

    index_bg = 0;

    /* index for scale factor */
    class_define_index(pba->index_bg_a, _TRUE_, index_bg, 1);

    /* - indices for H and its conformal-time-derivative */
    class_define_index(pba->index_bg_H, _TRUE_, index_bg, 1);
    class_define_index(pba->index_bg_H_prime, _TRUE_, index_bg, 1);

    /* - end of indices in the short vector of background values */
    pba->bg_size_short = index_bg;

    /* - index for rho_g (photon density) */
    class_define_index(pba->index_bg_rho_g, _TRUE_, index_bg, 1);

    /* - index for rho_b (baryon density) */
    class_define_index(pba->index_bg_rho_b, _TRUE_, index_bg, 1);

    /* - index for rho_cdm */
    class_define_index(pba->index_bg_rho_cdm, pba->has_cdm, index_bg, 1);

    /* - indices for ncdm. We only define the indices for ncdm1
       (density, pressure, pseudo-pressure), the other ncdm indices
       are contiguous */
    class_define_index(pba->index_bg_rho_ncdm1, pba->has_ncdm, index_bg, pba->N_ncdm);
    class_define_index(pba->index_bg_p_ncdm1, pba->has_ncdm, index_bg, pba->N_ncdm);
    class_define_index(pba->index_bg_pseudo_p_ncdm1, pba->has_ncdm, index_bg, pba->N_ncdm);

    /* - index for dcdm */
    class_define_index(pba->index_bg_rho_dcdm, pba->has_dcdm, index_bg, 1);

    /* - index for dr */
    class_define_index(pba->index_bg_rho_dr, pba->has_dr, index_bg, 1);

    /* - indices for scalar field */
    class_define_index(pba->index_bg_phi_scf, pba->has_scf, index_bg, 1);
    class_define_index(pba->index_bg_phi_prime_scf, pba->has_scf, index_bg, 1);
    class_define_index(pba->index_bg_V_scf, pba->has_scf, index_bg, 1);
    class_define_index(pba->index_bg_dV_scf, pba->has_scf, index_bg, 1);
    class_define_index(pba->index_bg_ddV_scf, pba->has_scf, index_bg, 1);
    class_define_index(pba->index_bg_rho_scf, pba->has_scf, index_bg, 1);
    class_define_index(pba->index_bg_p_scf, pba->has_scf, index_bg, 1);
    class_define_index(pba->index_bg_p_prime_scf, pba->has_scf, index_bg, 1);

    /* - index for Lambda */
    class_define_index(pba->index_bg_rho_lambda, pba->has_lambda, index_bg, 1);

    /* - index for fluid */
    class_define_index(pba->index_bg_rho_fld, pba->has_fld, index_bg, 1);
    class_define_index(pba->index_bg_w_fld, pba->has_fld, index_bg, 1);

    /* - index for ultra-relativistic neutrinos/species */
    class_define_index(pba->index_bg_rho_ur, pba->has_ur, index_bg, 1);

    /* - index for total density */
    class_define_index(pba->index_bg_rho_tot, _TRUE_, index_bg, 1);

    /* - index for total pressure */
    class_define_index(pba->index_bg_p_tot, _TRUE_, index_bg, 1);

    /* - index for derivative of total pressure */
    class_define_index(pba->index_bg_p_tot_prime, _TRUE_, index_bg, 1);

    /* - index for Omega_r (relativistic density fraction) */
    class_define_index(pba->index_bg_Omega_r, _TRUE_, index_bg, 1);

    /* - index interacting for dark radiation */
    class_define_index(pba->index_bg_rho_idr, pba->has_idr, index_bg, 1);

    /* - index for interacting dark matter */
    class_define_index(pba->index_bg_rho_idm_dr, pba->has_idm_dr, index_bg, 1);

    /* - put here additional ingredients that you want to appear in the
       normal vector */
    /* - indices for non dark energy scalar fields, by Gen Ye */
    class_define_index(pba->index_bg_phi_ndescf1, pba->has_ndescf, index_bg, pba->N_ndescf);
    class_define_index(pba->index_bg_phi_prime_ndescf1, pba->has_ndescf, index_bg, pba->N_ndescf);
    class_define_index(pba->index_bg_V_ndescf1, pba->has_ndescf, index_bg, pba->N_ndescf);
    class_define_index(pba->index_bg_dV_ndescf1, pba->has_ndescf, index_bg, pba->N_ndescf);
    class_define_index(pba->index_bg_ddV_ndescf1, pba->has_ndescf, index_bg, pba->N_ndescf);
    class_define_index(pba->index_bg_rho_ndescf1, pba->has_ndescf, index_bg, pba->N_ndescf);
    class_define_index(pba->index_bg_p_ndescf1, pba->has_ndescf, index_bg, pba->N_ndescf);
    class_define_index(pba->index_bg_p_prime_ndescf1, pba->has_ndescf, index_bg, pba->N_ndescf);
    /*    */

    /* - end of indices in the normal vector of background values */
    pba->bg_size_normal = index_bg;

    /* - indices in the long version : */

    /* -> critical density */
    class_define_index(pba->index_bg_rho_crit, _TRUE_, index_bg, 1);

    /* - index for Omega_m (non-relativistic density fraction) */
    class_define_index(pba->index_bg_Omega_m, _TRUE_, index_bg, 1);

    /* -> conformal distance */
    class_define_index(pba->index_bg_conf_distance, _TRUE_, index_bg, 1);

    /* -> angular diameter distance */
    class_define_index(pba->index_bg_ang_distance, _TRUE_, index_bg, 1);

    /* -> luminosity distance */
    class_define_index(pba->index_bg_lum_distance, _TRUE_, index_bg, 1);

    /* -> proper time (for age of the Universe) */
    class_define_index(pba->index_bg_time, _TRUE_, index_bg, 1);

    /* -> conformal sound horizon */
    class_define_index(pba->index_bg_rs, _TRUE_, index_bg, 1);

    /* -> density growth factor in dust universe */
    class_define_index(pba->index_bg_D, _TRUE_, index_bg, 1);

    /* -> velocity growth factor in dust universe */
    class_define_index(pba->index_bg_f, _TRUE_, index_bg, 1);

    /* -> put here additional quantities describing background */
    /*    */
    /*    */

    /* -> end of indices in the long vector of background values */
    pba->bg_size = index_bg;

    /* - now, indices in vector of variables to integrate.
       First {B} variables, then {C} variables. */

    index_bi = 0;

    /* -> index for conformal time in vector of variables to integrate */
    class_define_index(pba->index_bi_tau, _TRUE_, index_bi, 1);

    /* -> energy density in DCDM */
    class_define_index(pba->index_bi_rho_dcdm, pba->has_dcdm, index_bi, 1);

    /* -> energy density in DR */
    class_define_index(pba->index_bi_rho_dr, pba->has_dr, index_bi, 1);

    /* -> energy density in fluid */
    class_define_index(pba->index_bi_rho_fld, pba->has_fld, index_bi, 1);

    /* -> scalar field and its derivative wrt conformal time (Zuma) */
    class_define_index(pba->index_bi_phi_scf, pba->has_scf, index_bi, 1);
    class_define_index(pba->index_bi_phi_prime_scf, pba->has_scf, index_bi, 1);

    /* -> non dark energy scalar fields and its conformal derivative, by Gen Ye */
    class_define_index(pba->index_bi_phi_ndescf1, pba->has_ndescf, index_bi, pba->N_ndescf);
    class_define_index(pba->index_bi_phi_prime_ndescf1, pba->has_ndescf, index_bi, pba->N_ndescf);
    class_define_index(pba->index_bi_rho_ndescf1, pba->has_ndescf, index_bi, pba->N_ndescf);

    // log scale factor, needed in case of use adaptive step background sampling
    //class_define_index(pba->index_bi_loga, pba->adaptive_step_bk_sampling, index_bi, 1);
    /* end of modification by Gen Ye */

    /* End of {B} variables */
    pba->bi_B_size = index_bi;

    /* now continue with {C} variables */

    /* -> proper time (for age of the Universe) */
    class_define_index(pba->index_bi_time, _TRUE_, index_bi, 1);

    /* -> sound horizon */
    class_define_index(pba->index_bi_rs, _TRUE_, index_bi, 1);

    /* -> Second order equation for growth factor */
    class_define_index(pba->index_bi_D, _TRUE_, index_bi, 1);
    class_define_index(pba->index_bi_D_prime, _TRUE_, index_bi, 1);


    /* -> end of indices in the vector of variables to integrate */
    pba->bi_size = index_bi;

    /** - allocate background tables, moved here by Gen Ye */
    pba->bt_size = ppr->background_Nloga;

    class_alloc(pba->tau_table, pba->bt_size * sizeof(double), pba->error_message);
    class_alloc(pba->z_table, pba->bt_size * sizeof(double), pba->error_message);
    class_alloc(pba->loga_table, pba->bt_size * sizeof(double), pba->error_message);

    class_alloc(pba->d2tau_dz2_table, pba->bt_size * sizeof(double), pba->error_message);
    class_alloc(pba->d2z_dtau2_table, pba->bt_size * sizeof(double), pba->error_message);

    class_alloc(pba->background_table, pba->bt_size * pba->bg_size * sizeof(double),
                pba->error_message);
    class_alloc(pba->d2background_dloga2_table, pba->bt_size * pba->bg_size * sizeof(double),
                pba->error_message);

    return _SUCCESS_;

}

/**
 * This is the routine where the distribution function f0(q) of each
 * ncdm species is specified (it is the only place to modify if you
 * need a partlar f0(q))
 *
 * @param pbadist Input:  structure containing all parameters defining f0(q)
 * @param q       Input:  momentum
 * @param f0      Output: phase-space distribution
 */

int background_ncdm_distribution(
    void * pbadist,
    double q,
    double * f0
)
{
    struct background * pba;
    struct background_parameters_for_distributions * pbadist_local;
    int n_ncdm, lastidx;
    double ksi;
    double qlast, dqlast, f0last, df0last;
    double *param;
    /* Variables corresponding to entries in param: */
    //double square_s12,square_s23,square_s13;
    //double mixing_matrix[3][3];
    //int i;

    /** - extract from the input structure pbadist all the relevant information */
    pbadist_local = pbadist;          /* restore actual format of pbadist */
    pba = pbadist_local->pba;         /* extract the background structure from it */
    param = pba->ncdm_psd_parameters; /* extract the optional parameter list from it */
    n_ncdm = pbadist_local->n_ncdm;   /* extract index of ncdm species under consideration */
    ksi = pba->ksi_ncdm[n_ncdm];      /* extract chemical potential */

    /** - shall we interpolate in file, or shall we use analytical formula below? */

    /** - a) deal first with the case of interpolating in files */
    if(pba->got_files[n_ncdm] == _TRUE_)
    {

        lastidx = pbadist_local->tablesize - 1;
        if(q < pbadist_local->q[0])
        {
            //Handle q->0 case:
            *f0 = pbadist_local->f0[0];
        }
        else if(q > pbadist_local->q[lastidx])
        {
            //Handle q>qmax case (ensure continuous and derivable function with Boltzmann tail):
            qlast = pbadist_local->q[lastidx];
            f0last = pbadist_local->f0[lastidx];
            dqlast = qlast - pbadist_local->q[lastidx - 1];
            df0last = f0last - pbadist_local->f0[lastidx - 1];

            *f0 = f0last * exp(-(qlast - q) * df0last / f0last / dqlast);
        }
        else
        {
            //Do interpolation:
            class_call(array_interpolate_spline(
                           pbadist_local->q,
                           pbadist_local->tablesize,
                           pbadist_local->f0,
                           pbadist_local->d2f0,
                           1,
                           q,
                           &pbadist_local->last_index,
                           f0,
                           1,
                           pba->error_message),
                       pba->error_message,     pba->error_message);
        }
    }

    /** - b) deal now with case of reading analytical function */
    else
    {
        /**
           Next enter your analytic expression(s) for the p.s.d.'s. If
           you need different p.s.d.'s for different species, put each
           p.s.d inside a condition, like for instance: if (n_ncdm==2)
           {*f0=...}.  Remember that n_ncdm = 0 refers to the first
           species.
        */

        /**************************************************/
        /*    FERMI-DIRAC INCLUDING CHEMICAL POTENTIALS   */
        /**************************************************/

        *f0 = 1.0 / pow(2 * _PI_, 3) * (1. / (exp(q - ksi) + 1.) + 1. / (exp(q + ksi) + 1.));

        /**************************************************/

        /** This form is only appropriate for approximate studies, since in
            reality the chemical potentials are associated with flavor
            eigenstates, not mass eigenstates. It is easy to take this into
            account by introducing the mixing angles. In the later part
            (not read by the code) we illustrate how to do this. */

        if(_FALSE_)
        {

            /* We must use the list of extra parameters read in input, stored in the
               ncdm_psd_parameter list, extracted above from the structure
               and now called param[..] */

            /* check that this list has been read */
            class_test(param == NULL,
                       pba->error_message,
                       "Analytic expression wants to use 'ncdm_psd_parameters', but they have not been entered!");

            /* extract values from the list (in this example, mixing angles) */
            double square_s12 = param[0];
            double square_s23 = param[1];
            double square_s13 = param[2];

            /* infer mixing matrix */
            double mixing_matrix[3][3];
            int i;

            mixing_matrix[0][0] = pow(fabs(sqrt((1 - square_s12) * (1 - square_s13))), 2);
            mixing_matrix[0][1] = pow(fabs(sqrt(square_s12 * (1 - square_s13))), 2);
            mixing_matrix[0][2] = fabs(square_s13);
            mixing_matrix[1][0] = pow(fabs(sqrt((1 - square_s12) * square_s13 * square_s23) + sqrt(
                                               square_s12 * (1 - square_s23))), 2);
            mixing_matrix[1][1] = pow(fabs(sqrt(square_s12 * square_s23 * square_s13) - sqrt((
                                               1 - square_s12) * (1 - square_s23))), 2);
            mixing_matrix[1][2] = pow(fabs(sqrt(square_s23 * (1 - square_s13))), 2);
            mixing_matrix[2][0] = pow(fabs(sqrt(square_s12 * square_s23) - sqrt((1 - square_s12) *
                                           square_s13 * (1 - square_s23))), 2);
            mixing_matrix[2][1] = pow(sqrt((1 - square_s12) * square_s23) + sqrt(square_s12 * square_s13 *
                                      (1 - square_s23)), 2);
            mixing_matrix[2][2] = pow(fabs(sqrt((1 - square_s13) * (1 - square_s23))), 2);

            /* loop over flavor eigenstates and compute psd of mass eigenstates */
            *f0 = 0.0;
            for(i = 0; i < 3; i++)
            {

                *f0 += mixing_matrix[i][n_ncdm] * 1.0 / pow(2 * _PI_,
                        3) * (1. / (exp(q - pba->ksi_ncdm[i]) + 1.) + 1. / (exp(q + pba->ksi_ncdm[i]) + 1.));

            }
        } /* end of region not used, but shown as an example */
    }

    return _SUCCESS_;
}

/**
 * This function is only used for the purpose of finding optimal
 * quadrature weights. The logic is: if we can accurately convolve
 * f0(q) with this function, then we can convolve it accurately with
 * any other relevant function.
 *
 * @param pbadist Input:  structure containing all background parameters
 * @param q       Input:  momentum
 * @param test    Output: value of the test function test(q)
 */

int background_ncdm_test_function(
    void * pbadist,
    double q,
    double * test
)
{

    double c = 2.0 / (3.0 * _zeta3_);
    double d = 120.0 / (7.0 * pow(_PI_, 4));
    double e = 2.0 / (45.0 * _zeta5_);

    /** Using a + bq creates problems for otherwise acceptable distributions
        which diverges as \f$ 1/r \f$ or \f$ 1/r^2 \f$ for \f$ r\to 0 \f$*/
    *test = pow(2.0 * _PI_, 3) / 6.0 * (c * q * q - d * q * q * q - e * q * q * q * q);

    return _SUCCESS_;
}

/**
 * This function finds optimal quadrature weights for each ncdm
 * species
 *
 * @param ppr Input: precision structure
 * @param pba Input/Output: background structure
 */

int background_ncdm_init(
    struct precision *ppr,
    struct background *pba
)
{

    int index_q, k, tolexp, row, status, filenum;
    double f0m2, f0m1, f0, f0p1, f0p2, dq, q, df0dq, tmp1, tmp2;
    struct background_parameters_for_distributions pbadist;
    FILE *psdfile;

    pbadist.pba = pba;

    /* Allocate pointer arrays: */
    class_alloc(pba->q_ncdm, sizeof(double*)*pba->N_ncdm, pba->error_message);
    class_alloc(pba->w_ncdm, sizeof(double*)*pba->N_ncdm, pba->error_message);
    class_alloc(pba->q_ncdm_bg, sizeof(double*)*pba->N_ncdm, pba->error_message);
    class_alloc(pba->w_ncdm_bg, sizeof(double*)*pba->N_ncdm, pba->error_message);
    class_alloc(pba->dlnf0_dlnq_ncdm, sizeof(double*)*pba->N_ncdm, pba->error_message);

    /* Allocate pointers: */
    class_alloc(pba->q_size_ncdm, sizeof(int)*pba->N_ncdm, pba->error_message);
    class_alloc(pba->q_size_ncdm_bg, sizeof(int)*pba->N_ncdm, pba->error_message);
    class_alloc(pba->factor_ncdm, sizeof(double)*pba->N_ncdm, pba->error_message);

    for(k = 0, filenum = 0; k < pba->N_ncdm; k++)
    {
        pbadist.n_ncdm = k;
        pbadist.q = NULL;
        pbadist.tablesize = 0;
        /*Do we need to read in a file to interpolate the distribution function? */
        if((pba->got_files != NULL) && (pba->got_files[k] == _TRUE_))
        {
            psdfile = fopen(pba->ncdm_psd_files + filenum * _ARGUMENT_LENGTH_MAX_, "r");
            class_test(psdfile == NULL, pba->error_message,
                       "Could not open file %s!", pba->ncdm_psd_files + filenum * _ARGUMENT_LENGTH_MAX_);
            // Find size of table:
            for(row = 0, status = 2; status == 2; row++)
            {
                status = fscanf(psdfile, "%lf %lf", &tmp1, &tmp2);
            }
            rewind(psdfile);
            pbadist.tablesize = row - 1;

            /*Allocate room for interpolation table: */
            class_alloc(pbadist.q, sizeof(double)*pbadist.tablesize, pba->error_message);
            class_alloc(pbadist.f0, sizeof(double)*pbadist.tablesize, pba->error_message);
            class_alloc(pbadist.d2f0, sizeof(double)*pbadist.tablesize, pba->error_message);
            for(row = 0; row < pbadist.tablesize; row++)
            {
                status = fscanf(psdfile, "%lf %lf",
                                &pbadist.q[row], &pbadist.f0[row]);
                //        printf("(q,f0) = (%g,%g)\n",pbadist.q[row],pbadist.f0[row]);
            }
            fclose(psdfile);
            /* Call spline interpolation: */
            class_call(array_spline_table_lines(pbadist.q,
                                                pbadist.tablesize,
                                                pbadist.f0,
                                                1,
                                                pbadist.d2f0,
                                                _SPLINE_EST_DERIV_,
                                                pba->error_message),
                       pba->error_message,
                       pba->error_message);
            filenum++;
        }

        /* Handle perturbation qsampling: */
        if(pba->ncdm_quadrature_strategy[k] == qm_auto)
        {
            /** Automatic q-sampling for this species */
            class_alloc(pba->q_ncdm[k], _QUADRATURE_MAX_ * sizeof(double), pba->error_message);
            class_alloc(pba->w_ncdm[k], _QUADRATURE_MAX_ * sizeof(double), pba->error_message);

            class_call(get_qsampling(pba->q_ncdm[k],
                                     pba->w_ncdm[k],
                                     &(pba->q_size_ncdm[k]),
                                     _QUADRATURE_MAX_,
                                     ppr->tol_ncdm,
                                     pbadist.q,
                                     pbadist.tablesize,
                                     background_ncdm_test_function,
                                     background_ncdm_distribution,
                                     &pbadist,
                                     pba->error_message),
                       pba->error_message,
                       pba->error_message);
            pba->q_ncdm[k] = realloc(pba->q_ncdm[k], pba->q_size_ncdm[k] * sizeof(double));
            pba->w_ncdm[k] = realloc(pba->w_ncdm[k], pba->q_size_ncdm[k] * sizeof(double));


            if(pba->background_verbose > 0)
            {
                printf("ncdm species i=%d sampled with %d points for purpose of perturbation integration\n",
                       k + 1,
                       pba->q_size_ncdm[k]);
            }

            /* Handle background q_sampling: */
            class_alloc(pba->q_ncdm_bg[k], _QUADRATURE_MAX_BG_ * sizeof(double), pba->error_message);
            class_alloc(pba->w_ncdm_bg[k], _QUADRATURE_MAX_BG_ * sizeof(double), pba->error_message);

            class_call(get_qsampling(pba->q_ncdm_bg[k],
                                     pba->w_ncdm_bg[k],
                                     &(pba->q_size_ncdm_bg[k]),
                                     _QUADRATURE_MAX_BG_,
                                     ppr->tol_ncdm_bg,
                                     pbadist.q,
                                     pbadist.tablesize,
                                     background_ncdm_test_function,
                                     background_ncdm_distribution,
                                     &pbadist,
                                     pba->error_message),
                       pba->error_message,
                       pba->error_message);

            pba->q_ncdm_bg[k] = realloc(pba->q_ncdm_bg[k], pba->q_size_ncdm_bg[k] * sizeof(double));
            pba->w_ncdm_bg[k] = realloc(pba->w_ncdm_bg[k], pba->q_size_ncdm_bg[k] * sizeof(double));

            /** - in verbose mode, inform user of number of sampled momenta
                for background quantities */
            if(pba->background_verbose > 0)
            {
                printf("ncdm species i=%d sampled with %d points for purpose of background integration\n",
                       k + 1,
                       pba->q_size_ncdm_bg[k]);
            }
        }
        else
        {
            /** Manual q-sampling for this species. Same sampling used for both perturbation and background sampling, since this will usually be a high precision setting anyway */
            pba->q_size_ncdm_bg[k] = pba->ncdm_input_q_size[k];
            pba->q_size_ncdm[k] = pba->ncdm_input_q_size[k];
            class_alloc(pba->q_ncdm_bg[k], pba->q_size_ncdm_bg[k]*sizeof(double), pba->error_message);
            class_alloc(pba->w_ncdm_bg[k], pba->q_size_ncdm_bg[k]*sizeof(double), pba->error_message);
            class_alloc(pba->q_ncdm[k], pba->q_size_ncdm[k]*sizeof(double), pba->error_message);
            class_alloc(pba->w_ncdm[k], pba->q_size_ncdm[k]*sizeof(double), pba->error_message);
            class_call(get_qsampling_manual(pba->q_ncdm[k],
                                            pba->w_ncdm[k],
                                            pba->q_size_ncdm[k],
                                            pba->ncdm_qmax[k],
                                            pba->ncdm_quadrature_strategy[k],
                                            pbadist.q,
                                            pbadist.tablesize,
                                            background_ncdm_distribution,
                                            &pbadist,
                                            pba->error_message),
                       pba->error_message,
                       pba->error_message);
            for(index_q = 0; index_q < pba->q_size_ncdm[k]; index_q++)
            {
                pba->q_ncdm_bg[k][index_q] = pba->q_ncdm[k][index_q];
                pba->w_ncdm_bg[k][index_q] = pba->w_ncdm[k][index_q];
            }
            /** - in verbose mode, inform user of number of sampled momenta
                for background quantities */
            if(pba->background_verbose > 0)
            {
                printf("ncdm species i=%d sampled with %d points for purpose of background andperturbation integration using the manual method\n",
                       k + 1,
                       pba->q_size_ncdm[k]);
            }
        }

        class_alloc(pba->dlnf0_dlnq_ncdm[k],
                    pba->q_size_ncdm[k]*sizeof(double),
                    pba->error_message);


        for(index_q = 0; index_q < pba->q_size_ncdm[k]; index_q++)
        {
            q = pba->q_ncdm[k][index_q];
            class_call(background_ncdm_distribution(&pbadist, q, &f0),
                       pba->error_message, pba->error_message);

            //Loop to find appropriate dq:
            for(tolexp = _PSD_DERIVATIVE_EXP_MIN_; tolexp < _PSD_DERIVATIVE_EXP_MAX_; tolexp++)
            {

                if(index_q == 0)
                {
                    dq = MIN((0.5 - ppr->smallest_allowed_variation) * q,
                             2 * exp(tolexp) * (pba->q_ncdm[k][index_q + 1] - q));
                }
                else if(index_q == pba->q_size_ncdm[k] - 1)
                {
                    dq = exp(tolexp) * 2.0 * (pba->q_ncdm[k][index_q] - pba->q_ncdm[k][index_q - 1]);
                }
                else
                {
                    dq = exp(tolexp) * (pba->q_ncdm[k][index_q + 1] - pba->q_ncdm[k][index_q - 1]);
                }

                class_call(background_ncdm_distribution(&pbadist, q - 2 * dq, &f0m2),
                           pba->error_message, pba->error_message);
                class_call(background_ncdm_distribution(&pbadist, q + 2 * dq, &f0p2),
                           pba->error_message, pba->error_message);

                if(fabs((f0p2 - f0m2) / f0) > sqrt(ppr->smallest_allowed_variation)) break;
            }

            class_call(background_ncdm_distribution(&pbadist, q - dq, &f0m1),
                       pba->error_message, pba->error_message);
            class_call(background_ncdm_distribution(&pbadist, q + dq, &f0p1),
                       pba->error_message, pba->error_message);
            //5 point estimate of the derivative:
            df0dq = (+f0m2 - 8 * f0m1 + 8 * f0p1 - f0p2) / 12.0 / dq;
            //printf("df0dq[%g] = %g. dlf=%g ?= %g. f0 =%g.\n",q,df0dq,q/f0*df0dq,
            //Avoid underflow in extreme tail:
            if(fabs(f0) == 0.)
                pba->dlnf0_dlnq_ncdm[k][index_q] = -q; /* valid for whatever f0 with exponential tail in exp(-q) */
            else
                pba->dlnf0_dlnq_ncdm[k][index_q] = q / f0 * df0dq;
        }

        pba->factor_ncdm[k] = pba->deg_ncdm[k] * 4 * _PI_ * pow(pba->T_cmb * pba->T_ncdm[k] * _k_B_,
                              4) * 8 * _PI_ * _G_
                              / 3. / pow(_h_P_ / 2. / _PI_, 3) / pow(_c_, 7) * _Mpc_over_m_ * _Mpc_over_m_;

        /* If allocated, deallocate interpolation table:  */
        if((pba->got_files != NULL) && (pba->got_files[k] == _TRUE_))
        {
            free(pbadist.q);
            free(pbadist.f0);
            free(pbadist.d2f0);
        }
    }


    return _SUCCESS_;
}

/**
 * For a given ncdm species: given the quadrature weights, the mass
 * and the redshift, find background quantities by a quick weighted
 * sum over.  Input parameters passed as NULL pointers are not
 * evaluated for speed-up
 *
 * @param qvec     Input: sampled momenta
 * @param wvec     Input: quadrature weights
 * @param qsize    Input: number of momenta/weights
 * @param M        Input: mass
 * @param factor   Input: normalization factor for the p.s.d.
 * @param z        Input: redshift
 * @param n        Output: number density
 * @param rho      Output: energy density
 * @param p        Output: pressure
 * @param drho_dM  Output: derivative used in next function
 * @param pseudo_p Output: pseudo-pressure used in perturbation module for fluid approx
 *
 */

int background_ncdm_momenta(
    /* Only calculate for non-NULL pointers: */
    double * qvec,
    double * wvec,
    int qsize,
    double M,
    double factor,
    double z,
    double * n,
    double * rho, // density
    double * p,   // pressure
    double * drho_dM,  // d rho / d M used in next function
    double * pseudo_p  // pseudo-p used in ncdm fluid approx
)
{

    int index_q;
    double epsilon;
    double q2;
    double factor2;
    /** Summary: */

    /** - rescale normalization at given redshift */
    factor2 = factor * pow(1 + z, 4);

    /** - initialize quantities */
    if(n != NULL) *n = 0.;
    if(rho != NULL) *rho = 0.;
    if(p != NULL) *p = 0.;
    if(drho_dM != NULL) *drho_dM = 0.;
    if(pseudo_p != NULL) *pseudo_p = 0.;

    /** - loop over momenta */
    for(index_q = 0; index_q < qsize; index_q++)
    {

        /* squared momentum */
        q2 = qvec[index_q] * qvec[index_q];

        /* energy */
        epsilon = sqrt(q2 + M * M / (1. + z) / (1. + z));

        /* integrand of the various quantities */
        if(n != NULL) *n += q2 * wvec[index_q];
        if(rho != NULL) *rho += q2 * epsilon * wvec[index_q];
        if(p != NULL) *p += q2 * q2 / 3. / epsilon * wvec[index_q];
        if(drho_dM != NULL) *drho_dM += q2 * M / (1. + z) / (1. + z) / epsilon * wvec[index_q];
        if(pseudo_p != NULL) *pseudo_p += pow(q2 / epsilon, 3) / 3.0 * wvec[index_q];
    }

    /** - adjust normalization */
    if(n != NULL) *n *= factor2 / (1. + z);
    if(rho != NULL) *rho *= factor2;
    if(p != NULL) *p *= factor2;
    if(drho_dM != NULL) *drho_dM *= factor2;
    if(pseudo_p != NULL) *pseudo_p *= factor2;

    return _SUCCESS_;
}

/**
 * When the user passed the density fraction Omega_ncdm or
 * omega_ncdm in input but not the mass, infer the mass with Newton iteration method.
 *
 * @param ppr    Input: precision structure
 * @param pba    Input/Output: background structure
 * @param n_ncdm Input: index of ncdm species
 */

int background_ncdm_M_from_Omega(
    struct precision *ppr,
    struct background *pba,
    int n_ncdm
)
{
    double rho0, rho, n, M, deltaM, drhodM;
    int iter, maxiter = 50;

    rho0 = pba->H0 * pba->H0 *
           pba->Omega0_ncdm[n_ncdm]; /*Remember that rho is defined such that H^2=sum(rho_i) */
    M = 0.0;

    background_ncdm_momenta(pba->q_ncdm_bg[n_ncdm],
                            pba->w_ncdm_bg[n_ncdm],
                            pba->q_size_ncdm_bg[n_ncdm],
                            M,
                            pba->factor_ncdm[n_ncdm],
                            0.,
                            &n,
                            &rho,
                            NULL,
                            NULL,
                            NULL);

    /* Is the value of Omega less than a massless species?*/
    class_test(rho0 < rho, pba->error_message,
               "The value of Omega for the %dth species, %g, is less than for a massless species! It should be atleast %g. Check your input.",
               n_ncdm, pba->Omega0_ncdm[n_ncdm], pba->Omega0_ncdm[n_ncdm]*rho / rho0);

    /* In the strict NR limit we have rho = n*(M) today, giving a zeroth order guess: */
    M = rho0 / n; /* This is our guess for M. */
    for(iter = 1; iter <= maxiter; iter++)
    {

        /* Newton iteration. First get relevant quantities at M: */
        background_ncdm_momenta(pba->q_ncdm_bg[n_ncdm],
                                pba->w_ncdm_bg[n_ncdm],
                                pba->q_size_ncdm_bg[n_ncdm],
                                M,
                                pba->factor_ncdm[n_ncdm],
                                0.,
                                NULL,
                                &rho,
                                NULL,
                                &drhodM,
                                NULL);

        deltaM = (rho0 - rho) / drhodM; /* By definition of the derivative */
        if((M + deltaM) < 0.0) deltaM = -M / 2.0; /* Avoid overshooting to negative M value. */
        M += deltaM; /* Update value of M.. */
        if(fabs(deltaM / M) < ppr->tol_M_ncdm)
        {
            /* Accuracy reached.. */
            pba->M_ncdm[n_ncdm] = M;
            break;
        }
    }
    class_test(iter >= maxiter, pba->error_message,
               "Newton iteration could not converge on a mass for some reason.");
    return _SUCCESS_;
}

/**
 * Perform some check on the input background quantities, and send to
 * standard output some information about them
 *
 * @param ppr Input: pointer to precision structure
 * @param pba Input: pointer to initialized background structure
 * @return the error status
*/

int background_checks(
    struct precision* ppr,
    struct background* pba
)
{

    /** - define local variables */
    int n_ncdm;
    double rho_ncdm_rel, rho_nu_rel;
    double N_dark;
    double w_fld, dw_over_da, integral_fld;
    int filenum = 0;

    /** - control that cosmological parameter values make sense, otherwise inform user */

    /* H0 in Mpc^{-1} */
    /* Many users asked for this test to be supressed. It is commented out. */
    /*class_test((pba->H0 < _H0_SMALL_)||(pba->H0 > _H0_BIG_),
               pba->error_message,
               "H0=%g out of bounds (%g<H0<%g) \n",pba->H0,_H0_SMALL_,_H0_BIG_);*/

    /* consistency between h and H0 */
    class_test(fabs(pba->h * 1.e5 / _c_  / pba->H0 - 1.) > ppr->smallest_allowed_variation,
               pba->error_message,
               "inconsistency between Hubble and reduced Hubble parameters: you have H0=%f/Mpc=%fkm/s/Mpc, but h=%f",
               pba->H0, pba->H0 / 1.e5 * _c_, pba->h);

    /* T_cmb in K */
    /* Many users asked for this test to be supressed. It is commented out. */
    /*class_test((pba->T_cmb < _TCMB_SMALL_)||(pba->T_cmb > _TCMB_BIG_),
               pba->error_message,
               "T_cmb=%g out of bounds (%g<T_cmb<%g)",pba->T_cmb,_TCMB_SMALL_,_TCMB_BIG_);*/

    /* Omega_k */
    /* Many users asked for this test to be supressed. It is commented out. */
    /*class_test((pba->Omega0_k < _OMEGAK_SMALL_)||(pba->Omega0_k > _OMEGAK_BIG_),
               pba->error_message,
               "Omegak = %g out of bounds (%g<Omegak<%g) \n",pba->Omega0_k,_OMEGAK_SMALL_,_OMEGAK_BIG_);*/

    /* fluid equation of state */
    if(pba->has_fld == _TRUE_)
    {

        class_call(background_w_fld(pba, 0., &w_fld, &dw_over_da, &integral_fld), pba->error_message,
                   pba->error_message);

        class_test(w_fld >= 1. / 3.,
                   pba->error_message,
                   "Your choice for w(a--->0)=%g is suspicious, since it is bigger than -1/3 there cannot be radiation domination at early times\n",
                   w_fld);
    }

    /** - in verbose mode, send to standard output some additional information on non-obvious background parameters */
    if(pba->background_verbose > 0)
    {

        if(pba->has_ncdm == _TRUE_)
        {

            /* loop over ncdm species */
            for(n_ncdm = 0; n_ncdm < pba->N_ncdm; n_ncdm++)
            {

                /* inform if p-s-d read in files */
                if(pba->got_files[n_ncdm] == _TRUE_)
                {
                    printf(" -> ncdm species i=%d read from file %s\n", n_ncdm + 1,
                           pba->ncdm_psd_files + filenum * _ARGUMENT_LENGTH_MAX_);
                    filenum++;
                }

                /* inform the user also about the value of the ncdm
                   masses in eV and about */
                printf(" -> non-cold dark matter species with i=%d has m_i = %e eV (so m_i / omega_i =%e eV)\n",
                       n_ncdm + 1,
                       pba->m_ncdm_in_eV[n_ncdm],
                       pba->m_ncdm_in_eV[n_ncdm]*pba->deg_ncdm[n_ncdm] / pba->Omega0_ncdm[n_ncdm] / pba->h / pba->h);

                /* call this function to get rho_ncdm */
                background_ncdm_momenta(pba->q_ncdm_bg[n_ncdm],
                                        pba->w_ncdm_bg[n_ncdm],
                                        pba->q_size_ncdm_bg[n_ncdm],
                                        0.,
                                        pba->factor_ncdm[n_ncdm],
                                        0.,
                                        NULL,
                                        &rho_ncdm_rel,
                                        NULL,
                                        NULL,
                                        NULL);

                /* inform user of the contribution of each species to
                   radiation density (in relativistic limit): should be
                   between 1.01 and 1.02 for each active neutrino species;
                   evaluated as rho_ncdm/rho_nu_rel where rho_nu_rel is the
                   density of one neutrino in the instantaneous decoupling
                   limit, i.e. assuming T_nu=(4/11)^1/3 T_gamma (this comes
                   from the definition of N_eff) */
                rho_nu_rel = 56.0 / 45.0 * pow(_PI_, 6) * pow(4.0 / 11.0, 4.0 / 3.0) * _G_ / pow(_h_P_,
                             3) / pow(_c_, 7) *
                             pow(_Mpc_over_m_, 2) * pow(pba->T_cmb * _k_B_, 4);

                printf(" -> ncdm species i=%d sampled with %d (resp. %d) points for purpose of background (resp. perturbation) integration. In the relativistic limit it gives Delta N_eff = %g\n",
                       n_ncdm + 1,
                       pba->q_size_ncdm_bg[n_ncdm],
                       pba->q_size_ncdm[n_ncdm],
                       rho_ncdm_rel / rho_nu_rel);
            }
        }

        /* contribution of interacting dark radiation _idr to N_eff */
        if(pba->has_idr == _TRUE_)
        {
            N_dark = pba->Omega0_idr / 7.*8. / pow(4. / 11., 4. / 3.) / pba->Omega0_g;
            printf(" -> dark radiation Delta Neff %e\n", N_dark);
        }
    }

    return _SUCCESS_;
}

/**
 *  This function integrates the background over time, allocates and
 *  fills the background table
 *
 * @param ppr Input: precision structure
 * @param pba Input/Output: background structure
 */

int background_solve(
    struct precision *ppr,
    struct background *pba
)
{

    /** Summary: */

    /* in case of ndescf, rapid oscillation requires adaptive step size
     * as in CLASS 2.x we use generic integrator to perform integration step by step
     * function modified by Gen Ye */

    /** - define local variables */

    /* contains all quantities relevant for the integration algorithm */
    struct generic_integrator_workspace gi;
    /* parameters and workspace for the background_derivs function */
    struct background_parameters_and_workspace bpaw;
    /* a growing table (since the number of time steps is not known a priori) */
    growTable gTable;
    /* needed for growing table */
    double * pData;
    /* needed for growing table */
    void * memcopy_result;
    /* vector of quantities to be integrated */
    double * pvecback_integration;
    /* vector of all background quantities */
    double * pvecback;
    /* comoving radius coordinate in Mpc (equal to conformal distance in flat case) */
    double comoving_radius = 0.;
    /* conformal distance in Mpc (equal to comoving radius in flat case) */
    double conformal_distance;

    /* evolvers */
    extern int evolver_rk();
    extern int evolver_ndf15();
    int (*generic_evolver)() = evolver_ndf15;

    /* initial and final loga values */
    double loga_ini, loga_final;
    /* growth factor today */
    double D_today;
    /* indices for the different arrays */
    int index_loga, index_scf;
    /* what parameters are used in the output? */
    int * used_in_output;

    /* params for ndescf, by Gen Ye */
    int counter, sorted_idx_ndescf;
    double a, loga_start, loga_end, timescale, bk_integration_stepsize;
    double * loga_gt;

    /** - setup background workspace */
    bpaw.pba = pba;
    class_alloc(pvecback, pba->bg_size * sizeof(double), pba->error_message);
    bpaw.pvecback = pvecback;

    /** - allocate vector of quantities to be integrated */
    class_alloc(pvecback_integration, pba->bi_size * sizeof(double), pba->error_message);

    /** - impose initial conditions with background_initial_conditions() */
    class_call(background_initial_conditions(ppr, pba, pvecback, pvecback_integration, &(loga_ini)),
               pba->error_message,
               pba->error_message);

    /* with our conventions, loga is in fact log(a/a_0); we integrate until today, when log(a/a_0) = 0 */
    loga_final = 0.;

    if(pba->adaptive_step_bk_sampling == _TRUE_) // step by step integration added by Gen Ye
    {
        /** - initialize generic integrator with initialize_generic_integrator() */
        class_call(initialize_generic_integrator(pba->bi_size, &gi),
                   gi.error_message,
                   pba->error_message);

        /** - create a growTable with gt_init() */
        class_call(gt_init(&gTable),
                   gTable.error_message,
                   pba->error_message);

        /** - alloc array to temperaly store loga */
        class_alloc(loga_gt, sizeof(double) ,pba->error_message);

        /* each step start from loga_end of the previous step
         * thus we initially loga_end is the initial time */
        loga_end = loga_ini;

        /* in this case we do not know loga table in priori,
         * it is thus updated in pvecback_integration with index pba->index_bi_loga
         * and then stored in the growing table */
        //pvecback_integration[pba->index_bi_loga] = loga_ini;

        pba->bt_size = 0;

        /* points to the ndescf we are about to switch on fld approx */
        counter = 0;
        if(pba->has_ndescf == _TRUE_)
        {
            while((counter<pba->N_ndescf)&&(pba->ndescf_use_kg[pba->ndescf_sorted_index[counter]] == _FALSE_))
            {
                pba->ndescf_tau_fld[pba->ndescf_sorted_index[counter]] = 0.;
                counter++;
            }
        }

        /** - loop over integration steps: call background_functions(),
            find step size, save data in growTable with gt_add(),
            perform one step with generic_integrator(), store new value of loga */
        while(loga_end < loga_final)
        {
            /* each step start from loga_end of the previous step */
            loga_start = loga_end;

            /* -> find step size:
             *     1. call background_functions() to evaluate quantities may be needed
             *        note rho_ndescf is updated in this function, rather than in the integrator, if KG is solved
             *     2. call background_timescale to determine the minimal physical time scale
             *     3. check if we hit the ndescf fld approximation switches
             *     4. determine step size based on the minimal time scale */

            class_call(background_functions(pba, exp(loga_start), pvecback_integration, short_info, pvecback),
                       pba->error_message,
                       pba->error_message);

            class_call(background_timescale(loga_start, &bpaw, &(timescale), pba->error_message),
                       pba->error_message,
                       pba->error_message);

            bk_integration_stepsize = ppr->background_integration_stepsize * timescale;

            //printf("%g ",timescale);

            if(loga_start + bk_integration_stepsize + ppr->smallest_allowed_variation < loga_final)
            {
                loga_end = loga_start + bk_integration_stepsize;
            }
            else
            {
                loga_end = loga_final;
            }

            if((pba->has_ndescf == _TRUE_) && (counter < pba->N_ndescf))
            {
                sorted_idx_ndescf = pba->ndescf_sorted_index[counter];
                if(pba->ndescf_loga_fld[sorted_idx_ndescf] < loga_end+ppr->smallest_allowed_variation)
                    loga_end = pba->ndescf_loga_fld[sorted_idx_ndescf];
            }

            class_test(loga_end - loga_start < ppr->smallest_allowed_variation,
                       pba->error_message,
                       "integration step: change in loga =%e < machine precision : leads either to numerical error or infinite loop",
                       (loga_end - loga_start));

            /* -> save data in growTable */
            class_call(gt_add(&gTable, _GT_END_, (void *) pvecback_integration, sizeof(double)*pba->bi_size),
                       gTable.error_message,
                       pba->error_message);
            /* -> save log scale factor in growTable */
//            class_call(gt_add(&gTable, _GT_END_, (void *) &loga_start, sizeof(double)),
//                       gTable.error_message,
//                       pba->error_message);
            class_realloc(loga_gt, loga_gt, sizeof(double)*(pba->bt_size+1), pba->error_message);
            loga_gt[pba->bt_size] = loga_start;
            pba->bt_size++;

            /* -> perform one step */
            class_call_except(generic_integrator(background_derivs,
                                          loga_start,
                                          loga_end,
                                          pvecback_integration,
                                          &bpaw,
                                          ppr->tol_background_integration,
                                          ppr->smallest_allowed_variation,
                                          &gi),
                       gi.error_message,
                       pba->error_message,
                       // free local arrays upon failure to avoid memory leak in MCMC
                       free(pvecback);
                       free(pvecback_integration);
                       free(loga_gt);
                       cleanup_generic_integrator(&gi);
                       gt_free(&gTable);
                       // free the background module, it's not freed elsewhere if it returned error
                       background_free(pba));

            /*-> update the log scale factor here */
            //pvecback_integration[pba->index_bi_loga] = loga_end;

            /* switch on ndescf fld approx if conditons met */
            if((pba->has_ndescf == _TRUE_) && (counter < pba->N_ndescf))
            {
                if(loga_end == pba->ndescf_loga_fld[sorted_idx_ndescf])
                {
                    pba->ndescf_use_kg[sorted_idx_ndescf] = _FALSE_;
                    pba->ndescf_tau_fld[sorted_idx_ndescf] = pvecback_integration[pba->index_bi_tau];
                    // Begin preparing ndescf fld perturbation initial condition some before switching to fld approx
                    pba->ndescf_tau_avgfld[sorted_idx_ndescf] = pvecback_integration[pba->index_bi_tau] - 1.4 * _PI_ /
                            exp(loga_end) / pba->ndescf_m_fld[sorted_idx_ndescf];
                    counter++;
                }
            }
        }

        /** - save last data in growTable with gt_add() */
        class_call(gt_add(&gTable, _GT_END_, (void *) pvecback_integration, sizeof(double)*pba->bi_size),
                   gTable.error_message,
                   pba->error_message);
        class_realloc(loga_gt, loga_gt, sizeof(double)*(pba->bt_size+1), pba->error_message);
        loga_gt[pba->bt_size] = loga_final;
        pba->bt_size++;

        /** - clean up generic integrator with cleanup_generic_integrator() */
        class_call(cleanup_generic_integrator(&gi),
                   gi.error_message,
                   pba->error_message);

        /** - retrieve data stored in the growTable with gt_getPtr() */
        class_call(gt_getPtr(&gTable, (void**)&pData),
                   gTable.error_message,
                   pba->error_message);

        /** - reallocate background tables */
        class_realloc(pba->tau_table, pba->tau_table, pba->bt_size * sizeof(double), pba->error_message);
        class_realloc(pba->z_table, pba->z_table, pba->bt_size * sizeof(double), pba->error_message);
        class_realloc(pba->loga_table, pba->loga_table, pba->bt_size * sizeof(double), pba->error_message);

        class_realloc(pba->d2tau_dz2_table, pba->d2tau_dz2_table, pba->bt_size * sizeof(double), pba->error_message);
        class_realloc(pba->d2z_dtau2_table, pba->d2z_dtau2_table, pba->bt_size * sizeof(double), pba->error_message);

        class_realloc(pba->background_table, pba->background_table, pba->bt_size * pba->bg_size * sizeof(double),
                    pba->error_message);
        class_realloc(pba->d2background_dloga2_table, pba->d2background_dloga2_table, pba->bt_size * pba->bg_size * sizeof(double),
                    pba->error_message);

        /* never used in this case */
        class_alloc(used_in_output, pba->bt_size * sizeof(int), pba->error_message);
    }
    else
    {
        /** This is the original integration method */

        /* initialization of background tables moved to background_indices() by Gen Ye */

        class_alloc(used_in_output, pba->bt_size * sizeof(int), pba->error_message);

        /** - define values of loga at which results will be stored */
        for(index_loga = 0; index_loga < pba->bt_size; index_loga++)
        {
            pba->loga_table[index_loga] = loga_ini + index_loga * (loga_final - loga_ini) / (pba->bt_size - 1);
            used_in_output[index_loga] = 1;
        }

        /** - choose the right evolver */
        switch(ppr->background_evolver)
        {

        case rk:
            generic_evolver = evolver_rk;
            if(pba->background_verbose > 1)
            {
                printf("%s\n", "Chose rk as generic_evolver");
            }
            break;

        case ndf15:
            generic_evolver = evolver_ndf15;
            if(pba->background_verbose > 1)
            {
                printf("%s\n", "Chose ndf15 as generic_evolver");
            }
            break;
        }

        /** Perform piecewise integration with ndescf fld approx switches, by Gen Ye
          * same as vanilla CLASS when ndescf is absent */
        loga_start = loga_ini;
        loga_end = loga_ini;
        counter = 0;
        while(loga_start != loga_final)
        {
            if((pba->has_ndescf == _TRUE_) && (counter < pba->N_ndescf))
            {
                sorted_idx_ndescf = pba->ndescf_sorted_index[counter];
                if(pba->ndescf_loga_fld[sorted_idx_ndescf] < loga_start)
                {
                    //printf("Scalar field Nr. %d starts as fluid.\n". pba->ndescf_sorted_index[index_ndescf]);
                    pba->ndescf_use_kg[sorted_idx_ndescf] = _FALSE_;
                    pba->ndescf_tau_fld[sorted_idx_ndescf] = 0.;
                    counter ++;
                    continue;
                }
                else if(pba->ndescf_loga_fld[sorted_idx_ndescf] < loga_final)
                    loga_end = pba->ndescf_loga_fld[sorted_idx_ndescf];
                else loga_end = loga_final;
            }
            else
            {
                loga_end = loga_final;
            }
            //printf("start = %g,end=%g,use KG=%d\n",loga_s,loga_e,pba->ndescf_use_kg[pba->ndescf_sorted_index[index_ndescf]]);
            class_call_except(generic_evolver(background_derivs,
                                       loga_start,
                                       loga_end,
                                       pvecback_integration,
                                       used_in_output,
                                       pba->bi_size,
                                       &bpaw,
                                       ppr->tol_background_integration,
                                       ppr->smallest_allowed_variation,
                                       background_timescale, //'evaluate_timescale', required by evolver_rk but not by ndf15
                                       ppr->background_integration_stepsize,
                                       pba->loga_table,
                                       pba->bt_size,
                                       background_sources,
                                       NULL, //'print_variables' in evolver_rk could be set, but, not required
                                       pba->error_message),
                       pba->error_message,
                       pba->error_message,
                       // free local arrays if failed, to avoid memory leak in MCMC, by Gen Ye
                       free(pvecback);
                       free(pvecback_integration);
                       free(used_in_output);
                       // free the background module, it's not freed elsewhere if it returned error by Gen Ye
                       background_free(pba));
            if((pba->has_ndescf == _TRUE_) && (counter < pba->N_ndescf))
            {
                pba->ndescf_use_kg[sorted_idx_ndescf] = _FALSE_;
                pba->ndescf_tau_fld[sorted_idx_ndescf] = pvecback_integration[pba->index_bi_tau];
                // Begin preparing fld initial condition one oscillation before
                pba->ndescf_tau_avgfld[sorted_idx_ndescf] = pvecback_integration[pba->index_bi_tau] - 2 * _PI_ /
                        exp(loga_end) / pba->ndescf_m_fld[sorted_idx_ndescf];
            }
            loga_start = loga_end;
            counter ++;
        }
        /* end of modification by Gen Ye */
    }

    /** - recover some quantities today */
    /* -> age in Gyears */
    pba->age = pvecback_integration[pba->index_bi_time] / _Gyr_over_Mpc_;
    /* -> conformal age in Mpc */
    pba->conformal_age = pvecback_integration[pba->index_bi_tau];
    /* -> contribution of decaying dark matter and dark radiation to the critical density today: */
    if(pba->has_dcdm == _TRUE_)
    {
        pba->Omega0_dcdm = pvecback_integration[pba->index_bi_rho_dcdm] / pba->H0 / pba->H0;
    }
    if(pba->has_dr == _TRUE_)
    {
        pba->Omega0_dr = pvecback_integration[pba->index_bi_rho_dr] / pba->H0 / pba->H0;
    }
    /* contribution of ndescf to the critical density today, by Gen Ye */
    if(pba->has_ndescf == _TRUE_)
    {
        pba->Omega0_ndescf_tot_real = 0.;
        for(counter = 0; counter < pba->N_ndescf; counter++)
        {
            // ndescf_sorted_index should also sort ndescf_tau_fld, if not then something went wrong
            if(counter > 0)
            {
                class_test(pba->ndescf_tau_fld[pba->ndescf_sorted_index[counter - 1]] >
                           pba->ndescf_tau_fld[pba->ndescf_sorted_index[counter]],
                           pba->error_message,
                           "Something went wrong in the background module. ndescf_sorted_index does not sort ndescf_tau_fld.");
            }
            pba->Omega0_ndescf[counter] = pvecback_integration[pba->index_bi_rho_ndescf1 + counter] /
                                          pba->H0 / pba->H0 ;
            pba->Omega0_ndescf_tot_real += pba->Omega0_ndescf[counter];
            pba->m_ndescf_in_eV[counter] = pow(ddV_ndescf(pba, 0.,counter),0.5)/1.5638e29;
        }
    }
    /* -> scale-invariant growth rate today */
    D_today = pvecback_integration[pba->index_bi_D];

    /** - In a loop over lines, fill rest of background table for
          quantities that depend on numbers like "conformal_age" or
          "D_today" that were calculated just before */
    for(index_loga = 0; index_loga < pba->bt_size; index_loga++)
    {

        if(pba->adaptive_step_bk_sampling == _TRUE_)
        {
            /* Extract data from the growing table in case of step by step integration */
            // it is (pba->bi_size+1) because we also store loga
            //pba->tau_table[index_loga] = pData[index_loga * (pba->bi_size+1) + pba->index_bi_tau];
            pba->tau_table[index_loga] = pData[index_loga * pba->bi_size + pba->index_bi_tau];

            //pba->loga_table[index_loga] = pData[index_loga * (pba->bi_size+2) -1];
            pba->loga_table[index_loga] = loga_gt[index_loga];

            a = exp(pba->loga_table[index_loga]);

            pba->z_table[index_loga] = MAX(0., 1/a - 1.);

            /* -> compute all other quantities depending only on {B} variables.
               The value of {B} variables in pData are also copied to pvecback.*/
//            class_call(background_functions(pba, a, pData + index_loga * (pba->bi_size+1), long_info, pvecback),
//                       pba->error_message,
//                       pba->error_message);
            class_call(background_functions(pba, a, pData + index_loga * pba->bi_size, long_info, pvecback),
                       pba->error_message,
                       pba->error_message);

            /* -> write in the table */
            memcopy_result = memcpy(pba->background_table + index_loga * pba->bg_size, pvecback,
                                    pba->bg_size * sizeof(double));

            class_test(memcopy_result != pba->background_table + index_loga * pba->bg_size,
                       pba->error_message,
                       "cannot copy data back to pba->background_table");
        }

        pba->background_table[index_loga * pba->bg_size + pba->index_bg_D] *= 1. / D_today;

        conformal_distance = pba->conformal_age - pba->tau_table[index_loga];
        pba->background_table[index_loga * pba->bg_size + pba->index_bg_conf_distance] = conformal_distance;

        if(pba->sgnK == 0)
        {
            comoving_radius = conformal_distance;
        }
        else if(pba->sgnK == 1)
        {
            comoving_radius = sin(sqrt(pba->K) * conformal_distance) / sqrt(pba->K);
        }
        else if(pba->sgnK == -1)
        {
            comoving_radius = sinh(sqrt(-pba->K) * conformal_distance) / sqrt(-pba->K);
        }

        pba->background_table[index_loga * pba->bg_size + pba->index_bg_ang_distance] = comoving_radius /
                (1. + pba->z_table[index_loga]);
        pba->background_table[index_loga * pba->bg_size + pba->index_bg_lum_distance] = comoving_radius *
                (1. + pba->z_table[index_loga]);
    }

    /** - free the growTable with gt_free() by Gen Ye */
    if(pba->adaptive_step_bk_sampling == _TRUE_)
    {
        class_call(gt_free(&gTable),
                   gTable.error_message,
                   pba->error_message);

        free(loga_gt);
    }

    /** - fill tables of second derivatives (in view of spline interpolation) */
    class_call(array_spline_table_lines(pba->z_table,
                                        pba->bt_size,
                                        pba->tau_table,
                                        1,
                                        pba->d2tau_dz2_table,
                                        _SPLINE_EST_DERIV_,
                                        pba->error_message),
               pba->error_message,
               pba->error_message);

    class_call(array_spline_table_lines(pba->tau_table,
                                        pba->bt_size,
                                        pba->z_table,
                                        1,
                                        pba->d2z_dtau2_table,
                                        _SPLINE_EST_DERIV_,
                                        pba->error_message),
               pba->error_message,
               pba->error_message);

    class_call(array_spline_table_lines(pba->loga_table,
                                        pba->bt_size,
                                        pba->background_table,
                                        pba->bg_size,
                                        pba->d2background_dloga2_table,
                                        _SPLINE_EST_DERIV_,
                                        pba->error_message),
               pba->error_message,
               pba->error_message);

    /** - compute remaining "related parameters" */

    /**  - so-called "effective neutrino number", computed at earliest
        time in interpolation table. This should be seen as a
        definition: Neff is the equivalent number of
        instantaneously-decoupled neutrinos accounting for the
        radiation density, beyond photons */

    pba->Neff = (pba->background_table[pba->index_bg_Omega_r]
                 * pba->background_table[pba->index_bg_rho_crit]
                 - pba->background_table[pba->index_bg_rho_g])
                / (7. / 8.*pow(4. / 11., 4. / 3.) * pba->background_table[pba->index_bg_rho_g]);

    /** - send information to standard output */
    if(pba->background_verbose > 0)
    {
        printf(" -> age = %f Gyr\n", pba->age);
        printf(" -> conformal age = %f Mpc\n", pba->conformal_age);
        printf(" -> N_eff = %g (summed over all species that are non-relativistic at early times) \n",
               pba->Neff);
    }

    if(pba->background_verbose > 2)
    {
        if((pba->has_dcdm == _TRUE_) && (pba->has_dr == _TRUE_))
        {
            printf("    Decaying Cold Dark Matter details: (DCDM --> DR)\n");
            printf("     -> Omega0_dcdm = %f\n", pba->Omega0_dcdm);
            printf("     -> Omega0_dr = %f\n", pba->Omega0_dr);
            printf("     -> Omega0_dr+Omega0_dcdm = %f, input value = %f\n",
                   pba->Omega0_dr + pba->Omega0_dcdm, pba->Omega0_dcdmdr);
            printf("     -> Omega_ini_dcdm/Omega_b = %f\n", pba->Omega_ini_dcdm / pba->Omega0_b);
        }
        if(pba->has_scf == _TRUE_)
        {
            printf("    Scalar field details:\n");
            printf("     -> Omega_scf = %g, wished %g\n",
                   pba->background_table[(pba->bt_size - 1)*pba->bg_size + pba->index_bg_rho_scf] /
                   pba->background_table[(pba->bt_size - 1)*pba->bg_size + pba->index_bg_rho_crit], pba->Omega0_scf);
            if(pba->has_lambda == _TRUE_)
            {
                printf("     -> Omega_Lambda = %g, wished %g\n",
                       pba->background_table[(pba->bt_size - 1)*pba->bg_size + pba->index_bg_rho_lambda] /
                       pba->background_table[(pba->bt_size - 1)*pba->bg_size + pba->index_bg_rho_crit],
                       pba->Omega0_lambda);
            }
//            printf("     -> parameters: [lambda, alpha, A, B] = \n");
            printf("     -> parameters: [");
            for(index_scf = 0; index_scf < pba->scf_parameters_size - 1; index_scf++)
            {
                printf("%.3f, ", pba->scf_parameters[index_scf]);
            }
            printf("%.3f]\n", pba->scf_parameters[pba->scf_parameters_size - 1]);
        }
        /* ndescf details by Gen Ye */
        if(pba->has_ndescf == _TRUE_)
        {
            printf("    Additional scalar fields details:\n");
            printf("     -> Shooting guess Omega0_ndescf_tot = %g, actual evolution yields Omega0_ndescf_tot = %g.\n", pba->Omega0_ndescf_tot, pba->Omega0_ndescf_tot_real);
            for(counter = 0; counter < pba->N_ndescf; counter++)
            {
                index_scf = counter * pba->ndescf_parameters_size;
                printf("     -> Field Nr.%d:\n", counter);
                if(pba->ndescf_parameters[index_scf + pba->ndescf_parameters_size - 3]<0)
                {
                    printf("        - f_a = %.2e Mp = %.2e GeV\n", pba->ndescf_parameters[index_scf], pba->ndescf_parameters[index_scf]*2.43534e18);
                    printf("        - m_a = %.2e H0 = %.2e eV\n", pow(ddV_ndescf(pba, 0.,counter),0.5)/pba->H0, pba->m_ndescf_in_eV[counter]);
                }
                printf("        - Omega0 = %f, takes up %.2f%% of total matter\n", pba->Omega0_ndescf[counter], 100*pba->Omega0_ndescf[counter]/pba->background_table[(pba->bt_size - 1) * pba->bg_size + pba->index_bg_Omega_m]);
                if(pba->ndescf_tau_fld[counter]==0.)
                {
                    printf("        - Evolves as fluid with w = %.1f throughout\n", pba->ndescf_w_fld[counter]);
                }
                else if(pba->ndescf_tau_fld[counter]>pba->conformal_age)
                {
                    printf("        - Evolves as scalar field throughout\n");
                }
                else
                {
                    printf("        - Fluid approximation starts at z = %.1f, tau = %.1f Mpc\n", exp(-pba->ndescf_loga_fld[counter])-1, pba->ndescf_tau_fld[counter]);
                }
            }
        }
    }

    /**  - store information in the background structure */
    pba->Omega0_m = pba->background_table[(pba->bt_size - 1) * pba->bg_size + pba->index_bg_Omega_m];
    pba->Omega0_r = pba->background_table[(pba->bt_size - 1) * pba->bg_size + pba->index_bg_Omega_r];
    pba->Omega0_de = 1. - (pba->Omega0_m + pba->Omega0_r + pba->Omega0_k);

    free(pvecback);
    free(pvecback_integration);
    free(used_in_output);

    return _SUCCESS_;

}

/**
 * Assign initial values to background integrated variables.
 *
 * @param ppr                  Input: pointer to precision structure
 * @param pba                  Input: pointer to background structure
 * @param pvecback             Input: vector of background quantities used as workspace
 * @param pvecback_integration Output: vector of background quantities to be integrated, returned with proper initial values
 * @param loga_ini             Output: value of loga (in fact with our conventions log(a/a_0)) at initial time
 * @return the error status
 */

int background_initial_conditions(
    struct precision *ppr,
    struct background *pba,
    double * pvecback, /* vector with argument pvecback[index_bg] (must be already allocated, normal format is sufficient) */
    double * pvecback_integration, /* vector with argument pvecback_integration[index_bi] (must be already allocated with size pba->bi_size) */
    double * loga_ini
)
{

    /** Summary: */

    /** - define local variables */

    /* scale factor */
    double a;

    double rho_ncdm, p_ncdm, rho_ncdm_rel_tot = 0.;
    double f, Omega_rad, rho_rad;
    int counter, is_early_enough, n_ncdm;
    int index_ndescf = 0; // index for ndescf, by Gen Ye
    double phi, phi_prime; // variables for ndescf, by Gen Ye
    double scf_lambda;
    double rho_fld_today;
    double w_fld, dw_over_da_fld, integral_fld;

    /** - fix initial value of \f$ a \f$ */
    a = ppr->a_ini_over_a_today_default;

    /**  If we have ncdm species, perhaps we need to start earlier
         than the standard value for the species to be relativistic.
         This could happen for some WDM models.
    */

    if(pba->has_ncdm == _TRUE_)
    {

        for(counter = 0; counter < _MAX_IT_; counter++)
        {

            is_early_enough = _TRUE_;
            rho_ncdm_rel_tot = 0.;

            for(n_ncdm = 0; n_ncdm < pba->N_ncdm; n_ncdm++)
            {

                class_call(background_ncdm_momenta(pba->q_ncdm_bg[n_ncdm],
                                                   pba->w_ncdm_bg[n_ncdm],
                                                   pba->q_size_ncdm_bg[n_ncdm],
                                                   pba->M_ncdm[n_ncdm],
                                                   pba->factor_ncdm[n_ncdm],
                                                   1. / a - 1.0,
                                                   NULL,
                                                   &rho_ncdm,
                                                   &p_ncdm,
                                                   NULL,
                                                   NULL),
                           pba->error_message,
                           pba->error_message);
                rho_ncdm_rel_tot += 3.*p_ncdm;
                if(fabs(p_ncdm / rho_ncdm - 1. / 3.) > ppr->tol_ncdm_initial_w)
                {
                    is_early_enough = _FALSE_;
                }
            }
            if(is_early_enough == _TRUE_)
            {
                break;
            }
            else
            {
                a *= _SCALE_BACK_;
            }
        }
        class_test(counter == _MAX_IT_,
                   pba->error_message,
                   "Search for initial scale factor a such that all ncdm species are relativistic failed.");
    }

    /* Set initial values of {B} variables: */
    Omega_rad = pba->Omega0_g;
    if(pba->has_ur == _TRUE_)
    {
        Omega_rad += pba->Omega0_ur;
    }
    if(pba->has_idr == _TRUE_)
    {
        Omega_rad += pba->Omega0_idr;
    }
    rho_rad = Omega_rad * pow(pba->H0, 2) / pow(a, 4);
    if(pba->has_ncdm == _TRUE_)
    {
        /** - We must add the relativistic contribution from NCDM species */
        rho_rad += rho_ncdm_rel_tot;
    }
    if(pba->has_dcdm == _TRUE_)
    {
        /* Remember that the critical density today in CLASS conventions is H0^2 */
        pvecback_integration[pba->index_bi_rho_dcdm] =
            pba->Omega_ini_dcdm * pba->H0 * pba->H0 * pow(a, -3);
        if(pba->background_verbose > 3)
            printf("Density is %g. Omega_ini=%g\n", pvecback_integration[pba->index_bi_rho_dcdm],
                   pba->Omega_ini_dcdm);
    }

    if(pba->has_dr == _TRUE_)
    {
        if(pba->has_dcdm == _TRUE_)
        {
            /**  - f is the critical density fraction of DR. The exact solution is:
             *
             * `f = -Omega_rad+pow(pow(Omega_rad,3./2.)+0.5*pow(a,6)*pvecback_integration[pba->index_bi_rho_dcdm]*pba->Gamma_dcdm/pow(pba->H0,3),2./3.);`
             *
             * but it is not numerically stable for very small f which is always the case.
             * Instead we use the Taylor expansion of this equation, which is equivalent to
             * ignoring f(a) in the Hubble rate.
             */
            f = 1. / 3.*pow(a, 6) * pvecback_integration[pba->index_bi_rho_dcdm] * pba->Gamma_dcdm / pow(
                    pba->H0, 3) / sqrt(Omega_rad);
            pvecback_integration[pba->index_bi_rho_dr] = f * pba->H0 * pba->H0 / pow(a, 4);
        }
        else
        {
            /** There is also a space reserved for a future case where dr is not sourced by dcdm */
            pvecback_integration[pba->index_bi_rho_dr] = 0.0;
        }
    }

    if(pba->has_fld == _TRUE_)
    {

        /* rho_fld today */
        rho_fld_today = pba->Omega0_fld * pow(pba->H0, 2);

        /* integrate rho_fld(a) from a_ini to a_0, to get rho_fld(a_ini) given rho_fld(a0) */
        class_call(background_w_fld(pba, a, &w_fld, &dw_over_da_fld, &integral_fld), pba->error_message,
                   pba->error_message);

        /* Note: for complicated w_fld(a) functions with no simple
           analytic integral, this is the place were you should compute
           numerically the simple 1d integral [int_{a_ini}^{a_0} 3
           [(1+w_fld)/a] da] (e.g. with the Romberg method?) instead of
           calling background_w_fld */

        /* rho_fld at initial time */
        pvecback_integration[pba->index_bi_rho_fld] = rho_fld_today * exp(integral_fld);

    }

    /** - Fix initial value of \f$ \phi, \phi' \f$
     * set directly in the radiation attractor => fixes the units in terms of rho_ur
     *
     * TODO:
     * - There seems to be some small oscillation when it starts.
     * - Check equations and signs. Sign of phi_prime?
     * - is rho_ur all there is early on?
     */
    if(pba->has_scf == _TRUE_)
    {
        scf_lambda = pba->scf_parameters[0];
        if(pba->attractor_ic_scf == _TRUE_)
        {
            pvecback_integration[pba->index_bi_phi_scf] = -1 / scf_lambda *
                    log(rho_rad * 4. / (3 * pow(scf_lambda, 2) - 12)) * pba->phi_ini_scf;
            if(3.*pow(scf_lambda, 2) - 12. < 0)
            {
                /** - --> If there is no attractor solution for scf_lambda, assign some value. Otherwise would give a nan.*/
                pvecback_integration[pba->index_bi_phi_scf] = 1. / scf_lambda; //seems to do the work
                if(pba->background_verbose > 0)
                {
                    printf(" No attractor IC for lambda = %.3e ! \n ", scf_lambda);
                }
            }
            pvecback_integration[pba->index_bi_phi_prime_scf] = 2.*a * sqrt(V_scf(pba,
                    pvecback_integration[pba->index_bi_phi_scf])) * pba->phi_prime_ini_scf;
        }
        else
        {
            if(pba->background_verbose > 0)
            {
                printf("Not using attractor initial conditions\n");
            }
            /** - --> If no attractor initial conditions are assigned, gets the provided ones. */
            pvecback_integration[pba->index_bi_phi_scf] = pba->phi_ini_scf;
            pvecback_integration[pba->index_bi_phi_prime_scf] = pba->phi_prime_ini_scf;
        }
        class_test(!isfinite(pvecback_integration[pba->index_bi_phi_scf]) ||
                   !isfinite(pvecback_integration[pba->index_bi_phi_scf]),
                   pba->error_message,
                   "initial phi = %e phi_prime = %e -> check initial conditions",
                   pvecback_integration[pba->index_bi_phi_scf],
                   pvecback_integration[pba->index_bi_phi_scf]);
    }

    /** Set non dark energy scalar field initial conditions.
      * Currently must be specified by the last two entries of ndescf_parameters in the input file for each field. by Gen Ye*/
    if(pba->has_ndescf == _TRUE_)
    {
        for(counter = 0; counter < pba->N_ndescf; counter++)
        {
            index_ndescf += pba->ndescf_parameters_size;
            phi = pba->ndescf_parameters[index_ndescf - 2];
            phi_prime = pba->ndescf_parameters[index_ndescf - 1];
            pvecback_integration[pba->index_bi_phi_ndescf1 + counter] = phi;
            pvecback_integration[pba->index_bi_phi_prime_ndescf1 + counter] = phi_prime;
            pvecback_integration[pba->index_bi_rho_ndescf1 + counter] = phi_prime * phi_prime /
                    (2 * a * a) + V_ndescf(pba, phi, counter);
            if(log(a)>pba->ndescf_loga_fld[counter])
            {
                pba->ndescf_use_kg[counter] = _FALSE_;
            }
            else
            {
                pba->ndescf_use_kg[counter] = _TRUE_;
            }
        }
    }

    /* Infer pvecback from pvecback_integration */
    class_call(background_functions(pba, a, pvecback_integration, normal_info, pvecback),
               pba->error_message,
               pba->error_message);

    /* Just checking that our initial time indeed is deep enough in the radiation
       dominated regime */
    class_test(fabs(pvecback[pba->index_bg_Omega_r] - 1.) > ppr->tol_initial_Omega_r,
               pba->error_message,
               "Omega_r = %e, not close enough to 1. Decrease a_ini_over_a_today_default in order to start from radiation domination.",
               pvecback[pba->index_bg_Omega_r]);

    /** - compute initial proper time, assuming radiation-dominated
        universe since Big Bang and therefore \f$ t=1/(2H) \f$ (good
        approximation for most purposes) */

    class_test(pvecback[pba->index_bg_H] <= 0.,
               pba->error_message,
               "H = %e instead of strictly positive", pvecback[pba->index_bg_H]);

    pvecback_integration[pba->index_bi_time] = 1. / (2.* pvecback[pba->index_bg_H]);

    /** - compute initial conformal time, assuming radiation-dominated
        universe since Big Bang and therefore \f$ \tau=1/(aH) \f$
        (good approximation for most purposes) */
    pvecback_integration[pba->index_bi_tau] = 1. / (a * pvecback[pba->index_bg_H]);

    /** - compute initial sound horizon, assuming \f$ c_s=1/\sqrt{3} \f$ initially */
    pvecback_integration[pba->index_bi_rs] = pvecback_integration[pba->index_bi_tau] / sqrt(3.);

    /** - set initial value of D and D' in RD. D and D' need only be set up to an overall constant, since they will later be re-normalized. From Ma&Bertschinger, one can derive D ~ (ktau)^2 at early times, from which one finds D'/D = 2 aH (assuming aH=1/tau during RD) */
    pvecback_integration[pba->index_bi_D] = 1.;
    pvecback_integration[pba->index_bi_D_prime] = 2.*a * pvecback[pba->index_bg_H];

    /** - return the value finally chosen for the initial log(a) */
    *loga_ini = log(a);

    return _SUCCESS_;

}

/**
 * Find the time of radiation/matter equality and store characteristic
 * quantitites at that time in the background structure..
 *
 * @param ppr                  Input: pointer to precision structure
 * @param pba                  Input/Output: pointer to background structure
 * @return the error status
 */

int background_find_equality(
    struct precision *ppr,
    struct background *pba
)
{

    double Omega_m_over_Omega_r = 0.;
    int index_tau_minus = 0;
    int index_tau_plus = pba->bt_size - 1;
    int index_tau_mid = 0;
    double tau_minus, tau_plus, tau_mid = 0.;
    double * pvecback;

    /* first bracket the right tau value between two consecutive indices in the table */

    while((index_tau_plus - index_tau_minus) > 1)
    {

        index_tau_mid = (int)(0.5 * (index_tau_plus + index_tau_minus));

        Omega_m_over_Omega_r = pba->background_table[index_tau_mid * pba->bg_size + pba->index_bg_Omega_m]
                               / pba->background_table[index_tau_mid * pba->bg_size + pba->index_bg_Omega_r];

        if(Omega_m_over_Omega_r > 1)
            index_tau_plus = index_tau_mid;
        else
            index_tau_minus = index_tau_mid;

    }

    /* then get a better estimate within this range */

    tau_minus = pba->tau_table[index_tau_minus];
    tau_plus =  pba->tau_table[index_tau_plus];

    class_alloc(pvecback, pba->bg_size * sizeof(double), pba->error_message);

    while((tau_plus - tau_minus) > ppr->tol_tau_eq)
    {

        tau_mid = 0.5 * (tau_plus + tau_minus);

        class_call(background_at_tau(pba, tau_mid, long_info, inter_closeby, &index_tau_minus, pvecback),
                   pba->error_message,
                   pba->error_message);

        Omega_m_over_Omega_r = pvecback[pba->index_bg_Omega_m] / pvecback[pba->index_bg_Omega_r];

        if(Omega_m_over_Omega_r > 1)
            tau_plus = tau_mid;
        else
            tau_minus = tau_mid;

    }

    pba->a_eq = pvecback[pba->index_bg_a];
    pba->H_eq = pvecback[pba->index_bg_H];
    pba->z_eq = 1. / pba->a_eq - 1.;
    pba->tau_eq = tau_mid;

    if(pba->background_verbose > 0)
    {
        printf(" -> radiation/matter equality at z = %f\n", pba->z_eq);
        printf("    corresponding to conformal time = %f Mpc\n", pba->tau_eq);
    }

    free(pvecback);

    return _SUCCESS_;

}


/**
 * Subroutine for formatting background output
 *
 * @param pba                  Input: pointer to background structure
 * @param titles               Ouput: name of columns when printing the background table
 * @return the error status
 */

int background_output_titles(
    struct background * pba,
    char titles[_MAXTITLESTRINGLENGTH_]
)
{

    /** - Length of the column title should be less than _OUTPUTPRECISION_+6
        to be indented correctly, but it can be as long as . */
    int n;
    char tmp[40];

    class_store_columntitle(titles, "z", _TRUE_);
    class_store_columntitle(titles, "proper time [Gyr]", _TRUE_);
    class_store_columntitle(titles, "conf. time [Mpc]", _TRUE_);
    class_store_columntitle(titles, "H [1/Mpc]", _TRUE_);
    class_store_columntitle(titles, "comov. dist.", _TRUE_);
    class_store_columntitle(titles, "ang.diam.dist.", _TRUE_);
    class_store_columntitle(titles, "lum. dist.", _TRUE_);
    class_store_columntitle(titles, "comov.snd.hrz.", _TRUE_);
    class_store_columntitle(titles, "(.)rho_g", _TRUE_);
    class_store_columntitle(titles, "(.)rho_b", _TRUE_);
    class_store_columntitle(titles, "(.)rho_cdm", pba->has_cdm);
    if(pba->has_ncdm == _TRUE_)
    {
        for(n = 0; n < pba->N_ncdm; n++)
        {
            sprintf(tmp, "(.)rho_ncdm[%d]", n);
            class_store_columntitle(titles, tmp, _TRUE_);
            sprintf(tmp, "(.)p_ncdm[%d]", n);
            class_store_columntitle(titles, tmp, _TRUE_);
        }
    }
    class_store_columntitle(titles, "(.)rho_lambda", pba->has_lambda);
    class_store_columntitle(titles, "(.)rho_fld", pba->has_fld);
    class_store_columntitle(titles, "(.)w_fld", pba->has_fld);
    class_store_columntitle(titles, "(.)rho_ur", pba->has_ur);
    class_store_columntitle(titles, "(.)rho_idr", pba->has_idr);
    class_store_columntitle(titles, "(.)rho_idm_dr", pba->has_idm_dr);
    class_store_columntitle(titles, "(.)rho_crit", _TRUE_);
    class_store_columntitle(titles, "(.)rho_dcdm", pba->has_dcdm);
    class_store_columntitle(titles, "(.)rho_dr", pba->has_dr);

    /* Non dark energy scalar fields quantities, by Gen Ye */
    if(pba->has_ndescf == _TRUE_)
    {
        for(n = 0; n < pba->N_ndescf; n++)
        {
            sprintf(tmp, "(.)rho_ndescf[%d]", n);
            class_store_columntitle(titles, tmp, _TRUE_);
            sprintf(tmp, "(.)p_ndescf[%d]", n);
            class_store_columntitle(titles, tmp, _TRUE_);
            sprintf(tmp, "phi_ndescf[%d]", n);
            class_store_columntitle(titles, tmp, _TRUE_);
            sprintf(tmp, "phi'_ndescf[%d]", n);
            class_store_columntitle(titles, tmp, _TRUE_);
        }
    }

    class_store_columntitle(titles, "(.)rho_scf", pba->has_scf);
    class_store_columntitle(titles, "(.)p_scf", pba->has_scf);
    class_store_columntitle(titles, "(.)p_prime_scf", pba->has_scf);
    class_store_columntitle(titles, "phi_scf", pba->has_scf);
    class_store_columntitle(titles, "phi'_scf", pba->has_scf);
    class_store_columntitle(titles, "V_scf", pba->has_scf);
    class_store_columntitle(titles, "V'_scf", pba->has_scf);
    class_store_columntitle(titles, "V''_scf", pba->has_scf);

    class_store_columntitle(titles, "(.)rho_tot", _TRUE_);
    class_store_columntitle(titles, "(.)p_tot", _TRUE_);
    class_store_columntitle(titles, "(.)p_tot_prime", _TRUE_);

    class_store_columntitle(titles, "gr.fac. D", _TRUE_);
    class_store_columntitle(titles, "gr.fac. f", _TRUE_);

    return _SUCCESS_;
}

/**
 * Subroutine for writing the background output
 *
 * @param pba                  Input: pointer to background structure
 * @param number_of_titles     Input: number of background quantities to print at each time step
 * @param data                 Ouput: 1d array storing all the background table
 * @return the error status
 */

int background_output_data(
    struct background *pba,
    int number_of_titles,
    double *data
)
{

    int index_tau, storeidx, n;
    double *dataptr, *pvecback;

    /** Stores quantities */
    for(index_tau = 0; index_tau < pba->bt_size; index_tau++)
    {
        dataptr = data + index_tau * number_of_titles;
        pvecback = pba->background_table + index_tau * pba->bg_size;
        storeidx = 0;

        class_store_double(dataptr, 1. / pvecback[pba->index_bg_a] - 1., _TRUE_, storeidx);
        class_store_double(dataptr, pvecback[pba->index_bg_time] / _Gyr_over_Mpc_, _TRUE_, storeidx);
        class_store_double(dataptr, pba->conformal_age - pvecback[pba->index_bg_conf_distance], _TRUE_,
                           storeidx);
        class_store_double(dataptr, pvecback[pba->index_bg_H], _TRUE_, storeidx);
        class_store_double(dataptr, pvecback[pba->index_bg_conf_distance], _TRUE_, storeidx);
        class_store_double(dataptr, pvecback[pba->index_bg_ang_distance], _TRUE_, storeidx);
        class_store_double(dataptr, pvecback[pba->index_bg_lum_distance], _TRUE_, storeidx);
        class_store_double(dataptr, pvecback[pba->index_bg_rs], _TRUE_, storeidx);
        class_store_double(dataptr, pvecback[pba->index_bg_rho_g], _TRUE_, storeidx);
        class_store_double(dataptr, pvecback[pba->index_bg_rho_b], _TRUE_, storeidx);
        class_store_double(dataptr, pvecback[pba->index_bg_rho_cdm], pba->has_cdm, storeidx);
        if(pba->has_ncdm == _TRUE_)
        {
            for(n = 0; n < pba->N_ncdm; n++)
            {
                class_store_double(dataptr, pvecback[pba->index_bg_rho_ncdm1 + n], _TRUE_, storeidx);
                class_store_double(dataptr, pvecback[pba->index_bg_p_ncdm1 + n], _TRUE_, storeidx);
            }
        }
        class_store_double(dataptr, pvecback[pba->index_bg_rho_lambda], pba->has_lambda, storeidx);
        class_store_double(dataptr, pvecback[pba->index_bg_rho_fld], pba->has_fld, storeidx);
        class_store_double(dataptr, pvecback[pba->index_bg_w_fld], pba->has_fld, storeidx);
        class_store_double(dataptr, pvecback[pba->index_bg_rho_ur], pba->has_ur, storeidx);
        class_store_double(dataptr, pvecback[pba->index_bg_rho_idr], pba->has_idr, storeidx);
        class_store_double(dataptr, pvecback[pba->index_bg_rho_idm_dr], pba->has_idm_dr, storeidx);
        class_store_double(dataptr, pvecback[pba->index_bg_rho_crit], _TRUE_, storeidx);
        class_store_double(dataptr, pvecback[pba->index_bg_rho_dcdm], pba->has_dcdm, storeidx);
        class_store_double(dataptr, pvecback[pba->index_bg_rho_dr], pba->has_dr, storeidx);

        /* Non dark energy scalar fields quantities, by Gen Ye */
        if(pba->has_ndescf == _TRUE_)
        {
            for(n = 0; n < pba->N_ndescf; n++)
            {
                class_store_double(dataptr, pvecback[pba->index_bg_rho_ndescf1 + n], _TRUE_, storeidx);
                class_store_double(dataptr, pvecback[pba->index_bg_p_ndescf1 + n], _TRUE_, storeidx);
                class_store_double(dataptr, pvecback[pba->index_bg_phi_ndescf1 + n], _TRUE_, storeidx);
                class_store_double(dataptr, pvecback[pba->index_bg_phi_prime_ndescf1 + n], _TRUE_, storeidx);
            }
        }

        class_store_double(dataptr, pvecback[pba->index_bg_rho_scf], pba->has_scf, storeidx);
        class_store_double(dataptr, pvecback[pba->index_bg_p_scf], pba->has_scf, storeidx);
        class_store_double(dataptr, pvecback[pba->index_bg_p_prime_scf], pba->has_scf, storeidx);
        class_store_double(dataptr, pvecback[pba->index_bg_phi_scf], pba->has_scf, storeidx);
        class_store_double(dataptr, pvecback[pba->index_bg_phi_prime_scf], pba->has_scf, storeidx);
        class_store_double(dataptr, pvecback[pba->index_bg_V_scf], pba->has_scf, storeidx);
        class_store_double(dataptr, pvecback[pba->index_bg_dV_scf], pba->has_scf, storeidx);
        class_store_double(dataptr, pvecback[pba->index_bg_ddV_scf], pba->has_scf, storeidx);

        class_store_double(dataptr, pvecback[pba->index_bg_rho_tot], _TRUE_, storeidx);
        class_store_double(dataptr, pvecback[pba->index_bg_p_tot], _TRUE_, storeidx);
        class_store_double(dataptr, pvecback[pba->index_bg_p_tot_prime], _TRUE_, storeidx);

        class_store_double(dataptr, pvecback[pba->index_bg_D], _TRUE_, storeidx);
        class_store_double(dataptr, pvecback[pba->index_bg_f], _TRUE_, storeidx);
    }

    return _SUCCESS_;
}


/**
 * Subroutine evaluating the derivative with respect to loga
 * of quantities which are integrated (tau, t, etc).
 *
 * This is one of the few functions in the code which is passed to
 * the generic_integrator() routine.  Since generic_integrator()
 * should work with functions passed from various modules, the format
 * of the arguments is a bit special:
 *
 * - fixed input parameters and workspaces are passed through a generic
 * pointer. Here, this is just a pointer to the background structure
 * and to a background vector, but generic_integrator() doesn't know
 * its fine structure.
 *
 * - the error management is a bit special: errors are not written as
 * usual to pba->error_message, but to a generic error_message passed
 * in the list of arguments.
 *
 * @param loga                     Input: current value of log(a)
 * @param y                        Input: vector of variable
 * @param dy                       Output: its derivative (already allocated)
 * @param parameters_and_workspace Input: pointer to fixed parameters (e.g. indices)
 * @param error_message            Output: error message
 */

int background_derivs(
    double loga,
    double* y, /* vector with argument y[index_bi] (must be already allocated with size pba->bi_size) */
    double* dy, /* vector with argument dy[index_bi]
                                     (must be already allocated with
                                     size pba->bi_size) */
    void * parameters_and_workspace,
    ErrorMsg error_message
)
{

    /** Summary: */

    /** - define local variables */

    struct background_parameters_and_workspace * pbpaw;
    struct background * pba;
    double * pvecback, a, H, rho_M;
    int n_ndescf; // for ndescf components, by Gen Ye

    pbpaw = parameters_and_workspace;
    pba =  pbpaw->pba;
    pvecback = pbpaw->pvecback;

    /** - scale factor a (in fact, given our normalisation conventions, this stands for a/a_0) */
    a = exp(loga);

    //if(pba->adaptive_step_bk_sampling == _TRUE_) dy[pba->index_bi_loga]=1.;

    /** - calculate functions of \f$ a \f$ with background_functions() */
    class_call(background_functions(pba, a, y, normal_info, pvecback),
               pba->error_message,
               error_message);

    /** - Short hand notation for Hubble */
    H = pvecback[pba->index_bg_H];

    /** - calculate derivative of cosmological time \f$ dt/dloga = 1/H \f$ */
    dy[pba->index_bi_time] = 1. / H;

    /** - calculate derivative of conformal time \f$ d\tau/dloga = 1/aH \f$ */
    dy[pba->index_bi_tau] = 1. / a / H;

    class_test(pvecback[pba->index_bg_rho_g] <= 0.,
               error_message,
               "rho_g = %e instead of strictly positive", pvecback[pba->index_bg_rho_g]);

    /** - calculate detivative of sound horizon \f$ drs/dloga = drs/dtau * dtau/dloga = c_s/aH \f$*/
    dy[pba->index_bi_rs] = 1. / a / H / sqrt(3.*(1. + 3.*pvecback[pba->index_bg_rho_b] / 4. /
                           pvecback[pba->index_bg_rho_g])) * sqrt(1. - pba->K * y[pba->index_bi_rs] *
                                   y[pba->index_bi_rs]); // TBC: curvature correction

    /** - solve second order growth equation \f$ [D''(\tau)=-aHD'(\tau)+3/2 a^2 \rho_M D(\tau) \f$
        written as \f$ dD/dloga = D' / (aH) \f$ and \f$ dD'/dloga = -D' + (3/2) (a/H) \rho_M D \f$ */
    rho_M = pvecback[pba->index_bg_rho_b];
    if(pba->has_cdm == _TRUE_)
    {
        rho_M += pvecback[pba->index_bg_rho_cdm];
    }
    if(pba->has_idm_dr == _TRUE_)
    {
        rho_M += pvecback[pba->index_bg_rho_idm_dr];
    }

    dy[pba->index_bi_D] = y[pba->index_bi_D_prime] / a / H;
    dy[pba->index_bi_D_prime] = -y[pba->index_bi_D_prime] + 1.5 * a * rho_M * y[pba->index_bi_D] / H;

    if(pba->has_dcdm == _TRUE_)
    {
        /** - compute dcdm density \f$ d\rho/dloga = -3 \rho - \Gamma/H \rho \f$*/
        dy[pba->index_bi_rho_dcdm] = -3.*y[pba->index_bi_rho_dcdm] - pba->Gamma_dcdm / H *
                                     y[pba->index_bi_rho_dcdm];
    }

    if((pba->has_dcdm == _TRUE_) && (pba->has_dr == _TRUE_))
    {
        /** - Compute dr density \f$ d\rho/dloga = -4\rho - \Gamma/H \rho \f$ */
        dy[pba->index_bi_rho_dr] = -4.*y[pba->index_bi_rho_dr] + pba->Gamma_dcdm / H *
                                   y[pba->index_bi_rho_dcdm];
    }

    if(pba->has_fld == _TRUE_)
    {
        /** - Compute fld density \f$ d\rho/dloga = -3 (1+w_{fld}(a)) \rho \f$ */
        dy[pba->index_bi_rho_fld] = -3.*(1. + pvecback[pba->index_bg_w_fld]) * y[pba->index_bi_rho_fld];
    }

    if(pba->has_ndescf == _TRUE_)
    {
        /** - Non dark energy scalar field equation, by Gen Ye */
        for(n_ndescf = 0; n_ndescf < pba->N_ndescf; n_ndescf++)
        {
            if(pba->ndescf_use_kg[n_ndescf] == _TRUE_)
            {
                //printf("------->use KG\n");
                dy[pba->index_bi_phi_ndescf1 + n_ndescf] = y[pba->index_bi_phi_prime_ndescf1 + n_ndescf] / a / H;
                dy[pba->index_bi_phi_prime_ndescf1 + n_ndescf] = - 2 * y[pba->index_bi_phi_prime_ndescf1 + n_ndescf] - a *
                        dV_ndescf(pba, y[pba->index_bi_phi_ndescf1 + n_ndescf], n_ndescf) / H;
                // rho_ndescf is handled by background_functions() in this case
                dy[pba->index_bi_rho_ndescf1 + n_ndescf] = 0.;
            }
            else // fluid approximation for oscillating fields
            {
                //printf("------->use fld\n");
                /* fluid equation \f$ d\rho/dlna = -3(\rho+p) \f$ */
                dy[pba->index_bi_rho_ndescf1 + n_ndescf] = -3.*(1. + pba->ndescf_w_fld[n_ndescf]) * y[pba->index_bi_rho_ndescf1 +
                                                    n_ndescf];
                // we loose track of phi, phi' when using fld approx
                dy[pba->index_bi_phi_ndescf1 + n_ndescf] = 0.;
                dy[pba->index_bi_phi_prime_ndescf1 + n_ndescf] = 0.;
            }
        }
    }

    if(pba->has_scf == _TRUE_)
    {
        /** - Scalar field equation: \f$ \phi'' + 2 a H \phi' + a^2 dV = 0 \f$  (note H is wrt cosmological time)
            written as \f$ d\phi/dlna = phi' / (aH) \f$ and \f$ d\phi'/dlna = -2*phi' - (a/H) dV \f$ */
        dy[pba->index_bi_phi_scf] = y[pba->index_bi_phi_prime_scf] / a / H;
        dy[pba->index_bi_phi_prime_scf] = - 2 * y[pba->index_bi_phi_prime_scf] - a * dV_scf(pba,
                                          y[pba->index_bi_phi_scf]) / H ;
    }

    return _SUCCESS_;

}

/**
 * At some step during the integraton of the background equations,
 * this function extracts the qantities that we want to keep memory
 * of, and stores them in a row of the background table (as well as
 * extra tables: z_table, tau_table).
 *
 * This is one of the few functions in the code which is passed to the generic_integrator() routine.
 * Since generic_integrator() should work with functions passed from various modules, the format of the arguments
 * is a bit special:
 * - fixed parameters and workspaces are passed through a generic pointer.
 *   generic_integrator() doesn't know the content of this pointer.
 * - the error management is a bit special: errors are not written as usual to pba->error_message, but to a generic
 *   error_message passed in the list of arguments.
 *
 * @param loga                     Input: current value of log(a)
 * @param y                        Input: current vector of integrated quantities (with index_bi)
 * @param dy                       Input: current derivative of y w.r.t log(a)
 * @param index_loga               Input: index of the log(a) value within the background_table
 * @param parameters_and_workspace Input/output: fixed parameters (e.g. indices), workspace, background structure where the output is written...
 * @param error_message            Output: error message
 */

int background_sources(
    double loga,
    double * y,
    double * dy,
    int index_loga,
    void * parameters_and_workspace,
    ErrorMsg error_message
)
{

    struct background_parameters_and_workspace * pbpaw;
    struct background * pba;
    double a;
    double * bg_table_row;

    pbpaw = parameters_and_workspace;
    pba =  pbpaw->pba;

    /** - localize the row inside background_table where the current values must be stored */
    bg_table_row = pba->background_table + index_loga * pba->bg_size;

    /** - scale factor a (in fact, given our normalisation conventions, this stands for a/a_0) */
    a = exp(loga);

    /** - corresponding redhsift 1/a-1 */
    pba->z_table[index_loga] = MAX(0., 1. / a - 1.);

    /** - corresponding conformal time */
    pba->tau_table[index_loga] = y[pba->index_bi_tau];

    /** -> compute all other quantities depending only on a + {B} variables and get them stored
       in one row of background_table
       The value of {B} variables in pData are also copied to pvecback.*/
    class_call(background_functions(pba, a, y, long_info, bg_table_row),
               pba->error_message,
               pba->error_message);

    return _SUCCESS_;

}

/**
 * Evalute the typical timescale for the integration of he background
 * over loga=log(a/a_0). This is only required for rkck, but not for
 * the ndf15 evolver.
 *
 * The evolver will take steps equal to this value times
 * ppr->background_integration_stepsize.  Since our variable of
 * integration is loga, and the time steps are (delta a)/a, the
 * reference timescale is precisely one, i.e., the code will take some
 * steps such that (delta a)/a = ppr->background_integration_stepsize.
 *
 * The argument list is predetermined by the format of
 * generic_evolver; however in this particular case, they are never
 * used.
 *
 * This is one of the few functions in the code which is passed to the generic_integrator() routine.
 * Since generic_integrator() should work with functions passed from various modules, the format of the arguments
 * is a bit special:
 * - fixed parameters and workspaces are passed through a generic pointer (void *).
 *   generic_integrator() doesn't know the content of this pointer.
 * - the error management is a bit special: errors are not written as usual to pba->error_message, but to a generic
 *   error_message passed in the list of arguments.
 *
 * @param loga                     Input: current value of log(a/a_0)
 * @param parameters_and_workspace Input: fixed parameters (e.g. indices), workspace, approximation used, etc.
 * @param timescale                Output: perturbation variation timescale
 * @param error_message            Output: error message
 */

int background_timescale(
    double loga,
    void * parameters_and_workspace,
    double * timescale,
    ErrorMsg error_message
)
{
    /** Summary: */

    /** - define local variables */
    /* block added by Gen Ye */
    struct background_parameters_and_workspace * pbpaw;
    struct background * pba;
    double * pvecback;

    pbpaw = parameters_and_workspace;
    pba = pbpaw->pba;
    pvecback = pbpaw->pvecback;

    double scale_ndescf;
    int n_ndescf, sorted_n_ndescf;

    *timescale = 1.;

    if(pba->has_ndescf == _TRUE_)
    {
        for(n_ndescf = 0; n_ndescf < pba->N_ndescf; n_ndescf++)
        {
            sorted_n_ndescf = pba->ndescf_sorted_index[n_ndescf];
            if(loga < pba->ndescf_loga_fld[sorted_n_ndescf])
            {
                *timescale = MIN(pvecback[pba->index_bg_H] / pba->ndescf_m_fld[sorted_n_ndescf], *timescale);
                break;
            }
        }
    }

    return _SUCCESS_;
}

/**
 * Function outputting the fractions Omega of the total critical density
 * today, and also the reduced fractions omega=Omega*h*h
 *
 * It also prints the total budgets of non-relativistic, relativistic,
 * and other contents, and of the total
 *
 * @param pba                      Input: Pointer to background structure
 * @return the error status
 */

int background_output_budget(
    struct background* pba
)
{

    double budget_matter, budget_radiation, budget_other, budget_neutrino;
    int index_ncdm, index_ndescf; // add ndescf index, by Gen Ye

    budget_matter = 0;
    budget_radiation = 0;
    budget_other = 0;
    budget_neutrino = 0;

    //The name for the class_print_species macro can be at most 30 characters total
    if(pba->background_verbose > 1)
    {

        printf(" ---------------------------- Budget equation ----------------------- \n");

        printf(" ---> Nonrelativistic Species \n");
        class_print_species("Bayrons", b);
        budget_matter += pba->Omega0_b;
        if(pba->has_cdm == _TRUE_)
        {
            class_print_species("Cold Dark Matter", cdm);
            budget_matter += pba->Omega0_cdm;
        }
        if(pba->has_idm_dr == _TRUE_)
        {
            class_print_species("Interacting Dark Matter - DR ", idm_dr);
            budget_matter += pba->Omega0_idm_dr;
        }
        if(pba->has_dcdm == _TRUE_)
        {
            class_print_species("Decaying Cold Dark Matter", dcdm);
            budget_matter += pba->Omega0_dcdm;
        }


        printf(" ---> Relativistic Species \n");
        class_print_species("Photons", g);
        budget_radiation += pba->Omega0_g;
        if(pba->has_ur == _TRUE_)
        {
            class_print_species("Ultra-relativistic relics", ur);
            budget_radiation += pba->Omega0_ur;
        }
        if(pba->has_dr == _TRUE_)
        {
            class_print_species("Dark Radiation (from decay)", dr);
            budget_radiation += pba->Omega0_dr;
        }
        if(pba->has_idr == _TRUE_)
        {
            class_print_species("Interacting Dark Radiation", idr);
            budget_radiation += pba->Omega0_idr;
        }

        if(pba->N_ncdm > 0)
        {
            printf(" ---> Massive Neutrino Species \n");
        }
        if(pba->N_ncdm > 0)
        {
            for(index_ncdm = 0; index_ncdm < pba->N_ncdm; ++index_ncdm)
            {
                printf("-> %-26s%-4d Omega = %-15g , omega = %-15g\n", "Neutrino Species Nr.", index_ncdm + 1,
                       pba->Omega0_ncdm[index_ncdm], pba->Omega0_ncdm[index_ncdm]*pba->h * pba->h);
                budget_neutrino += pba->Omega0_ncdm[index_ncdm];
            }
        }

        if((pba->has_lambda == _TRUE_) || (pba->has_fld == _TRUE_) || (pba->has_scf == _TRUE_)
                || (pba->has_ndescf == _TRUE_)
                || (pba->has_curvature == _TRUE_))
        {
            printf(" ---> Other Content \n");
        }
        if(pba->has_lambda == _TRUE_)
        {
            class_print_species("Cosmological Constant", lambda);
            budget_other += pba->Omega0_lambda;
        }
        if(pba->has_fld == _TRUE_)
        {
            class_print_species("Dark Energy Fluid", fld);
            budget_other += pba->Omega0_fld;
        }
        if(pba->has_scf == _TRUE_)
        {
            class_print_species("Scalar Field", scf);
            budget_other += pba->Omega0_scf;
        }
        if(pba->has_ndescf == _TRUE_) // ndescf case by Gen Ye
        {
            class_print_species("Additional scalar fields", ndescf_tot_real)
            budget_other += pba->Omega0_ndescf_tot_real;
        }
        if(pba->has_curvature == _TRUE_)
        {
            class_print_species("Spatial Curvature", k);
            budget_other += pba->Omega0_k;
        }

        printf(" ---> Total budgets \n");
        printf(" Radiation                        Omega = %-15g , omega = %-15g \n", budget_radiation,
               budget_radiation * pba->h * pba->h);
        printf(" Non-relativistic                 Omega = %-15g , omega = %-15g \n", budget_matter,
               budget_matter * pba->h * pba->h);
        if(pba->N_ncdm > 0)
        {
            printf(" Neutrinos                        Omega = %-15g , omega = %-15g \n", budget_neutrino,
                   budget_neutrino * pba->h * pba->h);
        }
        if((pba->has_lambda == _TRUE_) || (pba->has_fld == _TRUE_) || (pba->has_scf == _TRUE_)
                || (pba->has_ndescf == _TRUE_)
                || (pba->has_curvature == _TRUE_))
        {
            printf(" Other Content                    Omega = %-15g , omega = %-15g \n", budget_other,
                   budget_other * pba->h * pba->h);
        }
        printf(" TOTAL                            Omega = %-15g , omega = %-15g \n",
               budget_radiation + budget_matter + budget_neutrino + budget_other,
               (budget_radiation + budget_matter + budget_neutrino + budget_other)*pba->h * pba->h);

        printf(" -------------------------------------------------------------------- \n");
    }

    return _SUCCESS_;
}

/**
 * Scalar field potential and its derivatives with respect to the field _scf
 * For Albrecht & Skordis model: 9908085
 * - \f$ V = V_{p_{scf}}*V_{e_{scf}} \f$
 * - \f$ V_e =  \exp(-\lambda \phi) \f$ (exponential)
 * - \f$ V_p = (\phi - B)^\alpha + A \f$ (polynomial bump)
 *
 * TODO:
 * - Add some functionality to include different models/potentials (tuning would be difficult, though)
 * - Generalize to Kessence/Horndeski/PPF and/or couplings
 * - A default module to numerically compute the derivatives when no analytic functions are given should be added.
 * - Numerical derivatives may further serve as a consistency check.
 *
 */

/**
 *
 * The units of phi, tau in the derivatives and the potential V are the following:
 * - phi is given in units of the reduced Planck mass \f$ m_{pl} = (8 \pi G)^{(-1/2)}\f$
 * - tau in the derivative is given in units of Mpc.
 * - the potential \f$ V(\phi) \f$ is given in units of \f$ m_{pl}^2/Mpc^2 \f$.
 * With this convention, we have
 * \f$ \rho^{class} = (8 \pi G)/3 \rho^{physical} = 1/(3 m_{pl}^2) \rho^{physical} = 1/3 * [ 1/(2a^2) (\phi')^2 + V(\phi) ] \f$
 and \f$ \rho^{class} \f$ has the proper dimension \f$ Mpc^-2 \f$.
*/

double V_e_scf(struct background *pba,
               double phi
              )
{
    double scf_lambda = pba->scf_parameters[0];
    //  double scf_alpha  = pba->scf_parameters[1];
    //  double scf_A      = pba->scf_parameters[2];
    //  double scf_B      = pba->scf_parameters[3];

    return  exp(-scf_lambda * phi);
}

double dV_e_scf(struct background *pba,
                double phi
               )
{
    double scf_lambda = pba->scf_parameters[0];
    //  double scf_alpha  = pba->scf_parameters[1];
    //  double scf_A      = pba->scf_parameters[2];
    //  double scf_B      = pba->scf_parameters[3];

    return -scf_lambda * V_scf(pba, phi);
}

double ddV_e_scf(struct background *pba,
                 double phi
                )
{
    double scf_lambda = pba->scf_parameters[0];
    //  double scf_alpha  = pba->scf_parameters[1];
    //  double scf_A      = pba->scf_parameters[2];
    //  double scf_B      = pba->scf_parameters[3];

    return pow(-scf_lambda, 2) * V_scf(pba, phi);
}


/** parameters and functions for the polynomial coefficient
 * \f$ V_p = (\phi - B)^\alpha + A \f$(polynomial bump)
 *
 * double scf_alpha = 2;
 *
 * double scf_B = 34.8;
 *
 * double scf_A = 0.01; (values for their Figure 2)
 */

double V_p_scf(
    struct background *pba,
    double phi)
{
    //  double scf_lambda = pba->scf_parameters[0];
    double scf_alpha  = pba->scf_parameters[1];
    double scf_A      = pba->scf_parameters[2];
    double scf_B      = pba->scf_parameters[3];

    return  pow(phi - scf_B,  scf_alpha) +  scf_A;
}

double dV_p_scf(
    struct background *pba,
    double phi)
{

    //  double scf_lambda = pba->scf_parameters[0];
    double scf_alpha  = pba->scf_parameters[1];
    //  double scf_A      = pba->scf_parameters[2];
    double scf_B      = pba->scf_parameters[3];

    return   scf_alpha * pow(phi -  scf_B,  scf_alpha - 1);
}

double ddV_p_scf(
    struct background *pba,
    double phi)
{
    //  double scf_lambda = pba->scf_parameters[0];
    double scf_alpha  = pba->scf_parameters[1];
    //  double scf_A      = pba->scf_parameters[2];
    double scf_B      = pba->scf_parameters[3];

    return  scf_alpha * (scf_alpha - 1.) * pow(phi -  scf_B,  scf_alpha - 2);
}

/* other potential functions, by Gen Ye */
/** phi2n \f$ V = V_0*(\phi/M_p)^2n-V_ads+V_lambda \f$ */
double V_phi2n_scf(
    struct background *pba,
    double phi)
{
    double V0 = pba->scf_parameters[0];
    double n2 = pba->scf_parameters[1] * 2;
    double Vads = pba->scf_parameters[2];
    return V0 * pow(phi, n2) - Vads;
}

double dV_phi2n_scf(
    struct background *pba,
    double phi)
{
    double V0 = pba->scf_parameters[0];
    double n2 = pba->scf_parameters[1] * 2;
    return n2 * V0 * pow(phi, n2 - 1);
}

double ddV_phi2n_scf(
    struct background *pba,
    double phi)
{
    double V0 = pba->scf_parameters[0];
    double n2 = pba->scf_parameters[1] * 2;
    if(n2 == 2) return 2 * V0;
    else return n2 * (n2 - 1) * V0 * pow(phi, n2 - 2);
}

/** alphaAdS */
double V_alphaads_scf(
    struct background *pba,
    double phi)
{
    double scf_V0 = pba->scf_parameters[0];
    double scf_gamma  = pba->scf_parameters[1];
    double scf_6alpha_root  = pba->scf_parameters[2];
    double scf_Vlambda = pba->scf_parameters[3];
    return scf_V0 * pow(expm1(-scf_gamma * tanh(phi / scf_6alpha_root)),
                        2) - scf_V0 + scf_Vlambda;
}

double dV_alphaads_scf(
    struct background *pba,
    double phi)
{
    double scf_V0 = pba->scf_parameters[0];
    double scf_gamma  = pba->scf_parameters[1];
    double scf_6alpha_root  = pba->scf_parameters[2];
    return (-2 * expm1(-(scf_gamma * tanh(phi / scf_6alpha_root))) * scf_gamma *
            scf_V0 * pow(cosh(phi / scf_6alpha_root),
                         -2)) / (exp(scf_gamma * tanh(phi / scf_6alpha_root)) * scf_6alpha_root);
}

double ddV_alphaads_scf(
    struct background *pba,
    double phi)
{
    double scf_V0 = pba->scf_parameters[0];
    double scf_gamma  = pba->scf_parameters[1];
    double scf_6alpha_root  = pba->scf_parameters[2];
    return (-2 * scf_gamma * scf_V0 * pow(cosh(phi / scf_6alpha_root),
                                          -4) * ((-1 + expm1(scf_gamma * tanh(phi / scf_6alpha_root))) * scf_gamma +
                                                  expm1(
                                                          scf_gamma * tanh(phi / scf_6alpha_root)) * sinh((2 * phi) / scf_6alpha_root))) /
           (exp(2 * scf_gamma * tanh(phi / scf_6alpha_root)) * pow(scf_6alpha_root, 2));
}

/** connection function contiuous up to second derivative
* RHS is constant zero
*/
double V_cnt_scf(
    struct background * pba,
    double phi,
    double x1,
    double x2,
    double V,
    double dV,
    double ddV)
{
    double a0 = -(20 * V * pow(x1, 2) * pow(x2, 3) - 8 * dV * pow(x1, 3) * pow(x2, 3) +
                  ddV * pow(x1, 4) * pow(x2, 3) - 10 * V * x1 * pow(x2, 4) + 10 * dV * pow(x1, 2) * pow(x2, 4) -
                  2 * ddV * pow(x1, 3) * pow(x2, 4) + 2 * V * pow(x2, 5) - 2 * dV * x1 * pow(x2, 5) +
                  ddV * pow(x1, 2) * pow(x2, 5)) / (2.*pow(x1 - x2, 5));
    double a1 = -(-60 * V * pow(x1, 2) * pow(x2, 2) + 24 * dV * pow(x1, 3) * pow(x2, 2) -
                  3 * ddV * pow(x1, 4) * pow(x2, 2) - 16 * dV * pow(x1, 2) * pow(x2, 3) +
                  4 * ddV * pow(x1, 3) * pow(x2, 3) - 10 * dV * x1 * pow(x2, 4) + ddV * pow(x1, 2) * pow(x2, 4) +
                  2 * dV * pow(x2, 5) - 2 * ddV * x1 * pow(x2, 5)) / (2.*pow(x1 - x2, 5));
    double a2 = -(x2 * (60 * V * pow(x1, 2) - 24 * dV * pow(x1, 3) + 3 * ddV * pow(x1,
                        4) + 60 * V * x1 * x2 -
                        12 * dV * pow(x1, 2) * x2 + 36 * dV * x1 * pow(x2, 2) - 8 * ddV * pow(x1, 2) * pow(x2, 2) +
                        4 * ddV * x1 * pow(x2, 3) + ddV * pow(x2, 4))) / (2.*pow(x1 - x2, 5));
    double a3 = -(-20 * V * pow(x1, 2) + 8 * dV * pow(x1, 3) - ddV * pow(x1, 4) - 80 * V * x1 * x2 +
                  32 * dV * pow(x1, 2) * x2 - 4 * ddV * pow(x1, 3) * x2 - 20 * V * pow(x2, 2) -
                  28 * dV * x1 * pow(x2, 2) + 8 * ddV * pow(x1, 2) * pow(x2, 2) - 12 * dV * pow(x2, 3) -
                  3 * ddV * pow(x2, 4)) / (2.*pow(x1 - x2, 5));
    double a4 = -(30 * V * x1 - 14 * dV * pow(x1, 2) + 2 * ddV * pow(x1,
                  3) + 30 * V * x2 - 2 * dV * x1 * x2 -
                  ddV * pow(x1, 2) * x2 + 16 * dV * pow(x2, 2) - 4 * ddV * x1 * pow(x2, 2) + 3 * ddV * pow(x2, 3)) /
                (2.*pow(x1 - x2, 5));
    double a5 = -(-12 * V + 6 * dV * x1 - ddV * pow(x1,
                  2) - 6 * dV * x2 + 2 * ddV * x1 * x2 - ddV * pow(x2, 2)) /
                (2.*pow(x1 - x2, 5));
    return a0 + a1 * phi + a2 * pow(phi, 2) + a3 * pow(phi, 3) + a4 * pow(phi, 4) + a5 * pow(phi, 5);
}

double dV_cnt_scf(
    struct background * pba,
    double phi,
    double x1,
    double x2,
    double V,
    double dV,
    double ddV)
{
    double a1 = -(-60 * V * pow(x1, 2) * pow(x2, 2) + 24 * dV * pow(x1, 3) * pow(x2, 2) -
                  3 * ddV * pow(x1, 4) * pow(x2, 2) - 16 * dV * pow(x1, 2) * pow(x2, 3) +
                  4 * ddV * pow(x1, 3) * pow(x2, 3) - 10 * dV * x1 * pow(x2, 4) + ddV * pow(x1, 2) * pow(x2, 4) +
                  2 * dV * pow(x2, 5) - 2 * ddV * x1 * pow(x2, 5)) / (2.*pow(x1 - x2, 5));
    double a2 = -(x2 * (60 * V * pow(x1, 2) - 24 * dV * pow(x1, 3) + 3 * ddV * pow(x1,
                        4) + 60 * V * x1 * x2 -
                        12 * dV * pow(x1, 2) * x2 + 36 * dV * x1 * pow(x2, 2) - 8 * ddV * pow(x1, 2) * pow(x2, 2) +
                        4 * ddV * x1 * pow(x2, 3) + ddV * pow(x2, 4))) / (2.*pow(x1 - x2, 5));
    double a3 = -(-20 * V * pow(x1, 2) + 8 * dV * pow(x1, 3) - ddV * pow(x1, 4) - 80 * V * x1 * x2 +
                  32 * dV * pow(x1, 2) * x2 - 4 * ddV * pow(x1, 3) * x2 - 20 * V * pow(x2, 2) -
                  28 * dV * x1 * pow(x2, 2) + 8 * ddV * pow(x1, 2) * pow(x2, 2) - 12 * dV * pow(x2, 3) -
                  3 * ddV * pow(x2, 4)) / (2.*pow(x1 - x2, 5));
    double a4 = -(30 * V * x1 - 14 * dV * pow(x1, 2) + 2 * ddV * pow(x1,
                  3) + 30 * V * x2 - 2 * dV * x1 * x2 -
                  ddV * pow(x1, 2) * x2 + 16 * dV * pow(x2, 2) - 4 * ddV * x1 * pow(x2, 2) + 3 * ddV * pow(x2, 3)) /
                (2.*pow(x1 - x2, 5));
    double a5 = -(-12 * V + 6 * dV * x1 - ddV * pow(x1,
                  2) - 6 * dV * x2 + 2 * ddV * x1 * x2 - ddV * pow(x2, 2)) /
                (2.*pow(x1 - x2, 5));
    return a1 + 2 * a2 * phi + 3 * a3 * pow(phi, 2) + 4 * a4 * pow(phi, 3) + 5 * a5 * pow(phi, 4);
}

double ddV_cnt_scf(
    struct background * pba,
    double phi,
    double x1,
    double x2,
    double V,
    double dV,
    double ddV)
{
    double a2 = -(x2 * (60 * V * pow(x1, 2) - 24 * dV * pow(x1, 3) + 3 * ddV * pow(x1,
                        4) + 60 * V * x1 * x2 -
                        12 * dV * pow(x1, 2) * x2 + 36 * dV * x1 * pow(x2, 2) - 8 * ddV * pow(x1, 2) * pow(x2, 2) +
                        4 * ddV * x1 * pow(x2, 3) + ddV * pow(x2, 4))) / (2.*pow(x1 - x2, 5));
    double a3 = -(-20 * V * pow(x1, 2) + 8 * dV * pow(x1, 3) - ddV * pow(x1, 4) - 80 * V * x1 * x2 +
                  32 * dV * pow(x1, 2) * x2 - 4 * ddV * pow(x1, 3) * x2 - 20 * V * pow(x2, 2) -
                  28 * dV * x1 * pow(x2, 2) + 8 * ddV * pow(x1, 2) * pow(x2, 2) - 12 * dV * pow(x2, 3) -
                  3 * ddV * pow(x2, 4)) / (2.*pow(x1 - x2, 5));
    double a4 = -(30 * V * x1 - 14 * dV * pow(x1, 2) + 2 * ddV * pow(x1,
                  3) + 30 * V * x2 - 2 * dV * x1 * x2 -
                  ddV * pow(x1, 2) * x2 + 16 * dV * pow(x2, 2) - 4 * ddV * x1 * pow(x2, 2) + 3 * ddV * pow(x2, 3)) /
                (2.*pow(x1 - x2, 5));
    double a5 = -(-12 * V + 6 * dV * x1 - ddV * pow(x1,
                  2) - 6 * dV * x2 + 2 * ddV * x1 * x2 - ddV * pow(x2, 2)) /
                (2.*pow(x1 - x2, 5));
    return 2 * a2 + 6 * a3 * phi + 12 * a4 * pow(phi, 2) + 20 * a5 * pow(phi, 3);
}

/** Fianlly we can obtain the overall potential \f$ V = V_p*V_e \f$
 */

double V_scf(
    struct background * pba,
    double phi)
{
    int psize = pba->scf_parameters_size; // convenience variable
    double ptype = pba->scf_parameters[psize - 3]; // last entry flags potential type

    /** loop over potential type, currently supported are:
     *       ptype<=0   -> original implementation in CLASS
     *       ptype=0.5 -> phi2n
     *       ptype=1.5 -> adsede-alphaAdS
     *       ptype=2.5 -> adsede-phi2n */
    class_test(ptype > 3,
               pba->error_message,
               "Invalid potential flag ptype=%g, only 0.5, 1.5, 2.5 is currently supported. Check your parameters_smg in the .ini file.\n",
               ptype);
    if(ptype <= 0) return V_e_scf(pba, phi) * V_p_scf(pba, phi);
    else if(ptype <= 1)
    {
        double Vlambda = pba->scf_parameters[3];
        return V_phi2n_scf(pba, phi) + Vlambda;
    }
    else if(ptype <= 2) return V_alphaads_scf(pba, phi);
    else if(ptype <= 3)
    {
        double V0 = pba->scf_parameters[0];
        double n2 = pba->scf_parameters[1] * 2.;
        double Vads = pba->scf_parameters[2];
        double Vlambda = pba->scf_parameters[3];
        double x1 = 0.9 * pow(Vads / V0, 1. / n2);
        double x2 = 1.22 * x1;
        double Vi = V_phi2n_scf(pba, x1);
        double dVi = dV_phi2n_scf(pba, x1);
        double ddVi = ddV_phi2n_scf(pba, x1);
        if(phi < x1) return V_phi2n_scf(pba, phi) + Vlambda;
        else if(phi < x2) return V_cnt_scf(pba, phi, x1, x2, Vi, dVi, ddVi) + Vlambda;
        else return Vlambda;
    }
}

double dV_scf(
    struct background * pba,
    double phi)
{
    int psize = pba->scf_parameters_size;
    double ptype = pba->scf_parameters[psize - 3];
    class_test(ptype > 3,
               pba->error_message,
               "Invalid potential flag ptype=%g, only 0.5, 1.5, 2.5 is currently supported. Check your parameters_smg in the .ini file.\n",
               ptype);
    if(ptype <= 0)
    {
        return dV_e_scf(pba, phi) * V_p_scf(pba, phi) + V_e_scf(pba, phi) * dV_p_scf(pba,
                phi);
    }
    else if(ptype <= 1) return dV_phi2n_scf(pba, phi);
    else if(ptype <= 2) return dV_alphaads_scf(pba, phi);
    else if(ptype <= 3)
    {
        double V0 = pba->scf_parameters[0];
        double n2 = pba->scf_parameters[1] * 2.;
        double Vads = pba->scf_parameters[2];
        double x1 = 0.9 * pow(Vads / V0, 1. / n2);
        double x2 = 1.22 * x1;
        double Vi = V_phi2n_scf(pba, x1);
        double dVi = dV_phi2n_scf(pba, x1);
        double ddVi = ddV_phi2n_scf(pba, x1);
        if(phi < x1) return dV_phi2n_scf(pba, phi);
        else if(phi < x2) return dV_cnt_scf(pba, phi, x1, x2, Vi, dVi, ddVi);
        else return 0.;
    }
}

double ddV_scf(
    struct background * pba,
    double phi)
{
    int psize = pba->scf_parameters_size;
    double ptype = pba->scf_parameters[psize - 3];
    class_test(ptype > 3,
               pba->error_message,
               "Invalid potential flag ptype=%g, only 0.5, 1.5, 2.5 is currently supported. Check your parameters_smg in the .ini file.\n",
               ptype);
    if(ptype <= 0)
    {
        return ddV_e_scf(pba, phi) * V_p_scf(pba, phi) + 2 * dV_e_scf(pba, phi) * dV_p_scf(pba,
                phi) + V_e_scf(pba, phi) * ddV_p_scf(pba, phi);
    }
    else if(ptype <= 1) return ddV_phi2n_scf(pba, phi);
    else if(ptype <= 2) return ddV_alphaads_scf(pba, phi);
    else if(ptype <= 3)
    {
        double V0 = pba->scf_parameters[0];
        double n2 = pba->scf_parameters[1] * 2.;
        double Vads = pba->scf_parameters[2];
        double x1 = 0.9 * pow(Vads / V0, 1. / n2);
        double x2 = 1.22 * x1;
        double Vi = V_phi2n_scf(pba, x1);
        double dVi = dV_phi2n_scf(pba, x1);
        double ddVi = ddV_phi2n_scf(pba, x1);
        if(phi < x1) return ddV_phi2n_scf(pba, phi);
        else if(phi < x2) return ddV_cnt_scf(pba, phi, x1, x2, Vi, dVi, ddVi);
        else return 0.;
    }
}

/* Put all non dark energy scalar field potentials here, use n_ndescf for labelling */
double V_ndescf(
    struct background * pba,
    double phi,
    int n_ndescf)
{
    int i = n_ndescf * pba->ndescf_parameters_size;
    // The last third parameter flags potential type
    double ptype = pba->ndescf_parameters[i + pba->ndescf_parameters_size - 3];
    double rlt;
    if(ptype <= 0)
    {
        double fa = pba->ndescf_parameters[i];
        double ma = pba->ndescf_parameters[i + 1];
        rlt = pow(ma, 2) * pow(fa, 2) * (1 - cos(phi / fa));
    }
    return rlt;

}

double dV_ndescf(
    struct background * pba,
    double phi,
    int n_ndescf)
{
    int i = n_ndescf * pba->ndescf_parameters_size;
    double ptype = pba->ndescf_parameters[i + pba->ndescf_parameters_size - 3];
    double rlt;
    if(ptype <= 0)
    {
        double fa = pba->ndescf_parameters[i];
        double ma = pba->ndescf_parameters[i + 1];
        rlt = pow(ma, 2) * fa * sin(phi / fa);
    }
    return rlt;
}

double ddV_ndescf(
    struct background * pba,
    double phi,
    int n_ndescf)
{
    int i = n_ndescf * pba->ndescf_parameters_size;
    double ptype = pba->ndescf_parameters[i + pba->ndescf_parameters_size - 3];
    double rlt;
    if(ptype <= 0)
    {
        double fa = pba->ndescf_parameters[i];
        double ma = pba->ndescf_parameters[i + 1];
        rlt = pow(ma, 2) * cos(phi / fa);
    }
    return rlt;
}
/* end of modification by Gen Ye */

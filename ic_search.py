# -*- coding: utf-8 -*-
"""
@author: Gen Ye
"""
from __future__ import print_function
from math import exp, expm1, log1p, tanh, atanh, sqrt, cosh, sinh
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import root_scalar, minimize_scalar, root


# the MCMC parameters are scf_f and ln(1+z_c)
# this code shoots for scalar field initial condition (phi_i,V_0) which
# reproduces scf_f at z_c
# not very precise but good enough

###################
# add support for supergravity potential
# V=V0*(exp(gamma*tanh(phi/root6alpha))-1)^2-(V0-V_Lambda)
###################

class scf_ic_finder:

    # input is the desired precision
    def __init__(
            self,
            precision=0.01,
            initial_time=-4.6,
            verbose=0):
        # model parameters
        self.smg_params = np.zeros(3)
        self.scf_params = np.zeros(4)
        # precision and initial time
        self.precision = precision
        self.t_i = initial_time
        # whether print results
        self.verbose = verbose
        # default cosmo parameters
        self.omega_m = 0.142
        self.T_T0 = 1
        # default MCMC parameters
        self.scf_f = 0.1
        self.z_c = expm1(8.2)
        # default scf parameters
        self.scf_initial = np.zeros(2)
        self.model = None
        self.model_set = False
        self.supported_models = ['phi2n', 'alpha_ads', 'phi2n_ads']
        self.supported_format = ['smg', 'ede']

    def set_model(self, model, format):
        if model in self.supported_models:
            self.model = model
        else:
            print('model %s not supported' % model)
            raise
        if format in self.supported_format:
            self.fmt = format
        else:
            print('format %s not supported' % format)
            raise
        self.model_set = True

    # input are standard LambdaCDM parameters: omega_m, T_cmb, H0
    def set_cosmo_params(
            self,
            omega_m=0.142,
            T0=2.7255,
            H0=70):
        self.omega_m = omega_m
        h = 0.01 * H0
        self.T_T0 = T0 / 2.7255
        # dedimensionalizing factor
        self.dd_fac = h * h * 8592.3

    def set_MCMC_params(
            self,
            scf_f=0.1,
            ln1plusz_c=8.2,
            smg0=0,  # three smg parameters
            smg1=0,
            smg2=0,
            v0=0,  # four potential parameters
            v1=0,
            v2=0,
            v3=0):
        self.scf_f = scf_f
        self.z_c = expm1(ln1plusz_c)
        self.smg_params[0] = smg0
        self.smg_params[1] = smg1
        self.smg_params[2] = smg2
        self.scf_params[0] = v0
        self.scf_params[1] = v1
        self.scf_params[2] = v2
        self.scf_params[3] = v3
        if self.verbose > 0:
            print("->Want scf_f = %.3f at z = %.1f" % (self.scf_f, self.z_c))

    def search(self):
        Hc2 = (self.omega_m * (self.z_c + 1) ** 3 + 0.000042 * (self.T_T0 * (self.z_c + 1)) ** 4) / 2998. / 2998.
        # get the initial guess for V_0 and phi_i
        if self.model == 'alpha_ads':
            # v0->V0, v1->gamma, v2->root6alpha, v3->Vlambda
            # In smg we force Lambda filling, thus v3=0
            self.scf_params[0] = 3 * self.scf_f * Hc2 * exp(-2 * self.scf_params[1])
            # initial point of root finding
            root_guess = -self.scf_params[2] * atanh(log1p(sqrt(0.5) * exp(self.scf_params[1])) / self.scf_params[1])
            # find the minimum of ddV
            phi_min = minimize_scalar(self.ddV, bracket=[-10, root_guess]).x
            self.scf_initial[0] = phi_min * 1.05
            ddV_min = self.ddV(phi_min)
        elif self.model == 'phi2n':  # use the same method as arXiv:1904. , adapted to implementation in EDE hi_CLASS, no shooting
            # v0->V0, v1->n, v2->Vads, v3->Vlambda
            rhoeq = 4124.0 / self.dd_fac
            zeq = 3394.8
            rhoc = rhoeq * (0.5 * (self.z_c / zeq) ** 3 + 0.5 * (self.z_c / zeq) ** 4)
            n = self.scf_params[1]
            phi_i = - sqrt(0.66667 * n * (2 * n - 1) * self.scf_f)
            v0 = self.scf_f * rhoc
            if self.fmt == 'smg':
                CLASS_params = np.array([self.smg_params[0], self.smg_params[1], self.smg_params[2], phi_i, 0.0, v0, n, 0, 0, 0.5])
            elif self.fmt == 'ede':
                CLASS_params = np.array(
                    [self.scf_params[0], self.scf_params[1], self.scf_params[2], self.scf_params[3], 0.5, self.scf_initial[0],
                     self.scf_initial[1]])
            return ', '.join(map(str, CLASS_params)), log1p(self.z_c)
        elif self.model == 'phi2n_ads':
            # v0->V0, v1->n, v2->Vads, v3->Vlambda
            ads = self.scf_params[2]
            n = self.scf_params[1]
            self.scf_params[0] = 3 ** (n + 1) * Hc2 / (2 * n * (2 * n - 1)) ** n / (self.scf_f + ads) ** (n - 1)
            self.scf_initial[0] = - sqrt(0.66667 * n * (2 * n - 1) * (self.scf_f + ads))
            self.scf_params[2] = 3. * Hc2 * ads

        if self.verbose > 0:
            print("->Initial guess V0 = " + str(self.scf_params[0]) + ", phi_i = " + str(self.scf_initial[0]) + ", Hc2 = " + str(Hc2))

        x0 = np.zeros(2)
        # V0
        x0[0] = self.scf_params[0]
        # phi_i
        x0[1] = self.scf_initial[0]
        ic_search_rlt = root(self.is_ic, x0, tol=self.precision)
        if not ic_search_rlt.success:
            return 'skip'
        self.scf_params[0] = ic_search_rlt.x[0]
        self.scf_initial[0] = ic_search_rlt.x[1]
        # Loop over all models, converting to parameters_smg in CLASS
        vflag = None
        if self.model == 'alpha_ads':
            vflag = 1.5
            print('Warning: alphaAdS not implemented in smg!')
        elif self.model == 'phi2n_ads':
            vflag = 2.5
        if vflag == None:
            print('model not recognized')
            raise
        if self.fmt == 'smg':
            CLASS_params = np.array(
                [self.smg_params[0], self.smg_params[2], self.scf_initial[0], self.scf_initial[1], self.scf_params[0],
                 self.scf_params[1]/ self.dd_fac, self.scf_params[2]/ self.dd_fac, self.scf_params[3], vflag])
        elif self.fmt == 'ede':
            CLASS_params = np.array(
                [self.scf_params[0], self.scf_params[1], self.scf_params[2], self.scf_params[3], vflag, self.scf_initial[0],
                 self.scf_initial[1]])

        # print(CLASS_params)
        # print(', '.join(map(str, CLASS_params)))
        return ', '.join(map(str, CLASS_params))

    def V(self, phi):
        V0 = self.scf_params[0]
        if self.model == 'alpha_ads':
            # V=V0*(exp(gamma*tanh(phi/root6alpha))-1)^2-(V0-V_Lambda)
            gamma = self.scf_params[1]
            root6alpha = self.scf_params[2]
            Vlambda = self.scf_params[3]
            return V0 * (expm1(gamma * tanh(-phi / root6alpha))) ** 2 - V0 + Vlambda
        elif self.model == 'phi2n':
            n = self.scf_params[1]
            return V0 * phi ** (2 * n)
        elif self.model == 'phi2n_ads':
            n = self.scf_params[1]
            Vads = self.scf_params[2]
            return V0 * phi ** (2 * n) - Vads

    def dV(self, phi):
        V0 = self.scf_params[0]
        if self.model == 'alpha_ads':
            gamma = self.scf_params[1]
            root6alpha = self.scf_params[2]
            return -2 * exp(- gamma * tanh(phi / root6alpha)) * expm1(
                -(gamma * tanh(phi / root6alpha))) * gamma * V0 / cosh(
                phi / root6alpha) / cosh(phi / root6alpha) / root6alpha
        elif self.model == 'phi2n_ads' or self.model=='phi2n':
            n = self.scf_params[1]
            Vads = self.scf_params[2]
            return 2 * n * V0 * phi ** (2 * n - 1)

    def ddV(self, phi):
        V0 = self.scf_params[0]
        if self.model == 'alpha_ads':
            gamma = self.scf_params[1]
            root6alpha = self.scf_params[2]
            return (-2 * gamma * V0 * pow(cosh(phi / root6alpha), -4) * (
                    (-1 + expm1(gamma * tanh(phi / root6alpha))) * gamma + expm1(
                gamma * tanh(phi / root6alpha)) * sinh((2 * phi) / root6alpha))) / (
                           exp(2 * gamma * tanh(
                               phi / root6alpha)) * root6alpha * root6alpha)
        elif self.model == 'phi2n_ads' or self.model == 'phi2n':
            n = self.scf_params[1]
            Vads = self.scf_params[2]
            return 2 * n * (2 * n - 1) * V0 * phi ** (2 * n - 2)

    def H2(self, t, scf_state):
        a_rel = (self.z_c + 1) * exp(-t)
        rho_m = 3 * self.omega_m * a_rel * a_rel * a_rel / 2998. / 2998.
        rho_r = 3 * 0.000042 * (a_rel * self.T_T0) ** 4 / 2998. / 2998.
        return (2 * (rho_r + rho_m + self.V(scf_state[0]))) / (
                6. - scf_state[1] * scf_state[1])

    # the EoM phi'' + (rho-P)/2/H^2*phi' + dV/dphi/H^2 = 0
    # x1 = phi,x2 = phi'
    def equation(self, t, scf_state):
        a_rel = (self.z_c + 1) * exp(-t)
        rho_m = 3 * self.omega_m * a_rel * a_rel * a_rel / 2998. / 2998.
        rho_r = 3 * 0.000042 * (a_rel * self.T_T0) ** 4 / 2998. / 2998.
        x_prime = np.zeros(2)
        x_prime[0] = scf_state[1]
        # rho - P = rho_m + 2/3 rho_r + 2V
        x_prime[1] = - (
                (rho_m + 2. / 3. * rho_r + 2 * self.V(scf_state[0])) / 2. * scf_state[
            1] + self.dV(
            scf_state[0])) / self.H2(
            t, scf_state)
        return x_prime

    def is_ic(self, x):
        rlt = np.zeros(2)
        phi_c = self.scf_initial[0]

        V0_original = self.scf_params[0]
        self.scf_params[0] = x[0]
        x0 = np.array([x[1], 0.])
        # time argument t = ln(a/a_c)
        scf_rlt = solve_ivp(
            self.equation, [
                self.t_i, 0], x0, method='RK23', t_eval=[0]).y.flatten()
        # rolling condition ddV_c=9H_c^2 is equivalent to phi=phi_c, exact ic at rlt[1]=1
        if self.model == 'alpha_ads':
            rlt[1] = scf_rlt[0] / phi_c
        elif self.model == 'phi2n_ads':
            rlt[1] = self.ddV(scf_rlt[0]) / self.H2(0, scf_rlt) / 9.
        # scalar field energy fraction rho_c=3fH_c^2, exact ic at rlt[0]=1
        rlt[0] = (scf_rlt[1] * scf_rlt[1] / 2. + self.V(scf_rlt[0]) / self.H2(0, scf_rlt)) / 3. / self.scf_f
        self.scf_params[0] = V0_original

        # exact ic at [0,0]
        return rlt - 1

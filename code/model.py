"""
Created on Mon Aug 19 2019

@author: benzekry
"""
import numpy as np
import pandas as pd
import carcinom_python.models_tumor_growth
import carcinom_python.utils
import scipy
import scipy.signal
import time
import pdb
#------------------------------------------------------------------------------
# Global parameters
#------------------------------------------------------------------------------
global_parameters                      = dict()
global_parameters['median_DT']         = 48 # (hours) median doubling time of neuroblastoma cells collected by Coline from meta-analysis in cellosaurus
#------------------------------------------------------------------------------
# Data
#------------------------------------------------------------------------------
class Data:
    '''
    Object containing data of a patient
    '''
    diam_1            = 146 # (mm) data of first patient from database
    diam_2            = 130 # (mm) data of first patient from database
    diam_3            = 115 # (mm) data of first patient from database
    S_diag_mm3        = 4/3*np.pi*diam_1/2*diam_2/2*diam_3/2 # (mm3)
    S_diag            = carcinom_python.utils.vol2cell(S_diag_mm3)
    SIOPEN            = 51
    LDH               = 1183
    def parse_df(self, df, idx):
        '''
        Parse data from dataframe to data object
        '''
        self.diam_1     = df.loc[idx, 'T1(mm)']
        self.diam_2     = df.loc[idx, 'T2(mm)']
        self.diam_3     = df.loc[idx, 'T3(mm)']
        if np.isnan(self.diam_3): # if only two diameters, assume the smallest of the first two as third value
            self.diam_3 = min(self.diam_1, self.diam_2)
        self.S_diag_mm3 = 4/3*np.pi*self.diam_1/2*self.diam_2/2*self.diam_3/2 # (mm3)
        self.S_diag     = carcinom_python.utils.vol2cell(self.S_diag_mm3) # (cells)
        self.SIOPEN     = df.loc[idx, 'SIOPEN']
        self.LDH        = df.loc[idx, 'LDH']
#------------------------------------------------------------------------------
# Parameters
#------------------------------------------------------------------------------
class Parameters:
    alpha0            = np.log(2)/global_parameters['median_DT']*24 # (day-1) proliferation rate at one cell
    # beta              = alpha0/np.log(carrying_capacity)
    mu                = 9.969e-12
    gamma             = 1
    sigma_SIOPEN      = 0.1 # standard deviation in error model of SIOPEN
    sigma_LDH         = 0.1 # standard deviation in error model of LDH
    visible_threshold = 72.2   # (mm)
    phi               = 1e-10
#------------------------------------------------------------------------------
# Model
#------------------------------------------------------------------------------
class Model:
    '''
    Model for description of SIOPEN and LDH in neuroblastoma from total and visible metastatic burden
    '''
    def __init__(self, parameters=Parameters, data=Data):
        self.parameters                 = parameters
        self.data                       = data
        self.param_opt_indiv_names_to_idx     = {
                                           'mu': 0,
                                           } # gives index of estimated parameter in the param_opt array. relative to param_opt_indiv
        self.param_opt_pop_names_to_idx = {'visible_threshold': 0} # relative to param_opt_pop
        self.set_param_opt_idx_to_names()
        self.set_T_diag(self.data)
    def set_param_opt_idx_to_names(self):
        self.param_opt_indiv_idx_to_names = {v: k for k, v in self.param_opt_indiv_names_to_idx.items()}
        self.param_opt_pop_idx_to_names   = {v: k for k, v in self.param_opt_pop_names_to_idx.items()}
    def set_param_pop(self, param_pop):
        '''
        Set population-level parameters in parameters object
        '''
        for param_name in self.param_opt_pop_names_to_idx.keys():
            param_idx = self.param_opt_pop_names_to_idx[param_name]
            setattr(self.parameters, param_name, param_pop[param_idx])
    def set_param_indiv(self, param_indiv):
        '''
        Set individual-level parameters in parameters object
        '''
        for param_name in self.param_opt_indiv_names_to_idx.keys():
            param_idx = self.param_opt_indiv_names_to_idx[param_name]
            setattr(self.parameters, param_name, param_indiv[param_idx])
    def set_T_diag(self, data):
        self.T_diag = self.time_to_size_PT(parameters=self.parameters, sizes= data.S_diag) # (day) taken from exponential model because gives "realistic" estimates (to be double-checked with the literature)
    #------------------------------------------------------------------------------
    # Estimation of tumor age
    #------------------------------------------------------------------------------
    def time_to_size_PT(self, parameters, sizes):
        alpha0 = parameters.alpha0
        times  = np.log(sizes)/alpha0
        return times
    #------------------------------------------------------------------------------
    # Metastatic burden at diagnosis
    #------------------------------------------------------------------------------
    def metastatic_burden(self,
                          parameters,
                          T_diag,
                          num=1001,
                          method='scipy_fft'):
        '''
        Metastatic burden at diagnosis through convolution formula (assumes no secondary dissemination)
        '''
        mu                 = parameters.mu
        gamma              = parameters.gamma
        time_vect          = np.linspace(start=0, stop=T_diag, num=num)
        dt                 = time_vect[1] - time_vect[0]
        f                  = self.PT_growth_model(parameters, time_vect)**gamma
        g                  = self.met_growth_model(parameters, time_vect)
        L                  = 2*num - 1
        if method is 'integral':
            integrand = lambda t: mu*(self.PT_growth_model(parameters, T_diag - t)**gamma)*PT_growth_model(parameters, t)
            M, err    = scipy.integrate.quad(integrand, 0, T_diag)
        elif method is 'convolve':
            convs = np.convolve(f, g)
            M     = mu*dt*convs[num]
        elif method is 'fft':
            Ff                 = np.fft.fft(f, L) # f has to be prolonged with zero values in order to compute convolution afterwards because support of f*g = supp(f) + supp(g), i.e. len(f*g) = len(f) + len(g) -1
            Fg                 = np.fft.fft(g, L) # same for g
            convs              = np.fft.ifft(Ff*Fg)
            M                  = mu*dt*convs[num]
        elif method is 'scipy_fft':
            convs              = scipy.signal.fftconvolve(f, g)
            M                  = mu*dt*convs[num]
        return M
        
    def metastatic_number(self,
                          parameters,
                          T):
        '''
        Computes number of metastases at time T
        N(T) = \mu\int_0^T S_p(t) dt
        '''
        mu     = parameters.mu
        gamma  = parameters.gamma
        Sp_f   = lambda t: mu*self.PT_growth_model(parameters, t)**gamma
        N, err = scipy.integrate.quad(Sp_f, 0, T)
        return N
        
    def met_growth_model(self, parameters, time_vect):
        '''
        Growth model for metastases
        '''
        # S = carcinom_python.models_tumor_growth.gompertz_function(parameters, time_vect, tI=0, VI=1, Vc=1) # (cells)
        S = np.exp(parameters.alpha0*time_vect)
        return S

    def PT_growth_model(self, parameters, time_vect):
        '''
        Growth model for primary tumor
        '''
        # Sp = carcinom_python.models_tumor_growth.gompertz_function(parameters, time_vect, tI=0, VI=1, Vc=1) # (cells)
        Sp = np.exp(parameters.alpha0*time_vect)
        return Sp
        
    def time_to_size_met(self, parameters, sizes):
        alpha0 = parameters.alpha0
        times  = np.log(sizes)/alpha0
        return times
    #------------------------------------------------------------------------------
    # Visible metastatic burden at diagnosis
    #------------------------------------------------------------------------------
    def metastatic_burden_visible(self,
                                  parameters,
                                  T,
                                  num=1001):
        '''
        Visible metastatic burden at time T
        M_vis (T)  = \mu \int_0^{T - \tau_vis} S_p(T - \tau_vis - t) S(t + \tau_{vis}) dt
                   (= \mu \int_{\tau_vis}^{T} S_p(T - t) S(t) dt, but cannot be used for discrete convolution calculation)
        '''
        mu                     = parameters.mu
        gamma                  = parameters.gamma
        # Time to reach visible threshold
        visible_threshold_cell = carcinom_python.utils.diam2cell(parameters.visible_threshold)
        tau_vis                = self.time_to_size_met(parameters, visible_threshold_cell)
        if tau_vis > T:
            Mvis = 0
            return Mvis
        # Functions for convolution
        time_vect              = np.linspace(start=0, stop=T-tau_vis, num=num)
        dt                     = time_vect[1] - time_vect[0]
        f                      = self.PT_growth_model(parameters, time_vect)**gamma
        g                      = self.met_growth_model(parameters, tau_vis + time_vect)
        L                      = 2*num - 1
        convs                  = scipy.signal.fftconvolve(f, g)
        Mvis                   = mu*dt*convs[num]
        return Mvis
        
    def metastatic_number_visible(self,
                                  parameters,
                                  T):
        '''
        Number of visible metastasis at time T
        N_vis(T) = \int_{S_vis}^{+\infty} \rho(T, s) ds
                 = N(T-tau_vis)
        '''
        # Time to reach visible threshold
        visible_threshold_cell = carcinom_python.utils.diam2cell(parameters.visible_threshold)
        tau_vis                = self.time_to_size_met(parameters, visible_threshold_cell)
        if tau_vis > T:
            N_vis = 0
            return N_vis
        N_vis = self.metastatic_number(parameters, T - tau_vis)
        return N_vis
        
    def model_SIOPEN(self,
                     parameters):
        '''
        Model SIOPEN = number of visible mets
        '''
        N_vis        = self.metastatic_number_visible(parameters, self.T_diag)
        model_SIOPEN = N_vis
        return model_SIOPEN
        
    def model_LDH(self,
                  parameters,
                  data):
        M            = self.metastatic_burden(parameters, self.T_diag)
        model_LDH    = parameters.phi*(M + data.S_diag)
        return model_LDH
#------------------------------------------------------------------------------
# Log-likelihood
#------------------------------------------------------------------------------
def loglik_SIOPEN_LDH(
                     param_opt,           # dict
                     data,                # data object
                     model,               # model object
                     bounds=None,              # array of (low, up) tuples of length P
                     output='full'):
    '''
    Opposite of loglikelihood to be minimized for fitting the model to SIOPEN and LDH data at diagnosis. Proportional error model
    SIOPEN = N_vis(\theta, S) * (1 + \sigma1 \varepsilon_1)
    LDH    = phi*(M(\theta, S) + \phi S)*(1 + \sigma2 \varepsilon_2)
    '''
    if not bounds is None:
        # Set objective value to inf if parameters are outside bounds
        low_bounds = [low for (low, up) in bounds]
        up_bounds  = [up for (low, up) in bounds]
        if any(param_opt < low_bounds) | any(param_opt > up_bounds):
            l = np.inf
            return l
    parameters = model.parameters
    # Update values of parameters that are to be estimated
    for param_name in model.param_opt_indiv_names_to_idx:
        param_value = param_opt[model.param_opt_indiv_names_to_idx[param_name]]
        setattr(parameters, param_name, param_value)
    # Compute loglikelihood
    sigma1       = parameters.sigma_SIOPEN
    sigma2       = parameters.sigma_LDH
    model_SIOPEN = model.model_SIOPEN(parameters)
    model_LDH    = model.model_LDH(parameters, data)
    M            = model.metastatic_burden(parameters, model.T_diag)
    if (sigma1*model_SIOPEN < 0) | (sigma2*model_LDH < 0):
        pdb.set_trace()
    if model_SIOPEN > 0: # if no visible metastasis, divide by SIOPEN to avoid divide by zero
        l_SIOPEN = (data.SIOPEN - model_SIOPEN)**2/(2*(sigma1*model_SIOPEN)**2) + np.log(sigma1*model_SIOPEN) + 1/2*np.log(2*np.pi)
    elif data.SIOPEN == 0:
        l_SIOPEN = 0
    else:
        l_SIOPEN = (data.SIOPEN - model_SIOPEN)**2/(2*(sigma1*data.SIOPEN)**2) + np.log(sigma1*data.SIOPEN) + 1/2*np.log(2*np.pi)
    l_LDH    = (data.LDH - model_LDH)**2/(2*(sigma2*model_LDH)**2) + np.log(sigma2*model_LDH) + 1/2*np.log(2*np.pi)
    l        =  l_SIOPEN + l_LDH
    if output is 'full':
        return l, l_SIOPEN, l_LDH, model_SIOPEN, M, model_LDH
    else:
        return l
#------------------------------------------------------------------------------
# Fit individuel
#------------------------------------------------------------------------------
def fit_SIOPEN_LDH(
                   data,              # data object
                   model,             # model object
                   param_opt_0,       # array
                   disp=True):
    '''
    Individual fit of SIOPEN and LDH data as visible and total metastatic burden, respectively
    '''
    model.set_T_diag(data)
    bounds = [(0, np.inf) for param in param_opt_0]
    res    = scipy.optimize.minimize(
                    loglik_SIOPEN_LDH,
                    x0=param_opt_0,
                    args=(data, model, bounds, 'simple'),
                    method='Nelder-Mead',
                    options={'disp':disp})
    param_fit = res.x
    l_fit     = res.fun
    return param_fit, l_fit
#------------------------------------------------------------------------------
# Fit population
# All patients share some parameters, which are still estimated
#------------------------------------------------------------------------------
def fit_SIOPEN_LDH_pop(
                       df,                # dataframe
                       model,             # model object
                       param_opt_pop_0,   # array
                       param_opt_indiv_0,  # array
                       disp_iter=True
                       ):
    '''
    Fit all patients simultaneously, with some parameters equal for all and some patient-specific
    '''
    N                  = len(df)
    param_opt_indivs_0 = np.tile(param_opt_indiv_0, N)
    param_opt_0        = np.concatenate((param_opt_pop_0, param_opt_indivs_0))
    bounds             = [(0, np.inf) for param in param_opt_0]
    if disp_iter:
        global Niter
        Niter              = 1 # for display
        callbackF_loc = lambda x: callbackF(x, df, model, bounds)
    else:
        callbackF_loc = lambda x: None
    res                = scipy.optimize.minimize(sum_logliks,
                                                x0=param_opt_0,
                                                args=(df, model, bounds),
                                                method='Nelder-Mead',
                                                callback=callbackF_loc,
                                                options={'disp':True,
                                                         'maxiter':40000}
                                                )
    params_fit = res.x
    l_pop_fit  = res.fun
    param_fit_pop, params_fit_indivs = parse_param_pop_all(params_fit, df, model)
    return param_fit_pop, params_fit_indivs, res
def sum_logliks(param_opt, # array of length P_pop + N*P_indiv where P_pop and P_indiv are the number of population- and individual-level parameters
                df,
                model,
                bounds=None
                ):
    '''
    Returns the sum of the individual logliks for some population-level parameters and an array of individual-level parameters
    '''
    if bounds is not None:
        # Set objective value to inf if parameters are outside bounds
        low_bounds = [low for (low, up) in bounds]
        up_bounds  = [up for (low, up) in bounds]
        if any(param_opt < low_bounds) | any(param_opt > up_bounds):
            loglik_pop = np.inf
            return loglik_pop
    # Parse parameter values
    param_opt_pop, param_opt_indivs_mat = parse_param_pop_all(param_opt, df, model)
    # Update values of parameters that are to be estimated
    parameters                 = model.parameters
    model.set_param_opt_idx_to_names()
    for idx, param_value in enumerate(param_opt_pop):
        # pdb.set_trace()
        param_name = model.param_opt_pop_idx_to_names[idx]
        setattr(parameters, param_name, param_value)
    model.parameters = parameters # updates value in model object so that it gets passed to individual likelihood maximization
    # Compute sum of logliks
    loglik_pop  = 0
    df_fit      = pd.DataFrame(index=df.index)
    for idx_enum, idx_patient in enumerate(df.index): # enumerate because index might have missing values for patients removed or if not continuous
        data                      = Data()
        data.parse_df(df, idx_patient)
        param_opt_indiv           = param_opt_indivs_mat[idx_enum, :]
        loglik_indiv              = loglik_SIOPEN_LDH(param_opt_indiv, data, model, output="single")
        loglik_pop               += loglik_indiv
    return loglik_pop
def parse_param_pop_all(param_pop_all, df, model):
    '''
    Parse parameter values given in a P_pop + N*P_indiv 1D array into a P_pop array and a NxP_indiv array
    '''
    N                          = len(df.index) # number of patients
    P_pop                      = len(model.param_opt_pop_names_to_idx)
    P_indiv                    = len(model.param_opt_indiv_names_to_idx)
    param_pop              = param_pop_all[0:P_pop]
    param_indivs_vect      = param_pop_all[P_pop:]
    param_indivs_mat       = param_indivs_vect.reshape(N, P_indiv)
    return param_pop, param_indivs_mat
def callbackF(x, df, model, bounds):
    '''
    Callback function for display during minimization process
    '''
    global Niter
    if (Niter % 1) == 0:
        print('Iter = {:4d}'.format(Niter))
        print('Loglik = {:4g}'.format(sum_logliks(x, df, model, bounds)))
    Niter += 1
#------------------------------------------------------------------------------
# Fit population with other method (min(sum(min)))
#------------------------------------------------------------------------------
def fit_SIOPEN_LDH_pop_sum_min(
   df,                # dataframe
   model,             # model object
   param_opt_pop_0,   # array
   param_opt_indiv_0,  # array
   disp_iter=True
   ):
   '''
   Fit population with other method (min(sum(min)))
   '''
   N                  = len(df)
   bounds             = [(0, np.inf) for param in param_opt_pop_0]
   if disp_iter:
       global Niter
       Niter              = 1 # for display
       callbackF_loc = lambda x: callbackF_sum_min(x, df, model, param_opt_indiv_0, bounds)
   else:
       callbackF_loc = lambda x: None
   res                = scipy.optimize.minimize(sum_min_logliks,
                                               x0=param_opt_pop_0,
                                               args=(df, model, param_opt_indiv_0, bounds),
                                               method='Nelder-Mead',
                                               callback=callbackF_loc,
                                               options={'disp':True
                                                        # 'maxiter':40000
                                                        }
                                               )
   param_fit_pop = res.x
   l_pop_fit     = res.fun
   return param_fit_pop, l_pop_fit
def sum_min_logliks(
        param_opt_pop,
        df,
        model,
        param_opt_indiv_0,
        bounds=None
        ):
    '''
    Returns sum of the minimal loglikelihood fitted individually
    '''
    # Set objective value to inf if parameters are outside bounds
    if bounds is not None:
       low_bounds = [low for (low, up) in bounds]
       up_bounds  = [up for (low, up) in bounds]
       if any(param_opt_pop < low_bounds) | any(param_opt_pop > up_bounds):
           loglik_pop = np.inf
           return loglik_pop
    # Update values of parameters that are to be estimated
    model.set_param_pop(param_opt_pop)
    # Compute sum of min of logliks
    loglik_pop  = 0
    df_fit      = pd.DataFrame(index=df.index)
    for idx_enum, idx_patient in enumerate(df.index): # enumerate because index might have missing values for patients removed or if not continuous
       data                      = Data()
       data.parse_df(df, idx_patient)
       params_fit, loglik_indiv  = fit_SIOPEN_LDH(data, model, param_opt_indiv_0, disp=False)
       loglik_pop               += loglik_indiv
       # print(idx_patient)
    return loglik_pop
def callbackF_sum_min(x, df, model, param_opt_indiv_0, bounds):
    '''
    Callback function for display during minimization process
    '''
    global Niter
    if (Niter % 1) == 0:
        print('Iter = {:4d}'.format(Niter))
        print('Loglik = {:4g}'.format(sum_min_logliks(x, df, model, param_opt_indiv_0, bounds)))
    Niter += 1

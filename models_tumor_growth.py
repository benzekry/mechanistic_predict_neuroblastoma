# -*- coding: utf-8 -*-
"""
Created on Mon Aug  3 17:16:10 2015

@author: benzekry
"""
import sys
sys.path.insert(1, 'levenberg')
sys.path.insert(1, 'models')
from numpy import *
import numpy as np
from numpy.linalg import norm
from scipy.integrate import odeint
import pdb
import warnings

# Global variables
"By default, volumes are assumed to be in mm3"
# ---> Specify here the volume of one cell if the models you use are based on this value
Vc = 1e-6 # (mm^3) volume of one cell based on 1 mm^3 = 10^6 cells
# ---> Specify here the initial volume (such as the number of injected cells)
V0 = 1. # (mm^3)
# ---> Specify here the in vitro doubling time
DT_invitro = array([35.9]) # (hours) doubling time in vitro of LLC cells (to be checked)
############################################################################################
### Class Model
############################################################################################
class Model:
    """
    The class Model contains the necessary attributes for a fit to be performed
    and the parameters to be estimated
    fit_method [string] = model_key of the optimization procedure to be used in the fit
    folder [string] = folder where the export is stored
    IC [string] = Indicates how to treat initial condition:
        - 'fixed' : the value at t = 0 is set to a fixed V0
        - 'V0' : the value at t = 0 is a free parameter V0
        - 'VI' : the initial condition is set to the first data point
    max_function_evals [double] = maximum number of function evaluations to be allowed in the fit
    max_iter [double] = maximum number of iterations to be allowed in the fit
    model_name [string] = model_key of the model to be displayed in output tables
    model_function = function of arguments
        - (param,time) if IC = 'fixed'
        - (param, time) if IC =  'V0'. CAUTION: param[-1] must be V0.
        - (param, time, tI, VI) if IC = 'VI'
    param0 [array] = array of the initial values of the parameters for the fit
    param_fixed [dict] = parameters not subject to the fit
    param_ap [array] = a priori parameters distribution for bayesian estimation
    param_ap_ind [dict] = parameters to be estimated using an a priori estimate
    units [list of strings (LaTeX allowed)] = units of the parameters
    xlabel [string] = label to be used in the plots for the x axis
    ylabel [string] = label to be used in the plots for the y axis
    """
    number_of_models = 0
    def __init__(self, model_key = 'model', model_name = '', param_names =[], units = ['-']*100,
        model_function = None, IC_tag = 'fixed', V0 = V0, param0 = None, low_bounds = [], up_bounds = []):
        """
        model_function = function of arguments (param, time, tI, VI)
        """
        Model.number_of_models +=1
        self.fit_method  = "Nelder-Mead"
        self.IC          = False
        self.param0      = param0
        self.param_names = param_names
        self.units       = units
        self.low_bounds  = low_bounds
        self.up_bounds   = up_bounds
        if IC_tag == 'fixed':
            self.folder         = model_key
            self.model_function = (lambda param, time: model_function(param, time, 0, V0))
            self.model_name     = model_name
        elif IC_tag == 'VI':
            self.folder         = model_key+'_VI'
            self.IC             = True
            self.model_function = model_function
            self.model_name     = model_name+' $V_I$'
        elif IC_tag == 'V0':
            self.folder = model_key+'_V0'
            self.model_function = (lambda param, time:  model_function(param, time, 0, param[-1]))
            self.param0 = append(param0, array([V0]))
            self.param_names = append(param_names, '$V_0$')
            self.units = append(units, '$mm^3$')
            self.model_name = model_name+' $V_0$'
            if len(low_bounds) > 0:
                self.low_bounds = append(low_bounds,0.)
            if len(up_bounds) > 0:
                self.up_bounds = append(up_bounds,100*V0)
        self.max_function_evals = 0 # 0 = no bound
        self.max_iter = 0 # 0 = no bound
        self.model_no_IC = ""
        self.param_fixed = {}
        self.param_ap = []
        self.param_ap_ind = []
        self.units = ['-']*100
        self.xlabel = 'Time (days)'
        self.ylabel = 'Volume (mm$^3$)'
        #Attributes for the L-M algorithm
        self.leven_parameters = None

###########################################################################
### Dynamic CC
###########################################################################
def dynamic_cc_function(param, time, tI, VI):
    """
    Model with dynamic carrying capacity inspired by [Hahnfeldt et al., Cancer Res, 1999]
    The form used here is due to B. Ribba (see ref in [Benzekry et al., PloS Comp Biol 2014])
    Initial carrying capacity KI has to be provided as the last entry of vector param
    """
    XI = array([VI, param[-1]])
    X = dynamic_cc_full(param, time, tI, XI)
    V = X[:,0]
    return V
def dynamic_cc_full(param, time, tI, XI):
    epsi = 1e-10
    if norm(time - tI)<epsi:
        X = XI
    else:
        time_old = time
        if abs(time[0]-tI) > epsi:
            time = append(array([tI]), time)
        X = odeint(dynamic_cc_rhs, XI, time, args = tuple(array([param])))
        if abs(time_old[0]-tI)>epsi:
            X = X[1:,:]
    return X
def dynamic_cc_rhs(X, t, param):
    a = param[0]
    b = param[1]
    d = param[2]
    V = X[0]
    K = X[1]
    dX = zeros(2)
    dX[0] = a*V*log(K/V)
    dX[1] = b*V**(2/3)
    return dX
class DynamicCC(Model):
    def __init__(self, IC_tag = 'fixed'):
        Model.__init__(self,
            IC_tag = IC_tag,
            model_key = 'dynamic_cc',
            model_function = (lambda param, time, tI, VI:
                dynamic_cc_function(append(param[0:2],array([10])), time, tI, VI)),
            model_name = 'Dynamic CC',
            param0 = [3., 0.5],
            #low_bounds = [0, 0],
            param_names = ['$a$', '$b$'], #, '$K_0$']
            units = ['$\\left[day^{-1}\\right]$','$\\left[mm^{-2}\\cdot day^{-1}\\right] $']#,'$\\left[mm^3\\right]$']
            )
        #self.max_function_evals = 1000
###########################################################################
### Exponential
###########################################################################
class Exponential(Model):
    """
    Exponential growth
    """
    def __init__(self, IC_tag = 'fixed'):
        Model.__init__(self,
            IC_tag         = IC_tag,
            model_key      = 'exponential',
            model_function = (lambda param, time, tI, VI: VI*exp(param[0]*(time-tI))),
            model_name     = 'Exponential',
            param0         = [0.1],
            param_names    = ['$a$'],
            units          = ['$\\left[day^{-1}\\right]$']
            )
def exponential_time_to_size(params, size):
    '''
    Computes time to reach a given size when growing following exponential growth starting from one cell
    '''
    alpha = params[0]
    t     = np.log(size)/alpha
    return t
###########################################################################
### Exponential-linear
###########################################################################
def exponential_linear_function(param, temps, tI, VI):
    """
    Exponential-linear  Switch time is determined by contuinity of derivative
    """
    lambda0 = param[0]; lambda1 = param[1];
    temps = temps - tI
    tau = (1/lambda0*log(lambda1/(VI*lambda0)))
    wtau = VI*exp(lambda0*tau)
    y = VI*exp(lambda0*temps)*(temps<tau)+(lambda1*(temps-tau)+wtau)*(temps>=tau)
    if isinf(y).any():
        warnings.warn('Exponential-linear model returned inf value');
    return y
class ExponentialLinear(Model):
    """
    Bi-phasic growth with first an exponential phase, followed by a linear one.
    Proposed by Simeoni et al. in Cancer Res., 2004
    """
    def __init__(self, IC_tag = 'fixed'):
        Model.__init__(self,
            IC_tag = IC_tag,
            model_key = 'exponential_linear',
            model_function = exponential_linear_function,
            model_name = 'Exponential-linear',
            param0 = [0.2, 100],
            low_bounds = [0., 0.],
            up_bounds = [10, 300],
            param_names = ['$a_0$','$a_1$'],
            units = ['$\\left[day^{-1}\\right]$', '$\\left[mm^{3}\\cdot day^{-1}\\right]$']
            )
###########################################################################
### Generalized logistic
###########################################################################
def generalized_logistic_function(param, temps, tI, VI):
    """
    Generalized logistic model
    dV/dt = aV(1-(V/K)^alpha), V(t = tI) = VI
    """
    a = param[0]
    K = param[1]
    alpha = param[2]
    temps = temps - tI
    V=VI*K/(VI**alpha+(K**alpha-VI**alpha)*exp(-a*alpha*temps))**(1./alpha)
    return V
class GeneralizedLogistic(Model):
    def __init__(self, IC_tag = 'fixed'):
        Model.__init__(self,
            IC_tag = IC_tag,
            model_key = 'generalized_logistic',
            model_function = generalized_logistic_function,
            model_name = 'Generalized logistic',
            param0 = [10, 10000, 0.01],
            param_names = ['$a$', '$K$', '$\\alpha$']
            )
###########################################################################
### Gompertz
###########################################################################
def gompertz_function(param, temps, tI, VI, Vc):
    """
    Gompertzian growth according to
    dV/dt = (alpha0-beta*ln(V/Vc))*V, V(t = tI) = VI
    hence alpha0 = relative growth rate at V = 1 cell (volume Vc)
    alpha0 = param(1) beta=param(2)
    >>> alpha0 = 1
    >>> beta = 0.1
    >>> Vc = 1
    >>> V0 = 1
    >>> temps = 0.1*arange(0,3000)
    >>> V = gompertz_function(array([alpha0,beta]),temps,Vc,V0)
    >>> norm(V[0] - 1)<1e-10
    True
    >>> norm(V[-1] - Vc*exp(alpha0/beta))<1e-6
    True
    """
    alpha0 = param[0]
    beta = param[1]
    temps = temps - tI
    V = Vc*(VI/Vc)**(exp(-beta*temps))*exp(alpha0/beta*(1-exp(-beta*temps)))
    return V
class Gompertz(Model):
    def __init__(self, IC_tag = 'fixed'):
        Model.__init__(self,
            IC_tag = IC_tag,
            model_key = 'gompertz',
            model_function = (lambda param, time, tI, VI: gompertz_function(param, time, tI, VI, Vc)),
            model_name = 'Gompertz',
            param0 = array([0.71, 0.7]),
            param_names = ['$\\alpha_0$','$\\beta$'],
            units = ['$\\left[day^{-1}\\right]$','$\\left[day^{-1}\\right]$']
            )
def doubling_time_gompertz(param, time, Vc, V0):
    """
    Computation of the doubling time of the Gompertz model
    param = [alpha0, beta] parameters of Gompertz growth, according to
    dV/dt=(alpha0 - beta*log(V/Vc))*V
    hence alpha0 = relative growth rate at V = 1 cell (volume Vc)
    If alpha0 and beta are in day^-1 then DT is in days
    >>> Vc = 1e-6
    >>> V0 = 1
    >>> alpha0 = 2.3
    >>> beta = 0.1
    >>> K = Vc*exp(alpha0/beta)
    >>> DT = doublingTimeGompertzAlpha0(array([alpha0, beta]), 0, Vc, V0)
    >>> norm(DT - 0.7846) < 1e-4
    True
    """
    alpha = param[0]
    beta = param[1]
    A = log((V0/Vc)**(exp(-beta*time)))-alpha/beta*exp(-beta*time)
    ind_nonAd = where((log(2)+A)/A < 0)[0]
    DT = -1/beta*log((log(2)+A)/A) # cf ComputationsGompertz.pdf
    if len(ind_nonAd) > 0:
        warnings.warn('Some doubling time(s) are set to infinite because \
        the double size is larger than the maximal reachable size')
        DT[ind_nonAd] = Inf
    return DT
def gompertz_time_to_size(params, size):
    '''
    Computes time to reach a given size when growing following the Gompertz law
    dV/dt = (alpha0 - beta*log(V))*V
    WARNING: this function with analytical formula is only valid in the case
    V(t=0) (= Vc) = 1
    '''
    alpha0 = params[0]
    beta   = params[1]
    if alpha0/beta <= log(size):
        t = np.inf
    else:
        t = -1/beta*log(1 - beta/alpha0*log(size))
    return t
###########################################################################
### Gomp-Exp
###########################################################################
def gompertz_exp_function(param, temps, tI, VI, DT_invitro, Vc, return_switch_time = False):
    """
    GompExp model
    param = [alpha0, beta] parameters of Gompertz growth, according to
    dV/dt=(alpha-beta*log(V/Vc))*V
    hence alpha = relative growth rate at V=1 cell (volume Vc)
    DT_invitro = in vitro doubling time in hours
    >>> alpha0 = 1
    >>> beta = 0.1
    >>> Vc = 1
    >>> V0 = 1
    >>> temps = 0.1*arange(0,100)
    >>> V, tau = gompertzExp_function(array([alpha0,beta]),DT_invitro,temps,Vc,V0, return_switch_time = True)
    >>> norm(tau - 6.9)< 1e-8
    True
    >>> norm(V[-1] - 142.66)< 1e-2
    True
    """
    DT_invitro_day = DT_invitro/24 # conversion in days
    SGR = log(2)/DT_invitro_day
    temps = temps - tI
    DTs = doubling_time_gompertz(param, temps, Vc, V0)
    if count_nonzero(DTs > DT_invitro_day) == 0:
        warnings.warn('Warning')
        print('All Gompertz-doubling times are smaller than the in vitro one. '
            'There will be only exponential phase')
    ind_exp = where((DTs < DT_invitro_day)*(DTs>0))[0]
    temps_exp = temps[ind_exp]
    temps_gomp = temps
    temps_gomp = delete(temps_gomp, ind_exp)
    V_exp = V0*exp(SGR*temps_exp)
    if len(temps_exp) == 0:
        tI_gomp = 0
        VI_gomp = V0
    else:
        tI_gomp = temps_exp[-1]
        VI_gomp = V_exp[-1]
    V_gomp = gompertz_function(param, temps_gomp, tI_gomp, VI_gomp, Vc)
    V=hstack([V_exp,V_gomp])
    if return_switch_time == True:
        return V, tI_gomp
    else:
        return V
class GompExp(Model):
    def __init__(self, IC_tag = 'fixed'):
        Model.__init__(self,
            IC_tag = IC_tag,
            model_key = 'gompertz_exp',
            model_function = (lambda param, time, tI, VI: gompertz_exp_function(
                param, time, tI, VI, DT_invitro, Vc)),
            model_name = 'Gomp-Exp',
            param0 = array([1, 0.1]),
            param_names = ['$\\alpha$','$\\beta$']
            )
###########################################################################
### Hahnfeldt
###########################################################################
def hahnfeldt_function(param, time, tI, VI):
    """
    Hanhnfeldt model of tumor growth under angiogenic signaling [Hahnfeldt et al., Cancer Res 1999]
    Initial carrying capacity KI has to be provided as the last entry of vector param
    """
    XI = array([VI, param[-1]])
    X = hahnfeldt_full(param, time, tI, XI)
    V = X[:,0]
    return V
def hahnfeldt_full(param, time, tI, XI):
    """
    Fulll Hahnfeldt returning V and K
    """
    epsi = 1e-10
    if norm(time - tI)<epsi:
        X = XI
    else:
        time_old = time
        if abs(time[0]-tI) > epsi:
            time = append(array([tI]), time)
        X = odeint(hahnfeldt_rhs, XI, time, args = tuple(array([param])))
        if abs(time_old[0]-tI)>epsi:
            X = X[1:,:]
    return X
def hahnfeldt_rhs(X, t, param):
    """
    Right hand side of the Hahnfeldt ODE model
    """
    a = param[0]
    b = param[1]
    d = param[2]
    V = X[0]
    K = X[1]
    dX = zeros(2)
    dX[0] = a*V*log(K/V)
    dX[1] = b*V - d*V**(2./3)*K
    return dX
class Hahnfeldt(Model):
    def __init__(self, IC_tag = 'fixed'):
        Model.__init__(self,
        IC_tag         = IC_tag,
        model_key      = 'hahnfeldt',
        model_function = (lambda param, time, tI, VI: hahnfeldt_function(
            append(param[0:2], array([0.0717, 300])), time, tI, VI)),
        model_name     = 'Hahnfeldt',
        param0         = [0.1, 15],
        param_names    = ['$a$', '$b$'],
        )
        self.max_function_evals = 1000
###########################################################################
### Logistic
###########################################################################
def logistic_function(param, time, tI, VI):
    time = time - tI
    V = VI*param[1]/((param[1]-VI)*exp(-param[0]*time)+VI)
    return V
class Logistic(Model):
    def __init__(self, IC_tag = 'fixed'):
        Model.__init__(self,
            model_key = 'logistic',
            model_function = logistic_function,
            model_name = 'Logistic',
            param0 = [1, 10000],
            param_names = ['$a$','$K$'],
            units = ['$\\left[day^{-1}\\cdot mm^{3(1-\\gamma)}\\right]$', '-']
            )
###########################################################################
### Power law
###########################################################################
def power_law_function(param, temps, tI, VI):
    """
    Growth model following
    dV/dt = a V^\gamma, V(t = tI) = VI
    """
    a = param[0];
    gamma = param[1];
    temps = temps - tI
    if (gamma<0) | (gamma>1) | (a<0) | (VI<0):
         warnings.warn('Warning')
         print('In the power law model, either gamma < 0 or gamma > 1 or a < 0 or V_I < 0.'+
            ' Model is set to zero')
         V = zeros(temps.shape)
    else:
        V = (VI**(1-gamma)+a*(1-gamma)*temps)**(1./(1-gamma));
    if isinf(V).any() | isnan(V).any():
        warnings.warn('Warning')
        print('Power law model has returned Inf or NaN value. Set to zero')
        V = zeros(temps.shape);
    return V
class PowerLaw(Model):
    def __init__(self, IC_tag = 'fixed'):
        Model.__init__(self,
            IC_tag = IC_tag,
            model_key = 'power_law',
            model_function = power_law_function,
            model_name = 'Power law',
            param0 = array([1, 2./3]),
            param_names = ['$\\alpha$','$\\gamma$'],
            units = ['$\\left[day^{-1}\\cdot mm^{3(1-\\gamma)}\\right]$','-']
            )
###########################################################################
### Von Bertalanffy
###########################################################################
def von_bertalanffy_function(param, temps, tI, VI):
    """
    Von Bertalanffy model of growth
    dV/dt = aV^\gamma - bV, V(t = tI) = VI
    """
    a = param[0]
    gamma = param[1]
    b = param[2]
    temps = temps - tI
    V = (a/b+(VI**(1-gamma)-a/b)*exp(-b*(1-gamma)*temps))**(1/(1-gamma))
    return V
class VonBertalanffy(Model):
    def __init__(self, IC_tag = 'fixed'):
        Model.__init__(self,
            IC_tag = IC_tag,
            model_key = 'von_bertalanffy',
            model_function = von_bertalanffy_function,
            model_name = 'Von Bertalanffy',
            param0 = [1., 2./3, 0.1],
            param_names = ['$a$','$\\gamma$','$b$']
        )
###########################################################################
### West
###########################################################################
class West(Model):
    def __init__(self, IC_tag = 'fixed'):
        Model.__init__(self,
            IC_tag = IC_tag,
            model_key = 'west',
            model_function = (lambda param, time, tI, VI: von_bertalanffy_function(
                array([param[0], 3./4, param[1]]), time, tI, VI)),
            model_name = 'West',
            param0 = [1., 0.1],
            param_names = ['$a$','$b$']
            )
############################################################################################
############################################################################################
############################################################################################
############################################################################################
class ErrorModel:
    """
    Defines an error model for uncertainty on the data to be used in the fit.
    Can be either a function of the data or a vector containing the variance of the error at each data point
    """
    def __init__(self, is_vector = False, vectors = array([])):
        # Default error model taken from Benzekry et al., Plos Comp Biol 2014 (caliper measurements of longitudinal tumor growth)
        Vm               = 83
        alpha            = 0.84
        self.is_function = True
        self.function    = (lambda V: Vm**alpha*(V<Vm) + V**alpha*(V >= Vm))
        if is_vector:
            self.is_function = False
            self.is_vector   = is_vector
            self.vectors     = vectors
        self.sigma = 0.21
############################################################################################
############################################################################################
############################################################################################
############################################################################################
class LevenbergParameters:
    """Contains the required parameters for the Levenberg-Marquardt algorithm"""
    def __init__(self):
        #Value of the stopping threshold condition for the norm of the gradient
        self.accuracy       = 1E-6
        #To fasten up the calculation, put here a function with input the point,
        #and which returns the Jacobian matrix of the functions ri.
        #Otherwise, just keep it this way and it will calculate it through finite differences.
        self.jacobian       = []
        #Value of the lambda.
        self.lambda0        = 100.
        #Norm of the discretization step.
        self.discretization = 1E-4
        #Stopping condition on the number of iterations
        self.max_iter       = 400.
        #Stopping condition on the number of function evaluations
        self.max_eval       = 1500.
        #Value of nu in the LM-algorithm
        self.nu             = 10.
############################################################################################
############################################################################################

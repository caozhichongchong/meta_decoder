'''
Created on Apr 12, 2012

@author: jonathanfriedman

Find strain proportions and genomes from allele counts.
Uses exhaustive enumeration of all possible genome configurations for each site.
'''
from numpy import (ones, log, zeros, r_, c_, array, shape, load)
from openopt import NLP
from Cython.Distutils.extension import warnings
from util import likelihhodRatioTest, logMeanExp

import numpy as np
import os
try:
    from . import objective_fun_full
    c_avail = True
except ImportError:
    c_avail = False
    warnings.warn('Cannot import c-compiled objective function.' +
                      '\nSlower Python implementation will be employed.')

def objective(params, counts, alleles, marginalize, e):
    '''
    '''
    if e is None:
        e = params[0]
        p = params[1:]
    else:
        p = params
    ## possible allele probabilities
    k = counts.shape[1]
    A = alleles.shape[0]
    f = zeros((A,k))
    for i,row in enumerate(alleles): 
        for j,a in enumerate(row): f[i,a] += p[j] 
    
    fe  = f*(1-e)+(1-f)*e/(k-1)
    lfe = log(fe)
    lls = lfe.dot(counts.T)

    if marginalize:
        temp = logMeanExp(lls)
    else:
        temp = lls.max(axis=0)
    return temp.sum()

def fitNstrains(n, counts, marginalize=False, p0=None, e0=5e-3, e=None,
                tol=1e-4, mat_path=None, **kwargs):
    """
    Find the maximum-likelihood relative-abundances of n strains given the observed allele counts.
    
    Parameters
    ----------
    n : int 
        number of strains to fit
    counts : array-like
        Number of counts of each allele.
        Rows represent loci and columns represent alleles.
    marginalize : bool
        Flag for marginalizing the genomes' compositions.
        If True sum over all possible configurations.
        If False take only the most likely configuration for each site.
    p0 : array
        Initial guess for strain proportions.
        If None is passed will be set to a random sample for a uniform Dirichlet.
    e0 : float
        Initial guess for the magnitude of sequencing noise e.
        Only used if e is None.
    e : float
        Magnitude of sequencing noise e.
        If None is passes will be determined by the optimization procedure.
    tol : float
        Default value for optimization tolerance "ftol", "gtol", and "xtol".
    mat_path : str
        Path to location of allele permutation matrix.
        If None is passes uses the default path of supplied presence matrices.
    kwargs :
        Additional keyword arguments to be passed to OpenOpt's NLP object.
    
    Returns
    -------
    e : float
        Fitted magnitude of sequencing noise e.
    fracs : array
        Fitted strain relative abundances.
    ll : float
        Log-likelihood of best fit.
    sol : OpenOpt Solution 
        Object holding the full solution details. 
    """
    from numpy.random.mtrand import dirichlet
    
    ## load presence matrix
    k = shape(counts)[1]
    if mat_path is None:
        mat_path = os.path.dirname(os.path.realpath(__file__)) + '/presence_matrices'
        presence_file = mat_path + '/strains_%d.alleles_%d.npy' %(n,k)
        alleles = (load(presence_file)).astype(int)
    else:
        alleles = (load(mat_path)).astype(int)
    
    ## initial proportions guess
    if p0 is None:
        p0 = dirichlet([1]*n)
    params0 = p0
    if e is None:
        params0 = r_[e0,params0]       
    
    if marginalize or not c_avail:
        objective_fun = objective
        objective_args = (counts,alleles, marginalize, e)
    else:
        objective_fun = objective_fun_full.objective
        objective_args = (counts,alleles)
    
    ## tolerances
    kwargs.setdefault('ftol',tol)
    kwargs.setdefault('xtol',tol)
    kwargs.setdefault('gtol',tol) 
    iprint = kwargs.pop('iprint', -1)
    
    ## make problem object
    prob = NLP(x0=params0, f=objective_fun, **kwargs)
    
    ## box bounds
    lb = 0*1e-3*ones(n)
    ub = 1*ones(n)
    if e is None:
        lb = r_[1e-3, lb]
        ub = r_[0.01, ub]
    prob.lb = lb
    prob.ub = ub
    
    ## sum to 1 constraint
    if e is None:
        Aeq = r_[0, ones(n)]
    else:
        Aeq = ones(n)
    beq = 1
    prob.Aeq = Aeq
    prob.beq = beq
    
    ## param order constraints
    A = zeros((n-1,n))
    for i in range(n-1):
        A[i,i] = -1
        A[i,i+1] = 1
    if e is None:
        A = c_[zeros(n-1),A]
    b = zeros(n-1)
    prob.A = A
    prob.b = b
    
    ## additional arguments for objective
    prob.args.f = objective_args
    
    #solve problem
    sol   = prob.maximize('ralg', iprint=iprint)
    if e is None:
        e     = sol.xf[0]
        fracs = sol.xf[1:]
    else:
        fracs = sol.xf
    ll    = sol.ff
    return e, fracs, ll, sol 

def fitStrains(counts, nmax=10, marginalize=False, iters=10, max_iters=20, verbose=True,
               break_condition='AIC', break_th=0.05, tol=1e-4, mat_path=None, **kwargs):
    """
    Find the number of strains and their relative-abundances that best fit the observed allele counts.
    Starting from a single stain, additional strains are being sequentially added until the breaking condition is met.
    
    Parameters
    ----------
    counts : array-like
        Number of counts of each allele.
        Rows represent loci and columns represent alleles.
    nmax : int 
        Maximal number of strains to fit.
    marginalize : bool
        Flag for marginalizing the genomes' compositions.
        If True sum over all possible configurations.
        If False take only the most likely configuration for each site.
    iters : int  
        Number of optimization iterations ran for each number of strain.
        Helps escape local optima.
    max_iters : int
        Maximal number of optimization iterations ran for each number of strain.
        Used as some optimization iterations do not converge.
    verbose : bool
        Indicates whether to print fit results for each number of strains.
    break_condition : 'AIC'| 'LRT' | 'll'
        Stopping condition for adding more strains.
        'AIC' - Stop when Akaike's information criterion increases.
        'LRT' - Stop when the log-ratio test p-value > break_th
        'll'  - Stop when the relative log-likelihood < break_th
    break_th : float
        Threshold used to determine whether to add aditional strains.
        Used only is break_condition is 'LRT' or 'll'.
    tol : float
        Default value for optimization tolerance "ftol", "gtol", and "xtol".
    mat_path : str
        Path to location of allele permutation matrix.
        If None is passes uses the default path of supplied presence matrices.
    kwargs :
        Additional keyword arguments to be passed to OpenOpt's NLP object.
    
    Returns
    -------
    fracs : dict
        Fitted strain relative abundances for each number of strains.
    e : dict
        Fitted magnitude of sequencing noise e for each number of strains.
    ll : dict
        Log-likelihood of best fit for each number of strains.
    dof : dict
        Number of degrees-of-freedom used in the fit for each number of strains.
    """
    if break_condition not in ('AIC', 'LRT', 'll'):
        msg = 'Encountered unsupported break condition "{}". \nOptimization will not break until {} strains are fitted.'.format(break_condition, nmax)
        warnings.warn(msg)
    ns = np.arange(1,nmax+1)
    m  = len(counts)
    ll    = {}
    fracs = {}
    e     = {}
    dof   = {}
    AIC = {}
    if verbose:
        print('\t'.join(['# strains', 'll', 'AIC', 'strain fracs'])) 
    for n in ns:
        if marginalize:
            dof[n] = n # number of fit parameters
        else:
            dof[n] = n + n*m # number of fit parameters
            if e is not None: dof[n]-=1
        
        lln_best = -np.inf
        fracsn_best = None
        # set initial condition to previous best answer + 0 proportion
        if n == 1:
            p0 = array([1])
        else:
            p0 = r_[fracs[n-1],0] 
        ## run a few optimization iterations to find global optimum
        i = 0 # total number of iterations
        j = 0 # number of feasible iterations
        while i < max_iters and j < iters:
            if i == 0: # run with previous best strains as initial condition
                en, fracsn, lln, sol = fitNstrains(n,counts, marginalize=marginalize, tol=tol, 
                                                   mat_path=mat_path, p0=p0, **kwargs)
            else: # run with random initial condition
                en, fracsn, lln, sol = fitNstrains(n,counts,marginalize=marginalize, tol=tol, 
                                                   mat_path=mat_path, **kwargs)
            i+=1 # all iters
            if not sol.isFeasible:
                continue
            j+=1 # feasible iters            
            if lln > lln_best:
                lln_best    = lln
                fracsn_best = fracsn
                en_best     = en 
        ll[n]    = lln_best
        fracs[n] = fracsn_best
        e[n]     = en_best
        AIC[n] = 2*(dof[n]-ll[n])
        if verbose:
            print(n, ll[n], AIC[n], fracs[n])
        ## stop adding strains if no improvement
        if n !=1:
            if break_condition=='AIC':
                if AIC[n] > AIC[n-1]: break
            elif break_condition=='LRT':
                p_val = likelihhodRatioTest(ll[n-1],ll[n],dof[n-1],dof[n])
                if p_val > break_th: break
            elif break_condition=='ll':
                ll_improve = np.abs((ll[n] - ll[n-1])/ll[n-1])
                if ll_improve <= break_th: break
    return fracs, e, ll, dof 

def genomes_given_fracs(counts, fracs, e, mat_path=None, alleles=None, return_all=False):
    """
    Find the maximum-likelihood genomes given the strains' relative-abundances
    and the allele counts matrix.
    
    .. note::
        Each position is the genomes is inferred only based on the allele counts in that position.
        (Unlike the strain fraction, which are constrained by all positions.) 
        Therefore, they can be inferred as accurately as the strain proportions. 
    
    Parameters
    ----------
    counts : array-like
        Number of counts of each allele.
        Rows represent loci and columns represent alleles.
    fracs : array
        Strain relative abundances.
    e : float
        Magnitude of sequencing noise e.
    mat_path : str
        Path to location of allele permutation matrix.
        If None is passes uses the default path of supplied presence matrices.
    alleles : iterable
        List of characters representing the different alleles.
        Order corresponds to order of counts array's columns.
        If None is passed integers will be used. 
    return_all : bool
        indicated whether to return the ll of all alleles for each strain in each position, or just the ML genomes.
    
    Returns
    -------
    genomes : array
        Inferred genomes.
        Rows represent loci and columns represent strains.
    """
    k = counts.shape[1]
    n = len(fracs)
    if mat_path is None:
        mat_path = os.path.dirname(os.path.realpath(__file__)) + '/presence_matrices'
        presence_file = mat_path + '/strains_%d.alleles_%d.npy' %(n,k)
        allele_perm = (load(presence_file)).astype(int)
    else:
        allele_perm = (load(mat_path)).astype(int)
    A = allele_perm.shape[0]
    f = zeros((A,k))
    for i,row in enumerate(allele_perm): 
        for j,a in enumerate(row): f[i,a] += fracs[j] 
    
    fe  = f*(1-e)+(1-f)*e/(k-1)
    lfe = log(fe)
    lls = lfe.dot(counts.T)
    
    if return_all:
        return lls
    else:
        llinds = lls.argmax(axis=0)
        genomes = allele_perm[llinds]
        if alleles is not None:
            if len(alleles)==k:
                genomes = genomes.astype(object)
                for i,a in enumerate(alleles):
                    genomes[genomes==i] = a
            else: 
                warnings.warn('Numer of passed alleles does not match shape of counts matrix. Ignoring passed alleles.')
        return genomes        

if __name__ == '__main__':   
    pass

    
    

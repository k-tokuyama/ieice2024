
import sys
import csv
import numpy as np
import math
import time

import sobol_seq

from scipy import integrate
#from scipy.integrate import quad
#from scipy.special import i0
#from scipy.stats import rice





def NoticeMessages():
    print('<< Notice >>')
    print('The list argument "param_sets" must contain the following 10 elements in its list.')
    print('the 10 elements -> 1 : skipping time for ppp')
    print('                   2 : skipping time for tcp')
    print('                   3 : moving velocity of a user')
    print('                   4 : intensity for ppp')
    print('                   5 : parent intensity for tcp')
    print('                   6 : daughter intensity for tcp')
    print('                   7 : daighter variance for tcp')
    print('                   8 : transmitting power for ppp (macro base stations)')
    print('                   9 : transmitting power for tcp (small base stations)')
    print('                   10: path-loss exponent')

def MeanPeriod(param_sets, flag_sobol, McNum_sobol=1024, DivNum=100):

    ### parameter input
    #s1, s2, veloc, lam1, lam2p, mb, sigma, P0, P1, beta = param_sets
    s1, s2, veloc, lam1, lam2p, mb, rd, P0, P1, beta = param_sets

    def output_AProb_pcp():
        
        ## random number generater for montecalro integration.
        ##   - lower bounds of integration variables must be 0
        ##   - uppber bounds are callable with the list: UB_l
        def rd_gen(UB_l, samplenum, flag_sobol=False):    ## UB_l: list of upper bounds of integration variable
            if not type(samplenum) == type(1):
                print('Error(rd_gen): argument samplenum must be Integer type.')
                sys.exit()
            if not type(UB_l) == type([1.0, 2.0, 3.0]):
                print('Error(rd_gen): argument UB_l must be List type.')
                sys.exit()
            varnum = len(UB_l)    ## number of integration variables
            if flag_sobol:
                return sobol_seq.i4_sobol_generate(varnum, samplenum) * np.array(UB_l)
            else:
                #np.random.seed(0)
                return np.random.uniform(0, 1, varnum*samplenum).reshape(samplenum, varnum) * np.array(UB_l)

        
        def P_bar(tier_i, tier_j):
            def P(tier_i):
                if tier_i == 1:
                    return P0
                elif tier_i == 2:
                    return P1
                else:
                    print('Error;def P(tier_i): tier_i must be either 1 or 2.')
                    sys.exit()
            return (P(tier_j)/P(tier_i))**(1/beta)

        #def func_fd_M(x, z):
        #    return rice.pdf(x, z/sigma, scale=sigma)
        #def func_Fd_M(r, z): 
        #    return rice.cdf(r, z/sigma, scale=sigma)
        def func_fd_M(r, z):
            if r<=max(rd - z, 0):
                return 2*r/rd**2
            elif abs(rd - z)<=r and r<=rd + z: 
                return 1/np.pi*np.arccos( round((r**2 + z**2 - rd**2)/(2*r*z), 8) )*2*r/rd**2
            else:
                return 0.0
        def func_Fd_M(r, z):
            def integd_x(L_x, z):
                return np.array([ x*np.arccos( round((x**2 + z**2 - rd**2)/(2*x*z), 8) ) for x in L_x])

            if min(r, abs(rd - z)) == min(r, rd + z):
                return (min(r, max(rd - z, 0))**2)/rd**2
            else:
                L_x = np.linspace(min(r, abs(rd - z)), min(r, rd + z), DivNum+1)
                return (min(r, max(rd - z, 0))**2 + 2/np.pi*integrate.trapz(integd_x(L_x, z), L_x))/rd**2
        
        def I00(r, z):
            return np.exp( -mb * func_Fd_M(P_bar(1, 1)*r, z) )
        def I01(r, z):
            return np.exp( -mb * func_Fd_M(P_bar(2, 1)*r, z) )
        def I10(r, z):
            return np.exp( -mb * func_Fd_M(P_bar(1, 2)*r, z) )
        def I11(r, z):
            return np.exp( -mb * func_Fd_M(P_bar(2, 2)*r, z) )


        def integ__H1(r):

            def integ__H1_1(r):
                def trapz__H1_1(r):
                    def func__H1_1(r, w):
                        return ( func_fd_M(r, w) * I11(r, w)  )*w*(1+w**2)
                    ### Trapezoidal integral over [0, 2*np.pi]
                    ## w = tan(psi)
                    w_UB = r + 0.5; w_LB = np.max([0.0, r - 0.5]);    ## the margin 0.5 could be managed depending on sigma: the daughter variance.
                    psi_range = np.linspace(np.arctan(w_LB), np.arctan(w_UB), DivNum+1)
                    return np.nanmean([func__H1_1(r, np.tan(psi)) for psi in psi_range])*(np.arctan(w_UB) - np.arctan(w_LB))
                return np.exp(-np.pi*lam1*P_bar(2, 1)**2*r**2) * trapz__H1_1(r)
            def integ__H1_2(r):
                def func__H1_2(r, w):
                    return (1 - I11(r, w))*w*(1+w**2)
                ### Trapezoidal integral over [0, 2*np.pi]
                ## w = tan(psi)
                psi_range = np.linspace(0, np.pi/2, DivNum+1)
                trapz__H1_2 = np.nanmean([func__H1_2(r, np.tan(psi)) for psi in psi_range])*np.pi/2
                return np.exp( -2*np.pi*lam2p*trapz__H1_2 )

            return integ__H1_1(r) * integ__H1_2(r) * (1+r**2)
        


        ### montecalro integration
        UB_phi, UB_vphi, UB_theta, UB2_theta = np.pi/2, np.pi/2, 2*np.pi, np.pi
        UB_phi_alt = np.pi/4


        mcint__H1 = UB_phi_alt * np.nanmean([integ__H1(np.tan(phi_list[0])) for phi_list in rd_gen([UB_phi_alt], McNum_sobol, flag_sobol)])

        return 2*np.pi*mb*lam2p*mcint__H1

    
    #start = time.time()
    AProb_pcp = output_AProb_pcp()
    #finish = time.time()

    AProb_ppp = 1.0 - AProb_pcp

    ### test print ###
    #print('AProb_pcp(small) = {}, elapsed time = {}'.format(AProb_pcp, finish - start))

    return s1*AProb_ppp + s2*AProb_pcp






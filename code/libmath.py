#!/usr/bin/env python
#This script contains functions to evaluate <\chi_m(x)|\chi_n(x)> and <\chi_m(x)|x|\chi_n(x)>
#Warning: only 
from math import exp, sqrt, pi
import numpy as np
import numpy.polynomial.polynomial
import scipy.integrate
import mpmath
from mpmath import mp, mpf
mp.dps = 50
#mp.prec = 80
mp.pretty = False

#Enable mpmath?
#If not using mpmath the summation of Gaussian integrals are not stable as there are lots of very large numbers and small numbers
b_mp = True

#Note we must apply mpf to all input variables like frequencies and bionmial coefficients

#if (b_mp):
#    cache_fac = [mpf(1)]
#    cache_fac2 = [mpf(1),mpf(1)]
#else:
#    cache_fac = [1]
#    cache_fac2 = [1,1]
cache_fac = [1]
cache_fac2 = [1,1]
cache_bi = {}

def factorial(n):
    '''
    Factorial
    '''
    if (len(cache_fac) < n+1):
        n0 = len(cache_fac)
        if (not b_mp):
            for i in range(n0, n+1):
                cache_fac.append(cache_fac[-1] * i)
        else:
            for i in range(n0, n+1):
                cache_fac.append(cache_fac[-1] * mpf(i))
            

    if (n < -1):
        raise ValueError("n < 0")
    elif (n == -1):
        return 1
    return cache_fac[n]

def fac2(n):
    '''
    Double factorial
    '''
    if (len(cache_fac2) < n+1):
        n0 = len(cache_fac2)
        for i in range(n0, n+1):
            cache_fac2.append(cache_fac2[-2] * i)

    if (n < -1):
        raise ValueError("n < -1")
    elif (n == -1):
        return 1
    return cache_fac2[n]

def binomial(n,m):
    '''
    Binomial coef
    '''
    if (not ( m >= 0 and n >= m)):
        raise ValueError("n / m not valid")
    if ((n,m) not in cache_bi):
        cache_bi[(n,m)] = factorial(n) / factorial(m) / factorial(n - m)
#       print([factorial(n), factorial(m), factorial(n-m), cache_bi[(n,m)]])

    return cache_bi[(n,m)]

def Rational(n,m):
    '''
    rational number (just use float to replace)
    '''
    if (not b_mp): #This gives large error
        return 1.0 * n / m
    else:
        return mpf(n) / mpf(m)
#   return 1.0 * n / m
    


def create_hermite_polynomial_1_to_n(n):
    '''
    Create an array, each includes the coefficient of Hermite polynomial from order 0 to n
    Coefficients in order of 1, x, x^2, ...
    '''
    if (b_mp):
        l0 = [mpf(1)]
    else:
        l0 = [1]
    list_coef = [l0]
    for i in range(1, n+1):
        l1 = [0 for x in range(i+1)]
        for j in range(i):
            if (j > 0):
                l1[j-1] += j * l0[j]
            l1[j+1] += -2 * l0[j]
        l1 = [x * -1 for x in l1]
        list_coef.append(l1)
        l0 = l1
    
    return list_coef

cache_hermite = []

def get_hermite_polynomial(n):
    '''
    Get the Hermite polynomial coefficients in 1,x,x^2 ... order
    '''
    if (n < 0):
        raise ValueError("Order must be 0 or positive")
    global cache_hermite
    if (len(cache_hermite) < n+1):
        cache_hermite = create_hermite_polynomial_1_to_n(n)

    return cache_hermite[n]

def eval_hermite_polynomial(n, x):
    '''
    Get the value of n-th Hermite polynomial at x
    '''
    y = numpy.polynomial.polynomial.polyval(x, get_hermite_polynomial(n))
    return y



#These functions used to calculate <f1|x|f2>
#Where f1 and f2 are terms in quantum harmonic oscillator eigenstates at different equalibrium geometry
#(sqrt(m\omega)*(-R+x))^n exp(-(x-R)^2/2*sqrt(m\omega)) * sqrt(sqrt(m\omega / pi) / 2^n / n!)
#For R1 and R2
#
#By 2Sa = R1 + R2, 2Sb = R1 - R2
#For m\omega = 1 on both side, this can be redunced to 
#(-Sb+y)^n(Sb+y)^m(Sa+y)e^(-Sb^2)e^(-y^2)  / sqrt(pi * 2^(m+n) n! m!)
#Where y = x - Sa


#For m\omega is different for two sides, this can be reduced to
#x = (z+ Sa - Lb / La * Sb)
#2La = m1\omega1 + m2\omega2
#2Lb = m1\omega1 - m2\omega2
#
#(La *z + m1\omega_1 Sb)^m * (La*z - m2\omega_2 Sb)^n * exp(-La * z^2) * constant
#Constant = sqrt(1/(2^(m+n) * pi * m! * n!)) * La^(-m-n) * (m1\omega_1)^(n/2+1/4) * (m2\omega_2)^(m/2+1/4) * exp(-m1\omega1m2\omega2 Sb^2/La)

def calc_gaussian_x_alpha(n, m):
    '''
    Calculate the Gaussian integral from -infty to infty 
    Equations see above
    Return coefficients, meaning of coefficients see below
    '''
    #Each array 
    #Index is order of z
    #Store: a list of [coef, order of Sb, order of La, order of M2w2, order of M1w1], all for z^i
    list_order_Sa = [[] for x in range(n+m+1)]
    list_order_y = [[] for x in range(n+m+1+1)]
    
    #Order of y
    for ixn in range(0, n+1):
        for ixm in range(0, m+1):
            coef = binomial(n,ixn) * binomial(m, ixm)*(-1)**(n-ixn) #* (-sb)**(n-ixn)*sb**(m-ixm)
            t = [coef, n+m-ixm-ixn, ixm+ixn,n-ixn, m-ixm]
            if ((ixn + ixm) % 2 == 0): #Even ,for constant term in the *x part
                list_order_Sa[ixn + ixm].append(t)
            else: #Odd, for x term in *x part
                list_order_y[ixn + ixm + 1].append(t)
#Convert coefficient to integrals (except sqrt(pi/La) and La^-k part)
#Note here c0 is square of the actual coefficient
    c0 = 2**(m+n)*factorial(m)*factorial(n)

    #Coefficient is converted into 
    for l0 in (list_order_Sa, list_order_y):
        for i in range(len(l0)):
            k = i // 2
            for c in l0[i]:
                #print(c)
                if (c[0] != 0):
                    c[0] = (c[0] ,fac2(2*k-1)**2, 2**(k*2) * c0)
                else:
                    c[0] = (0, 0, 0)
        
#   for x in list_order_y:
#       print(x)
    return list_order_Sa, list_order_y

cache_gauss_x1 = {}
n_cache_gauss_miss = [0,0]

def calc_gaussian_x_alpha_coef(n, m):
    '''
    Wrapper of calc_gaussian_x_alpha, 
    Cached calculated values
    '''
    n_cache_gauss_miss[1] += 1
    if ((n,m) not in cache_gauss_x1):
        cache_gauss_x1[(n,m)] = calc_gaussian_x_alpha(n, m)
        n_cache_gauss_miss[0] += 1

    return cache_gauss_x1[(n,m)]

def report_cache_gauss_miss():
    '''
    Report percentage of miss
    '''
    if (n_cache_gauss_miss[1] > 0):
        print("Gaussian integral coefficient miss %i/%i (%.6f)" % (n_cache_gauss_miss[0], n_cache_gauss_miss[1], 
            n_cache_gauss_miss[0] * 1.0 / n_cache_gauss_miss[1]))


def calc_gaussian_x_alpha_val_part(n, m, R1, R2, momega1, momega2, order_x=1):
    '''
    Calculate the number of Gaussian integral (with or without  x) of one term in Hermite polynomial 

    except the constant exp part
    except constant 1/sqrt(2^n_level*n_level!) where n_level is quantum number (not order!)
    with additional 1/sqrt(2^n*n!) where n is order

    R, momega : position and m*\omega for the oscilator
    n : order of x^n in the Hermite polynomial of oscillator 1
    m : of oscillator 2
    '''
    list_order_Sa, list_order_y = calc_gaussian_x_alpha_coef(n, m)
    Sa = (R1 + R2) / 2.0
    Sb = (R1 - R2) / 2.0
    Lmomegaa = (momega1 + momega2) / 2.0
    Lmomegab = (momega1 - momega2) / 2.0
    if (b_mp):
        Sa = mpf(Sa)
        Sb = mpf(Sb)
        Lmomegab = mpf(Lmomegab)
        Lmomegaa = mpf(Lmomegaa)

    #print(list_order_y)
    t = 0.0

#Convert to float
    n = float(n)
    m = float(m)

    dic_list_coef = {
            0 : [(1, list_order_Sa),],
            1 : [((Sa + Sb * Lmomegab / Lmomegaa), list_order_Sa),
                (1, list_order_y),],
            }
    
    for term0, list_order in dic_list_coef[order_x]:
        for i, list_o in enumerate(list_order):
            if (i % 2 != 0):
                continue
            k = i // 2
            for (c, num, den), oSb, oLa, o2, o1 in list_o:
                if (num != 0):
#                   print("z La Sb m2w2 m1w1 coef num den", i, oLa, oSb, o2, o1, c, num, den)
                    if (b_mp):
                        t1 = mpmath.power(Lmomegaa, -k-Rational(1,2))  \
                            * term0 \
                            *mpmath.power(Sb,oSb) \
                            *mpmath.power(Lmomegaa, oLa) \
                            *mpmath.power(momega1, o1) \
                            *mpmath.power(momega2, o2) \
                            * c * Rational(int(num), int(den)) **Rational(1,2)

                    else:
                        t1 = Lmomegaa ** (-k-Rational(1,2))  \
                            * term0 \
                            *Sb**oSb*Lmomegaa**oLa* momega1**o1 * momega2**o2 *\
                            c* Rational(int(num), int(den)) **Rational(1,2)
#                   print([t1])
                    t += t1

    if (not b_mp):
        t = t * Lmomegaa**(-m-n) * momega1 ** Rational(2*n+1,4) * momega2 ** Rational(2*m+1,4)
    else:
        t = t * mpmath.power(Lmomegaa, -m-n) * mpmath.power(momega1, Rational(2*n+1,4)) \
                * mpmath.power(momega2, Rational(2*m+1,4))


#   print([t])
    return t


def calc_gaussian_x_alpha_val(n, m, R1, R2, momega1, momega2, order_x=1):
    '''
    Calculate the number of Gaussian integral (with x) of one term in Hermite polynomial (include the constant exp part)
    except constant 1/sqrt(2^n_level*n_level!) where n_level is quantum number (not order!)
    with additional 1/sqrt(2^n*n!) where n is order
    '''
    t = calc_gaussian_x_alpha_val_part(n, m, R1, R2, momega1, momega2, order_x)
    Sb = (R1 - R2) / 2.0
    if (b_mp):
        Sb = mpf(Sb)
    t = t * exp(-2*momega1*momega2*Sb**2/(momega1 + momega2))
    return t

def overlap_x_gaussian(ni, nf, oi, of, Ri, mi, freqi, mf, freqf, order_x):
    '''
    Calculate <\chi_fm|Q-Q0|\chi_in> where Q0 is equalibrium geomeotry of chi_f, only one Gaussian term in Hermite equation
    Wrapper of calc_gaussian_x_alpha_val for two quantum oscialltor


    Note :  in calc_gaussian_x_alpha_val:
    except constant 1/sqrt(2^n_level*n_level!) where n_level is quantum number (not order!)
    with additional 1/sqrt(2^n*n!) where n is order
    '''
    if (b_mp):
        Ri = mpf(Ri)
        freqi = mpf(freqi)
        freqf = mpf(freqf)
    R1 = Ri
    R2 = 0
    momega1 = mi * freqi
    momega2 = mf * freqf
    t = calc_gaussian_x_alpha_val(oi, of, R1, R2, momega1, momega2, order_x)
    t *= 2**(-(ni+nf-oi-of)/2.0) * sqrt(1.0 * factorial(oi)/factorial(ni)*factorial(of)/factorial(nf))

    return t

def overlap_x_quantum_harmonic(ni, nf, Ri, mi, freqi, mf, freqf, order_x=1, debug_print=False):
    '''
    Calculate <\chi_fm|Q-Q0|\chi_in> where Q0 is equalibrium geomeotry of chi_f
    This is equavlent to <\chi_fm|Q|\chi_in> with origin shifted so chi_f is located at zero

    Or <\chi_fm|\chi_in> overlap integral

    Note :  in calc_gaussian_x_alpha_val:
    except constant 1/sqrt(2^n_level*n_level!) where n_level is quantum number (not order!)
    with additional 1/sqrt(2^n*n!) where n is order

    :param order_x: 1 to calculate <\chi|Q|\chi>, 0 to calculate <\chi|\chi>
    '''
    if (b_mp):
        Ri = mpf(Ri)
        freqi = mpf(freqi)
        freqf = mpf(freqf)
    R1 = Ri
    R2 = 0
    momega1 = mi * freqi
    momega2 = mf * freqf
        
    list_ci = get_hermite_polynomial(ni)
    list_cf = get_hermite_polynomial(nf)
    t = 0
    list_sum = []
    for oi, ci in enumerate(list_ci):
        for of, cf in enumerate(list_cf):
            if (ci == 0 or cf == 0):
                continue
#Integral and coef
            t0 = calc_gaussian_x_alpha_val(oi, of, R1, R2, momega1, momega2, order_x)
            if (b_mp):
                t1 = mpmath.power(2, -mpf(ni+nf-oi-of)/2.0) * mpmath.sqrt(mpf(1.0) * factorial(oi)/factorial(ni)*factorial(of)/factorial(nf))
            else:
                t1 = 2**(-(ni+nf-oi-of)/2.0) * sqrt(1.0 * factorial(oi)/factorial(ni)*factorial(of)/factorial(nf))
#MPMATH version
#           t += mpf(t0) * mpf(t1) * mpf(ci) * mpf(cf)
#Normal version
            t += t0 * t1 * ci * cf
#           print([t0,t1,ci,cf])
            list_sum.append(t0 * t1 * ci * cf)
            if (debug_print):
                print("%3i %3i %16.7g %16.7g %16.7g %16.7g %25.18g"  % (oi, of, ci, cf, t0, t0*t1, ci*cf*t1*t0))
#Another summation to avoid numerical errors ; not that useful
#   list_sum.sort(key=lambda x:abs(x), reverse=True)
#   t2 = sum(list_sum)
#   print("Summation different: ", t,t2)
            
    return t


def overlap_x_quantum_harmonic_ladder(ni, nf, Ri, mi, freqi, mf, freqf, order_x=1, debug_print=False):
    '''
    Same as overlap_x_quantum_harmonic, but always use ladder operator to convert to two overlap integral
    This is much more stable than directly integrate <chi|x|chi>
    '''
    if (order_x == 0 or nf == 0):
        return overlap_x_quantum_harmonic(ni, nf, Ri, mi, freqi, mf, freqf, order_x, debug_print)
    
    if (order_x == 1):
        ta = overlap_x_quantum_harmonic(ni, nf-1, Ri, mi, freqi, mf, freqf, order_x=0, debug_print=debug_print)
        tc = overlap_x_quantum_harmonic(ni, nf+1, Ri, mi, freqi, mf, freqf, order_x=0, debug_print=debug_print)

        mfreqf = mf * freqf
        ca = 1 / sqrt(mfreqf * 2) * sqrt(nf)
        cc = 1 / sqrt(mfreqf * 2) * sqrt(nf+1)
        t = ca * ta + cc * tc
        return t

def overlap_x_quantum_harmonic_ladder_hr(nf, Ri, mf, freqf, order_x=1, debug_print=False):
    '''
    Using Huang-Rhys factor to compute the <0|n> integral and <0|x|n> integral
    initial/final m and frequency are approximate to be same to applied equation 
    |<0|n>|^2 =e^-S S^n / n!
    '''
    S=  freqf * Ri ** 2 / 2
    if (order_x == 0):
        return exp(-S/2) * S**(nf/2.0) / sqrt(factorial(nf))
    
    if (order_x == 1):
        if (nf >= 1):
            ta = overlap_x_quantum_harmonic_ladder_hr(nf-1, Ri, mf, freqf, order_x=0, debug_print=debug_print)
        else:
            ta = 0
        tc = overlap_x_quantum_harmonic_ladder_hr(nf+1, Ri, mf, freqf, order_x=0, debug_print=debug_print)

        mfreqf = mf * freqf
        ca = 1 / sqrt(mfreqf * 2) * sqrt(nf)
        cc = 1 / sqrt(mfreqf * 2) * sqrt(nf+1)
        t = ca * ta + cc * tc
        return t


def eval_gaussian(x, n_level, n_order, R, m, freq):
    '''
    Eval x^n exp(-x^2) style term in eigenstates
    '''
    x = x - R
    mfreq = m * freq
    y = 1/ sqrt(1.0 * factorial(n_level) * 2**n_level) * (mfreq / pi) **0.25 * exp(-mfreq*x**2/2) * (x * sqrt(mfreq))**n_order
#   print(1/ sqrt(1.0 * factorial(n_level) * 2**n_level))
#   y = (mfreq / pi) **0.25 * exp(-mfreq*x**2/2) * (x * sqrt(mfreq))**n_order
#   print(y)
    return y

def eval_overlap_x_gaussian(x, ni, nf, oi, of, Ri, mi, freqi, mf, freqf, order_x):
    '''
    Eval <x^n exp(-x^2)|x|x^m exp(-x^2)> style term in harmonic eigenstates value at given x
    '''
    y =  eval_gaussian(x, nf, of, 0, mf, freqf) * x ** order_x * eval_gaussian(x, ni, oi, Ri, mi, freqi)
#   print(x, y)
    return y

def calc_gaussian_x_alpha_val_num(n, m, Ri, momega1, momega2, order_x=1):
    '''
    '''
    y, abserr = scipy.integrate.quad(eval_overlap_x_gaussian, -np.inf, np.inf, args=(n, m, n, m, Ri, 1, momega1, 1, momega2, order_x),limit=500, epsrel=1e-14, epsabs=1e-14)
#   print(y, abserr)
    return y

def overlap_x_gaussian_num(ni, nf, oi, of, Ri, mi, freqi, mf, freqf, order_x):
    y, abserr = scipy.integrate.quad(eval_overlap_x_gaussian, -np.inf, np.inf, args=(ni, nf, oi, of, Ri, mi, freqi, mf, freqf, order_x),limit=500, epsrel=1e-14, epsabs=1e-14)
#   print(y, abserr)
    return y


def eval_quantum_harmonic(x, n, R, m, freq):
    '''
    Evaluate quantum harmonic eigenstate at x with n/R/m/freq
    '''
    x = x - R
    y = 1/ sqrt(factorial(n) * 2**n) * (m * freq / pi) **0.25 * exp(-m*freq*x**2/2) * eval_hermite_polynomial(n, x * sqrt(m*freq))
    return y

def eval_overlap_x(x, ni, nf, Ri, mi, freqi, mf, freqf, order_x):
    '''
    Evaluete  chi_fm(x) * x  * chi_in(x)
    '''
    y =  eval_quantum_harmonic(x, nf, 0, mf, freqf) * x**order_x * eval_quantum_harmonic(x, ni, Ri, mi, freqi)
#   print(x, y)
    return y

def overlap_x_quantum_harmonic_num(ni, nf, Ri, mi, freqi, mf, freqf, order_x=1):
    '''
    Numerical integration of <\chi_jm|Q-Q0|\chi_in> where Q0 is equalibrium geomeotry of chi_i
    This is equavlent to <\chi_jm|Q|\chi_in> with origin shifted so chi_i is located at zero
    '''
    y, abserr = scipy.integrate.quad(eval_overlap_x, -np.inf, np.inf, args=(ni, nf, Ri, mi, freqi, mf, freqf, order_x),limit=500)
#   print(y, abserr)
    return y

def overlap_x_quantum_harmonic_num2(ni, nf, Ri, mi, freqi, mf, freqf, order_x=1, debug_print=False):
    '''
    Same as overlap_x_quantum_harmonic_num, but summation after integration
    '''
    R1 = Ri
    R2 = 0
    momega1 = mi * freqi
    momega2 = mf * freqf
    list_ci = get_hermite_polynomial(ni)
    list_cf = get_hermite_polynomial(nf)
    t = 0
    for oi, ci in enumerate(list_ci):
        for of, cf in enumerate(list_cf):
            if (ci == 0 or cf == 0):
                continue
            t0 = calc_gaussian_x_alpha_val_num(oi, of, Ri, momega1, momega2, order_x)
            t1 = 2**(-(ni+nf-oi-of)/2.0) * sqrt(1.0 * factorial(oi)/factorial(ni)*factorial(of)/factorial(nf))
            t += t1 * ci * cf * t0
            if (debug_print):
                print("%3i %3i %16.7g %16.7g %16.7g %16.7g %16.7g"  % (oi, of, ci, cf, t0, t0*t1, ci*cf*t1*t0))
            
    return t


######################
# Python 3
# 
# takes tight-binding parameters of the 2D lattice (including mu)
# takes local self-energy given in two files for real and imaginary parts 
# computes local Green's function on the lattice
# prints out hybridization function for the next DMFT iterations in two files, re and im parts separately
#
# usage: 
# $ python get_Delta.py [-param=value .... or just -param if param is supposed to be true], 

import numpy
import scipy
import itertools
import sys
import os
import time
from numpy import cos, sin
from numba import njit
import cubepy as cp

print("Welcome to get_Delta_serial.py!!!!")

   

#@njit
#def optimized_inv(M):
#    #returns the inverse of M, but only 0,0 element is actually computed
#    #we assume real symmetric matrix, and only use M_ij with i<=j
#    #res = numpy.zeros((3,3),dtype=numpy.complex_)
#    det = M[0,0]*M[1,1]*M[2,2] + 2*M[0,1]*M[1,2]*M[0,2]\
#        - M[0,0]*M[1,2]**2 - M[0,1]**2*M[2,2] - M[0,2]**2*M[1,1]
#        #- M[0,0]*M[1,2]*M[2,1] - M[0,1]*M[1,0]*M[2,2] - M[0,2]*M[1,1]*M[2,0]
#    #res[0,0] = (M[1,1]*M[2,2]-M[1,2]**2)/det
#    return (M[1,1]*M[2,2]-M[1,2]**2)/det

def write_output(fnr, fni, xs, ys):
    #take function y given on a real grid x and printout to two separate files real and imaginary part of y
    r_data = numpy.array(list(zip(xs,ys.real)))
    i_data = numpy.array(list(zip(xs,ys.imag)))

    numpy.savetxt(fnr, r_data)
    numpy.savetxt(fni, i_data)

def load_input(fnr, fni):
    #load complex function y of real argument x written in two separate files for real and imaginary part of y
    i_data = numpy.loadtxt(fni)
    r_data = numpy.loadtxt(fnr)
    assert numpy.all(i_data[:,0]==r_data[:,0]), "real and imaginary part of Sigma must be given on the identical frequency grid"
    xs = i_data[:,0]
    ys = r_data[:,1] + 1j*i_data[:,1]
    return xs, ys


# default parameters
params = {
  'fnis': 'res/imsigma.dat', #name of the input file for the imaginary part of Sigma
  'fnrs': 'res/resigma.dat', #name of the input file for the real part of Sigma
  'fnig': 'res/ImGlattloc.dat', #name of the output file for the imaginary part of Delta
  'fnrg': 'res/ReGlattloc.dat', #name of the output file for the real part of Delta
  'fnid': 'ImDelta.dat', #name of the output file for the imaginary part of Delta
  'fnrd': 'ReDelta.dat', #name of the output file for the real part of Delta
  'model': 'Emery', #'Hubbard' 
  'fnparam': 'tb_params.py', #file where the tight-binding parameters are stored, should be relevant for the model used
  'check': False, #in the commandline args use 0/1 for False/True
  'test': False, #if True, printout a Sigma to test the code,
#  'deal_with_Delta_shift': False,
  'use_Hubbard_DOS': True,
  'fndos': 'DOS.dat'  
}

# read the commandline arguments to set parameters
for s in sys.argv:
    if s[0]!='-': continue
    if "=" not in s:
        key, val = s[1:], True
    else:
        key, val = s[1:].split("=")    
    if type(params[key])==bool:
        params[key] = bool(int(val))
    elif type(params[key])==int:
        params[key] = int(val)
    else: params[key] = val

print("params:",params)

# if test mode - printout a simple self-energy Sigma(w) = i Gamma, to be used as input
if params['test']:
    print("printing out test Sigma", end="...")
    test_ws = numpy.linspace(-5,5,1000,endpoint=True)    
    test_Sigma = -0.05j*numpy.ones((len(test_ws),))
    write_output(params['fnrs'], params['fnis'], test_ws, test_Sigma)
    print("done")

# take the tb parameters from file
print("taking parameters from ",params['fnparam'])
exec(open(params['fnparam']).read())
print("done")

# load the self-energy
print("about to load Sigma from ",params['fnrs'], params['fnis'],end="...")
ws, Sigmaw = load_input(params['fnrs'], params['fnis'])
nw = len(ws)
print("done")

def get_Glattloc_Emery(nw, ws, Sigmaw, mu, epsd, epsp, tpd, tpp, tppp):
    # get the Glatt loc
    Glattloc = numpy.zeros(len(ws),dtype=numpy.complex_)
    for wi,(w,Sigma) in enumerate(zip(ws,Sigmaw)):
        #if wi!=500: continue 
        print("w=",w, " Sigma=", Sigma)

        low = [0.0, 0.0]
        high = [numpy.pi, numpy.pi]

        if True: #True is a bit faster, but False also works
            #precomputations
            w_plus_mu = w+mu
            minus_epsp_plus_w_plus_mu = -epsp + w_plus_mu
            minus_epsd_plus_w_plus_mu_minus_Sigma = -epsd + w_plus_mu - Sigma
            tpdsquarred = tpd**2

            def comp(kx, ky):
                coskx = numpy.cos(kx)
                cosky = numpy.cos(ky)
                sinkxhalfsquarred = numpy.sin(kx/2)**2
                sinkyhalfsquarred = numpy.sin(ky/2)**2
                two_tppp_coskx = 2 * tppp * coskx
                two_tppp_cosky = 2 * tppp * cosky
                
                numerator = ((minus_epsp_plus_w_plus_mu - two_tppp_coskx) * (minus_epsp_plus_w_plus_mu - two_tppp_cosky) -
                  16 * tpp**2 * sinkxhalfsquarred * sinkyhalfsquarred)
                denominator = ((2 * tpdsquarred * (-1 + coskx) + minus_epsd_plus_w_plus_mu_minus_Sigma  * (minus_epsp_plus_w_plus_mu - two_tppp_coskx)) *
                  (minus_epsp_plus_w_plus_mu - two_tppp_cosky) +
                  4 * tpdsquarred * (epsp - w_plus_mu - 2 * tpp + 2 * (tpp + tppp) * coskx) * sinkyhalfsquarred +
                  16 * tpp * (-tpdsquarred + tpp * (-minus_epsd_plus_w_plus_mu_minus_Sigma )) * sinkxhalfsquarred * sinkyhalfsquarred)
                #res = numerator / denominator
                #print("res.shape:",res.shape, "kx.shape:",kx.shape)
                return numerator / denominator#res
        else:
            w_plus_mu = w+mu
            w_plus_mu_minus_epsp = w_plus_mu-epsp 
            minus_two_tpd = -2*tpd
            two_tppp = 2*tppp
            M00 = w_plus_mu-epsd-Sigma
            def comp(kx, ky):
                coskx = numpy.cos(kx)
                cosky = numpy.cos(ky)
                sinkxhalf = numpy.sin(kx/2)
                sinkyhalf = numpy.sin(ky/2)                      
                M01 = minus_two_tpd*sinkxhalf
                M02 = -minus_two_tpd*sinkyhalf
                M11 = w_plus_mu_minus_epsp - two_tppp*coskx
                M12 = 4*tpp*sinkxhalf*sinkyhalf
                M12squared = M12**2   
                M22 = w_plus_mu_minus_epsp - two_tppp*cosky
                det = M00*M11*M22 + 2*M01*M12*M02 - M00*M12squared - M01**2*M22 - M02**2*M11
                return (M11*M22-M12squared)/det

        def f1(kx,ky): return comp(kx,ky).real
        def f2(kx,ky): return comp(kx,ky).imag

        resr,errorr = cp.integrate(f1, low, high, itermax=400)
        resi,errori = cp.integrate(f2, low, high, itermax=400)
        relerrorr = errorr/resr
        relerrori = errori/resi
        res = resr + 1j*resi
        res *= 4 # symmetry!
        res /= pow(2.0*numpy.pi,2)
        Glattloc[wi] = res

        print("Gloc=",Glattloc[wi], " relerrorr=", relerrorr, " relerrori=", relerrori)

    return Glattloc

def get_Glattloc_Hubbard(nw, ws, Sigmaw, mu, epsd, epsp, tpd, tpp, tppp):
    #here we are keeping the top band of H0, without Edc shift!!!
    def get_H(kx,ky):
        H = numpy.zeros((3,3),dtype=numpy.float_)    
        H[0,1:] = [2*tpd*numpy.sin(kx/2.),-2.*tpd*numpy.sin(ky/2.)]
        H[1,2] = -4.*tpp*numpy.sin(kx/2.)*numpy.sin(ky/2.)
        H += numpy.conjugate(numpy.transpose(H))
        numpy.fill_diagonal(H, [epsd, epsp+2*tppp*numpy.cos(kx),epsp+2*tppp*numpy.cos(ky)])
        return H

    low = [0.0, 0.0]
    high = [numpy.pi, numpy.pi]

    Glattloc = numpy.zeros(len(ws),dtype=numpy.complex_)
    for wi,(w,Sigma) in enumerate(zip(ws,Sigmaw)):   
        print("w=",w, " Sigma=", Sigma) 
        def comp(kx, ky):
            Gk = numpy.zeros((len(kx),1),dtype=numpy.complex_)
            for ki, ([the_kx],[the_ky]) in enumerate(zip(kx,ky)):
                H = get_H(the_kx,the_ky)
                Es,_ = numpy.linalg.eigh(H)
                Gk[ki,0] = 1/( w + mu - Es[2] - Sigma )
            return Gk

        def f1(kx,ky): return comp(kx,ky).real
        def f2(kx,ky): return comp(kx,ky).imag

        resr,errorr = cp.integrate(f1, low, high, itermax=400)
        resi,errori = cp.integrate(f2, low, high, itermax=400)
        relerrorr = errorr/resr
        relerrori = errori/resi
        res = resr + 1j*resi
        res *= 4 # symmetry!
        res /= pow(2.0*numpy.pi,2)
        Glattloc[wi] = res

        print("Gloc=",Glattloc[wi], " relerrorr=", relerrorr, " relerrori=", relerrori)

    return Glattloc

def get_Glattloc_Hubbard_using_DOS():    
    cmd = f"hilb -v -G -x {mu} -d {params['fndos']} {params['fnrs']} {params['fnis']} {params['fnrg']} {params['fnig']}"
    print(cmd)    
    os.system(cmd)
    _, Glattlocw = load_input(params['fnrg'], params['fnig'])
    return Glattlocw

get_Glattloc = {
    'Emery': lambda: get_Glattloc_Emery(nw, ws, Sigmaw, mu, epsd, epsp, tpd, tpp, tppp),
    'Hubbard': [
        lambda: get_Glattloc_Hubbard(nw, ws, Sigmaw, mu, epsd, epsp, tpd, tpp, tppp), 
        get_Glattloc_Hubbard_using_DOS
    ][params['use_Hubbard_DOS']]
}[params['model']]

print("about to get Glattloc")
t1 = time.time()
Glattlocw = get_Glattloc()
t2 = time.time()

if params['check'] or params['test']:
    print("   Glattloc norm (should be 1):",-(1/numpy.pi)*numpy.trapz(Glattlocw.imag, ws))
print("   get_Glattloc took",t2-t1,"seconds")

print("done")

if params['check'] or params['test']:
    print("printing out Glattloc",end="...")
    write_output(params['fnrg'], params['fnig'], ws, Glattlocw)
    print("done")

#if params['deal_with_Delta_shift']:
#    #get eps from param.eps
#    with open("param.eps", "r") as f:
#        epsd_plus_shift = float(f.read())
#
#    print("DELTA_SHIFT_BUSINESS: read epsd_plus_shift from param.eps:",epsd_plus_shift)
#else:
#    epsd_plus_shift = epsd

print("about to get Delta",end="...")
Deltaw = ws + mu - epsd - Sigmaw - 1.0/Glattlocw
print("done")

#if params['deal_with_Delta_shift']:
#    #extract shift and modify param.eps for the next iteration
#    shift = 0.5*(Deltaw.real[0]+Deltaw.real[-1])
#    print("DELTA_SHIFT_BUSINESS: got new shift from Delta:", shift)
#
#    with open("param.eps", "w") as f:
#        f.write("%g"%(epsd+shift))

print("about to printout Delta",end="...")
write_output(params['fnrd'], params['fnid'], ws, Deltaw)
print("ALL DONE")

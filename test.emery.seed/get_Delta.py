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
import time
from numpy import cos, sin
import cubepy as cp

from mpi4py import MPI
mpi = MPI.COMM_WORLD
mpirank = mpi.Get_rank()
mpisize = mpi.Get_size()

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
  'fnid': 'ImDelta.dat', #name of the output file for the imaginary part of Delta
  'fnrd': 'ReDelta.dat', #name of the output file for the real part of Delta
  'nk': 256, #number of k-points in a given direction, total number of points = nk**2
  'model': 'Emery', #'Hubbard' 
  'fnparam': 'tb_params.py', #file where the tight-binding parameters are stored, should be relevant for the model used
  'check': False, #in the commandline args use 0/1 for False/True
  'test': False, #if True, printout a Sigma to test the code
  'mpi_parallel': True
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

if mpirank==0: print("params:",params)
mpi_parallel = params['mpi_parallel']

# if test mode - printout a simple self-energy Sigma(w) = i Gamma, to be used as input
if params['test']:
    if mpirank==0: print("printing out test Sigma", end="...")
    test_ws = numpy.linspace(-5,5,1000,endpoint=True)    
    test_Sigma = -0.05j*numpy.ones((len(test_ws),))
    if mpirank==0: write_output(params['fnrs'], params['fnis'], test_ws, test_Sigma)
    if mpirank==0: print("done")
mpi.barrier()

# take the tb parameters from file
if mpirank==0: print("taking parameters from ",params['fnparam'])
exec(open(params['fnparam']).read())
if mpirank==0: print("done")

# load the self-energy
if mpirank==0: print("about to load Sigma from ",params['fnrs'], params['fnis'],end="...")
ws, Sigmaw = load_input(params['fnrs'], params['fnis'])
nw = len(ws)
if mpirank==0: print("done")

#prepare k-grid
nk = params['nk']
ks = numpy.linspace(numpy.pi/nk,2.*numpy.pi,nk,endpoint=False)


def get_Glattloc_Emery(nw, ws, Sigmaw, nk, ks, mu, epsd, epsp, tpd, tpp, tppp):
    I = numpy.eye(3) #3x3 unit matrix
    Idd = numpy.outer([1,0,0],[1,0,0]) #3x3 matrix with only the [0,0] component equal 1

    # implementing here arxiv:2308.09528, supp.mat. Eq.12; the hk matrix is purely real, and we are using that
    def hk(kx,ky):
        h = numpy.zeros((3,3),dtype=numpy.double)
        h[0,0] = epsd
        h[0,1] = h[1,0] = 2*tpd*sin(kx/2)
        h[0,2] = h[2,0] = -2*tpd*sin(ky/2)
        h[1,1] = epsp + 2*tppp*cos(kx)
        h[1,2] = h[2,1] = -4*tpp*sin(kx/2)*sin(ky/2)
        h[2,2] = epsp + 2*tppp*cos(ky)
        return h

    # get the Glatt loc
    Glattloc = numpy.zeros(len(ws),dtype=numpy.complex_)
    for wi,(w,Sigma) in enumerate(zip(ws,Sigmaw)):
        #if wi!=1000: continue
        if mpi_parallel:
            if wi % mpisize != mpirank: continue
        print("rank=", mpirank, " w=",w, " Sigma=", Sigma)

        low = [0.0, 0.0]
        high = [numpy.pi, numpy.pi]

        def comp(kx, ky):
            numerator = ((-epsp + mu + w - 2 * tppp * numpy.cos(kx)) * (-epsp + mu + w - 2 * tppp * numpy.cos(ky)) -
              16 * tpp**2 * numpy.sin(kx / 2)**2 * numpy.sin(ky / 2)**2)
            denominator = ((2 * tpd**2 * (-1 + numpy.cos(kx)) + (-epsd + mu - Sigma + w) * (-epsp + mu + w - 2 * tppp * numpy.cos(kx))) *
              (-epsp + mu + w - 2 * tppp * numpy.cos(ky)) +
              4 * tpd**2 * (epsp - mu - 2 * tpp - w + 2 * (tpp + tppp) * numpy.cos(kx)) * numpy.sin(ky / 2)**2 +
              16 * tpp * (-tpd**2 + tpp * (epsd - mu + Sigma - w)) * numpy.sin(kx / 2)**2 * numpy.sin(ky / 2)**2)
            res = numerator / denominator
            #print("res.shape:",res.shape, "kx.shape:",kx.shape)
            return res

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

def get_Glattloc_Hubbard(ws, Sigmaw, nk, ks, mu, t, tp, tpp):
    # prepare the tight-binding h(k)
    #implementing here arxiv:1705.08332, Eq.3
    epsk =   2*t*(cos(ks)[:,None]+cos(ks)[None,:])\
           + 4*tp*cos(ks)[:,None]*cos(ks)[None,:]\
           + 2*tpp*(cos(2*ks)[:,None]+cos(2*ks)[None,:])

    # get the Glatt loc
    Glattloc = numpy.zeros((len(ws),),dtype=numpy.complex_)
    for wi,(w,Sigma) in enumerate(zip(ws,Sigmaw)):
        if mpi_parallel:
            if wi % mpisize != mpirank: continue
        Glattloc[wi] = numpy.sum(1/( w + mu - epsk - Sigma ))
    Glattloc /= nk**2
    return Glattloc

get_Glattloc = {
  'Emery': lambda: get_Glattloc_Emery(nw, ws, Sigmaw, nk, ks, mu, epsd, epsp, tpd, tpp, tppp),
  'Hubbard': lambda: get_Glattloc_Hubbard(ws, Sigmaw, nk, ks, mu, t, tp, tpp)
}[params['model']]

if mpirank==0: print("about to get Glattloc")
t1 = time.time()
Glattlocw = get_Glattloc()
Glattlocw = mpi.reduce(Glattlocw)
t2 = time.time()

if mpirank!=0: quit() #we've done the parallel part, now only the master node does the work

if params['check'] or params['test']:
    print("   Glattloc norm (should be 1):",-(1/numpy.pi)*numpy.trapz(Glattlocw.imag, ws))
print("   get_Glattloc took",t2-t1,"seconds")

print("done")

if params['check'] or params['test']:
    print("printing out Glattloc",end="...")
    write_output('res/ReGlattloc.dat', 'res/ImGlattloc.dat', ws, Glattlocw)
    print("done")

print("about to get Delta",end="...")
Deltaw = ws + mu - Sigmaw - 1.0/Glattlocw
print("done")

print("about to printout Delta",end="...")
write_output(params['fnrd'], params['fnid'], ws, Deltaw)
print("ALL DONE")

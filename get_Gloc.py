######################
# Python 3
# 
# takes tight-binding parameters of the 2D lattice (including mu)
# takes local self-energy given in two files for real and imaginary parts 
# computes local Green's function on the lattice
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

def get_params_from_command_line(params): #params is dict
    #takes command line arguments in the form -key=val or -key. 
    #for each key a default value needs to be provided, and it will set the type of the parameter
    #bool parameters are set with integer (0=False, 1=True), but True can also be set without =1
    # -key=0 (False)
    # -key=1 or -key (True)
    
    for s in sys.argv:
        if s[0]!='-': continue
        if "=" not in s:
            key, val = s[1:], True
        else:
            key, val = s[1:].split("=")    

        assert key in params.keys(), "key %s not in params! default value must be provided!"

        if type(params[key])==bool:
            params[key] = bool(int(val))
        elif type(params[key])==str:
            params[key] = val
        elif type(params[key]) in [int, float]: params[key] = type(params[key])(val)
        else: params[key] = type(params[key])(eval(val))

        print("params[%s]: %s"%(key,str(params[key])))
    return params

# read the commandline arguments to set parameters
params = get_params_from_command_line({
  'fnis': 'res/imsigma.dat', #name of the input file for the imaginary part of Sigma
  'fnrs': 'res/resigma.dat', #name of the input file for the real part of Sigma
  'model': 'Emery', #'Hubbard' 
  'wstep': 0.,
  'fnparam': 'tb_params.py', #file where the tight-binding parameters are stored, should be relevant for the model used
  'test': False, #if True, printout a Sigma to test the code
  'skip': [],
  'components': [(0,0),(1,1)]#[(0,0),(0,1),(0,2),(1,1),(1,2),(2,2)]
})

if mpirank==0: print("params:",params)

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

wstep=params['wstep']
wis_to_do = [0]
for wi,w in enumerate(ws):
    if w>ws[wis_to_do[-1]]+wstep: wis_to_do.append(wi)
wis_to_do= numpy.array(wis_to_do)
if mpirank==0:
    print("wstep: %g, doing this many w-points: %d/%d"%(wstep, len(wis_to_do),len(ws)))

def get_Glattloc_Emery(nw, ws, Sigmaw, mu, epsd, epsp, tpd, tpp, tppp, component):
    # get the Glatt loc
    Glattloc = numpy.zeros(len(wis_to_do),dtype=numpy.complex_)
    #for wi,(w,Sigma) in enumerate(zip(ws,Sigmaw)):
    for wii,wi in enumerate(wis_to_do):
        if wii % mpisize != mpirank: continue        
        w = ws[wi]
        Sigma = Sigmaw[wi]

        print("rank=", mpirank, " w=",w, " Sigma=", Sigma)

        low = [0.0, 0.0]
        if component==(0,0):
            high = [numpy.pi, numpy.pi]
        else:
            high = [2*numpy.pi, 2*numpy.pi]

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
            M01squared = M01**2
            M02squared = M02**2
            M12squared = M12**2   
            M22 = w_plus_mu_minus_epsp - two_tppp*cosky
            det = M00*M11*M22 + 2*M01*M12*M02 - M00*M12squared - M01squared*M22 - M02squared*M11            
            return { 
               (0,0): lambda: M11*M22 - M12squared,
               (0,1): lambda: M02*M12 - M01*M22,
               (0,2): lambda: M01*M12 - M02*M11,
               (1,1): lambda: M00*M22 - M02squared,
               (1,2): lambda: M02*M01 - M00*M12,
               (2,2): lambda: M00*M11 - M01squared     
            }[component]()/det

        def f1(kx,ky): return comp(kx,ky).real
        def f2(kx,ky): return comp(kx,ky).imag

        resr,errorr = cp.integrate(f1, low, high, itermax=400)
        resi,errori = cp.integrate(f2, low, high, itermax=400)
        relerrorr = errorr/resr
        relerrori = errori/resi
        res = resr + 1j*resi
        if component==(0,0): res *= 4 # symmetry!
        res /= (2.*numpy.pi)**2
        Glattloc[wii] = res

        print("Gloc=",Glattloc[wii], " relerrorr=", relerrorr, " relerrori=", relerrori)

    return Glattloc

for component in params['components']:
    if component in params['skip']:
        print("skipping component:",component)
        continue
    orbs = ['d','px','py']
    s = "%s.%s"%(orbs[component[0]],orbs[component[1]])
    t1 = time.time() 
    Glattlocw = get_Glattloc_Emery(nw, ws, Sigmaw, mu, epsd, epsp, tpd, tpp, tppp, component)
    Glattlocw = mpi.reduce(Glattlocw)
    t2 = time.time()
    if mpirank==0: 
       print("   get_Glattloc component %s took %g seconds"%(s,t2-t1))
       write_output('res/ReGlattloc.%s.dat'%s, 'res/ImGlattloc.%s.dat'%s, ws[wis_to_do], Glattlocw)


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
import h5py
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
  'test': False, #if True, printout a Sigma to test the code,
  'nk': 32,
  'skip': []
})

if mpirank==0: print("params:",params)


def get_k_points(equidistant_checkpoints=False):
    nk = params['nk']
    ks = numpy.linspace(0,2*numpy.pi,nk, endpoint=False)

    dk = 2.0*numpy.pi/nk
    dk_diag = numpy.sqrt(2.0)*dk if not equidistant_checkpoints else dk
    #print ks
    klist = []
    xs= []
    path_covered = 0.0

    for i in range(nk//2):
        klist.append((ks[i],ks[i]))
        xs.append(path_covered)
        path_covered += dk_diag

    for i in reversed(range(1,nk//2+1)):
        klist.append((ks[i],ks[nk//2]))
        xs.append(path_covered)
        path_covered += dk

    for i in range(nk//2):
        klist.append((ks[i],ks[nk//2-i]))
        xs.append(path_covered)
        path_covered += dk_diag

    for i in reversed(range(nk//2+1)):
        klist.append((ks[i],ks[0]))
        xs.append(path_covered)
        path_covered += dk
    xtcs=[xs[0],xs[nk//2],xs[nk],xs[3*nk//2],xs[-1]]
    xtcs_labels = [r"$(0,0)$",r"$(\pi,\pi)$",r"$(0,\pi)$",r"$(\pi,0)$",r"$(0,0)$"]
    return numpy.array(xs),numpy.array(klist), xtcs, xtcs_labels

xs, klist, xtcs, xtcs_labels = get_k_points()

if mpirank==0: print("klist:",klist)


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

def get_Gkw_Emery(nw, ws, Sigmaw, mu, epsd, epsp, tpd, tpp, tppp, component):
    # get the Glatt loc
    Gkw = numpy.zeros((len(klist),len(wis_to_do)), dtype=numpy.complex_)
    #for wi,(w,Sigma) in enumerate(zip(ws,Sigmaw)):
    for wii,wi in enumerate(wis_to_do):
        if wii % mpisize != mpirank: continue        
        w = ws[wi]
        Sigma = Sigmaw[wi]

        print("rank=", mpirank, " w=",w, " Sigma=", Sigma)

        low = [0.0, 0.0]
        #if component==(0,0):
        #high = [numpy.pi, numpy.pi]
        #else:
        high = [2*numpy.pi, 2*numpy.pi]

        w_plus_mu = w+mu
        w_plus_mu_minus_epsp = w_plus_mu-epsp 
        minus_two_tpd = -2*tpd
        two_tppp = 2*tppp
        M00 = w_plus_mu-epsd-Sigma
        def comp(kx, ky):                     
            #print("kx.shape before:",kx.shape)
            #print("ky.shape before:",ky.shape)

            kx = kx[:,0]
            ky = ky[:,0]
            #print("kx.shape:",kx.shape)
            #print("ky.shape:",ky.shape)
            coskx = numpy.cos(kx)
            cosky = numpy.cos(ky)
            sinkxhalf = numpy.sin(kx/2)
            sinkyhalf = numpy.sin(ky/2)                      

            hk = numpy.zeros((3,3,len(kx)),dtype=numpy.double)
            hk[0,0,:] = epsd
            hk[0,1,:] = hk[1,0,:] = 2*tpd*sinkxhalf
            hk[0,2,:] = hk[2,0,:] = -2*tpd*sinkyhalf
            hk[1,1,:] = epsp + 2*tppp*coskx
            hk[1,2,:] = hk[2,1,:] = -4*tpp*sinkxhalf*sinkyhalf
            hk[2,2,:] = epsp + 2*tppp*cosky  
            Pk = numpy.zeros((3,3,len(kx)),dtype=numpy.complex_)
            Pkinv = numpy.zeros((3,3,len(kx)),dtype=numpy.complex_)
            for ki in range(len(kx)):
                _, Pk[:,:,ki] = numpy.linalg.eigh(hk[:,:,ki])  
                #Pkinv[:,:,ki] = numpy.conjugate(Pk[:,:,ki].T) #should be fine
                Pkinv[:,:,ki] = numpy.linalg.inv(Pk[:,:,ki])
            M01 = -hk[0,1,:]
            M02 = -hk[0,2,:]
            M11 = w_plus_mu-hk[1,1,:]
            M12 = -hk[1,2,:]
            M01squared = M01**2
            M02squared = M02**2
            M12squared = M12**2   
            M22 = w_plus_mu - hk[2,2,:]
            det = M00*M11*M22 + 2*M01*M12*M02 - M00*M12squared - M01squared*M22 - M02squared*M11            
            Gk00 = (M11*M22 - M12squared)/det
            Gk01 = (M02*M12 - M01*M22)/det
            Gk02 = (M01*M12 - M02*M11)/det
            Gk11 = (M00*M22 - M02squared)/det
            Gk12 = (M02*M01 - M00*M12)/det
            Gk22 = (M00*M11 - M01squared)/det

            Gk = numpy.zeros((3,3,len(kx)),dtype=numpy.complex_)
            Gk[0,0,:]=Gk00
            Gk[0,1,:]= Gk[1,0,:] = Gk01
            Gk[0,2,:]= Gk[2,0,:] = Gk02
            Gk[1,1,:]=Gk11
            Gk[1,2,:]=Gk[2,1,:]=Gk12
            Gk[2,2,:]=Gk22

            res = [] 
            for ki in range(len(kx)):
                res.append([ (Pkinv[:,:,ki] @ Gk[:,:,ki] @ Pk[:,:,ki])[component[0],component[1]] ])
            res = numpy.array(res)
            #print("res.shape:",res.shape)
            return res 

        kxs = numpy.array([[k[0]] for k in klist])
        kys = numpy.array([[k[1]] for k in klist])
        Gkw[:,wii] = comp(kxs,kys)[:,0]
    return Gkw

for alpha in [0,1,2]:
    if alpha in params['skip']: 
        print("skipping:",alpha)
        continue
    s = "alpha%d"%(alpha)
    t1 = time.time() 
    Gkw = get_Gkw_Emery(nw, ws, Sigmaw, mu, epsd, epsp, tpd, tpp, tppp, (alpha,alpha))
    Gkw = mpi.reduce(Gkw)
    t2 = time.time()
    if mpirank==0: 
       print("   get_Glattloc component %s took %g seconds"%(s,t2-t1))
       f = h5py.File("res/Gkw.%s.h5"%s, 'w')
       f.create_dataset("Akw", data=-1/numpy.pi*Gkw.imag)
       f.create_dataset("ReGkw", data=Gkw.real)
       f.create_dataset("ImGkw", data=Gkw.imag)
       f.create_dataset("ws", data=ws[wis_to_do])
       f.create_dataset("klist", data=klist)
       f.create_dataset("xs", data=xs)
       f.create_dataset("xtcs", data=xtcs)
       f.create_dataset("xtcs_labels",data=xtcs_labels)
       f.close()


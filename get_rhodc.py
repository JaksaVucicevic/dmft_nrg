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
import os
import time
from numpy import cos, sin
import cubepy as cp

from mpi4py import MPI
mpi = MPI.COMM_WORLD
mpirank = mpi.Get_rank()
mpisize = mpi.Get_size()

def write_real_output(fn, xs, ys):
    #take function y given on a real grid x and printout to two separate files real and imaginary part of y
    data = numpy.array(list(zip(xs,ys)))
    numpy.savetxt(fn, data)


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

def load_real_input(fn):
    #load complex function y of real argument x written in two separate files for real and imaginary part of y
    data = numpy.loadtxt(fn)
    xs = data[:,0]
    ys = data[:,1]
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
  'fnPhi': 'Phi.dat', #only needed if model is Hubbard
  'fnis': 'res/imsigma.dat', #name of the input file for the imaginary part of Sigma
  'fnrs': 'res/resigma.dat', #name of the input file for the real part of Sigma
  'model': 'auto',#'Emery', #'Hubbard' #'auto'
  'wstep': 0.00005,
  'fnparam': 'tb_params.py', #file where the tight-binding parameters are stored, should be relevant for the model used
  'test': False, #if True, printout a Sigma to test the code
  'skip': [],
  'components_and_sym_prefs': [
    ((0, 1, 0, 1), 4),
    #((0, 1, 0, 1), 2),
    #((0, 1, 1, 0), 2),
    ((0, 1, 1, 1), 4),
    ((0, 1, 1, 2), 8),
    #((0, 1, 1, 2), 4),
    #((0, 1, 2, 1), 4),
    ((1, 1, 1, 1), 1),
    ((1, 1, 1, 2), 4),
    ((1, 2, 1, 2), 4),
    #((1, 2, 1, 2), 2),
    #((1, 2, 2, 1), 2)
  ],
  'T': float('nan'),
  'hbar': 1,
  'c': 1,
  'intnFp_threshold': 0.0005
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


# prepare the grid; take into account the desired number of points, and the temeprature that set's the energy window

if params['T']!=params['T']:    
    #take the value from param.loop
    try:
        with open("param.loop","r") as f:
            lines = f.readlines()
        found = False
        for line in lines: 
            if line[:2]=="T=":
                params['T'] = float(line[2:].strip(" \n"))
                found=True
                break
        if not found: raise 
        print("loaded T from param.loop:",params['T'])
    except:
        print("ERROR: could not read T from param.loop")
        quit()
beta = 1.0/params['T']
intnFp_threshold = params['intnFp_threshold']

nFp = -beta * numpy.exp(beta*ws)/(1.+numpy.exp(beta*ws))**2 
nFp[nFp!=nFp]=0.

write_real_output('res/nFp.dat', ws, nFp)

#get wi_start
wi_start = None
for wi,w in enumerate(ws):
    if wi>3:
        intnFp = -numpy.trapz(nFp[:wi],ws[:wi])
        #print("w: %g intnFp=%g"%(w,intnFp))
        if intnFp>intnFp_threshold:
            wi_start = wi; break
if wi_start is None:
    print("ERROR, could not find wi_start")
else: print("w_start:",ws[wi_start])
#get wi_end
for wi,w in reversed(list(enumerate(ws))):
    if wi<len(ws)-3:
        intnFp = -numpy.trapz(nFp[:wi],ws[:wi])
        #print("w: %g intnFp=%g"%(w,intnFp))
        if intnFp<1-intnFp_threshold:
            wi_end = wi; break
if wi_end is None:
    print("ERROR, could not find wi_end")
else: print("w_end:",ws[wi_end])

wstep=params['wstep']
wis_to_do = [wi_start]

for wi in range(wi_start,wi_end+1):
    #print("wi:",wi)
    w = ws[wi]
    #print("w:",w)
    #print("wis_to_do:",wis_to_do)
    if w>ws[wis_to_do[-1]]+wstep: 
        wis_to_do.append(wi)
wis_to_do= numpy.array(wis_to_do)
ws_to_do = ws[wis_to_do]

if mpirank==0:
    print("wstep: %g, doing this many w-points: %d/%d, in the range %.3f,%.3f"%(wstep, len(wis_to_do),len(ws),ws_to_do[0],ws_to_do[-1]))

if params['model']=='auto':
    params['model']=['Emery','Hubbard'][int(os.path.exists(params['fnPhi']))]
    print("automatic model determination based on the existence of Phi file: model=",params['model'])

if params['model']=='Emery':
    def get_integrandw_Emery(component):
        # get the Glatt loc
        integrandw = numpy.zeros(len(wis_to_do), dtype=numpy.float_)
        l1,l2,l3,l4 = component
        print("doing component: ",l1,l2,l3,l4)
        #for wi,(w,Sigma) in enumerate(zip(ws,Sigmaw)):
        for wii,wi in enumerate(wis_to_do):
            if wii % mpisize != mpirank: continue        
            w = ws[wi]
            Sigma = Sigmaw[wi]

            #print("rank=", mpirank, " w=",w, " Sigma=", Sigma)

            low = [0.0, 0.0]
            high = [2*numpy.pi, 2*numpy.pi]

            w_plus_mu = w+mu
            w_plus_mu_minus_epsp = w_plus_mu-epsp 
            two_tpd = 2*tpd
            two_tppp = 2*tppp
            M00 = w_plus_mu-epsd-Sigma
            def Gk(kx, ky, l, lp):  
                if l>lp: l,lp=lp,l     #we know that the matrix will be symmetric      
                coskx = numpy.cos(kx)
                cosky = numpy.cos(ky)
                sinkxhalf = numpy.sin(kx/2)
                sinkyhalf = numpy.sin(ky/2)                      
                M01 = -two_tpd*sinkxhalf
                M02 = two_tpd*sinkyhalf
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
                }[(l,lp)]()/det

            def vk(kx,ky,l,lp):
                if l>lp: l,lp=lp,l    
                return {
                  (0,1): lambda: tpd*numpy.cos(kx/2),
                  (1,1): lambda: -2*tppp*numpy.sin(kx),
                  (1,2): lambda: -2*tpp*cos(kx/2)*sin(ky/2)   
                }[(l,lp)]()


            def f(kx,ky): 
                return vk(kx,ky,l1,l2)*vk(kx,ky,l3,l4)*(Gk(kx,ky,l2,l3).imag)*(Gk(kx,ky,l4,l1).imag)

            res,error = cp.integrate(f, low, high, itermax=400)
            relerror = error/res        
            res /= (2.*numpy.pi)**2

            integrandw[wii] = -res*nFp[wi] #!!!!!!!!!!!!! why do we need a minus sign here???
            #print("integrandw=",integrandw[wii], " relerror=", relerror, "nFp[wii]", nFp[wi], "res:",res)

        return integrandw

    ti =  time.time() 
    sigma_dc = 0.
    for component, sym_pref in params['components_and_sym_prefs']:
        if component in params['skip']:
            print("skipping component:",component)
            continue
        s = "%d%d%d%d"%component
        t1 = time.time() 
        integrandw = get_integrandw_Emery(component)
        integrandw = mpi.reduce(integrandw)
        integral = numpy.trapz(integrandw,ws_to_do)
        sigma_dc += integral*sym_pref    
        t2 = time.time()
        if mpirank==0: 
           print("   get_rhodc component %s took %g seconds, partial_sigma_dc=%.2le "%(s,t2-t1, integral*sym_pref))
           write_real_output('res/sigmadc_integrandw.%s.dat'%s, ws[wis_to_do], integrandw)
    tf = time.time() 
    print("    full calculation took %g seconds, total sigma_dc=%.2le, total rho_dc=%.2le"%(tf-ti, sigma_dc, 1.0/sigma_dc))
elif params['model']=='Hubbard':
    def get_integrandw_Hubbard():
        #load Phi
        ws_Phi, Phiw = load_real_input(params['fnPhi'])

        # get the Glatt loc
        integrandw = numpy.zeros(len(wis_to_do), dtype=numpy.float_)
        #for wi,(w,Sigma) in enumerate(zip(ws,Sigmaw)):
        for wii,wi in enumerate(wis_to_do):
            if wii % mpisize != mpirank: continue        
            w = ws[wi]
            Sigma = Sigmaw[wi]

            #print("rank=", mpirank, " w=",w, " Sigma=", Sigma)

            low = [0.0, 0.0]
            high = [2*numpy.pi, 2*numpy.pi]
            
            def Geps(eps):  
                return 1.0/(w+mu-eps-Sigma)

            res = numpy.trapz( Phiw*(Geps(ws_Phi).imag)**2, ws_Phi )

            integrandw[wii] = -res*nFp[wi] #!!!!!!!!!!!!! why do we need a minus sign here???
            #print("integrandw=",integrandw[wii], " relerror=", relerror, "nFp[wii]", nFp[wi], "res:",res)

        return integrandw

    ti =  time.time() 
    integrandw = get_integrandw_Hubbard()
    integrandw = mpi.reduce(integrandw)
    sigma_dc = numpy.trapz(integrandw,ws_to_do)
    tf = time.time()     
    if mpirank==0: 
        print("    full calculation took %g seconds, total sigma_dc=%.2le, rho_dc=%.2le"%(tf-ti, sigma_dc, 1.0/sigma_dc))
        write_real_output('res/sigmadc_integrandw.dat', ws_to_do, integrandw)
else:
    print("ERROR: unknown model",params['model'])
    quit()
if mpirank==0: 
    with open("res/rhodc.txt","w") as f:
        f.write("%g"%(1./sigma_dc))


import sys
import numpy
import itertools

# default parameters
params = {
  'fnPhi': 'Phi.dat', #name of the input file for the imaginary part of Sigma
  'fnparam': 'tb_params.py', #file where the tight-binding parameters are stored, should be relevant for the model used
  'alpha': 2, #in the commandline args use 0/1 for False/True 
  'nk': 4000,
  'margine': 0.02,
  'nbins_factor': 0.1  
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

globals().update(params)

# take the tb parameters from file
print("taking parameters from ",fnparam)
exec(open(fnparam).read())
print("done")

ks = numpy.linspace(0, 2*numpy.pi, nk, endpoint=False)
dk = ks[1]-ks[0]

nbins = int(numpy.sqrt(nk**2.0)*nbins_factor)

Hk = numpy.zeros((nk,nk,3,3),dtype=numpy.float_)    
Hk[:,:,0,0] = epsd
Hk[:,:,1,1] = epsp+2*tppp*numpy.cos(ks)[:,None]
Hk[:,:,2,2] = epsp+2*tppp*numpy.cos(ks)[None,:]
Hk[:,:,0,1] = 2*tpd*numpy.sin(ks/2.)[:,None]
Hk[:,:,0,2] = -2.*tpd*numpy.sin(ks/2.)[None,:]
Hk[:,:,1,2] = -4.*tpp*numpy.sin(ks/2.)[:,None]*numpy.sin(ks/2.)[None,:]
Hk += numpy.conjugate(numpy.transpose(Hk,axes=(0,1,3,2)))

print("prepared Hk")

Ek = numpy.zeros((nk,nk,3),dtype=numpy.float_)    
for kxi,kyi in itertools.product(range(nk),repeat=2):
    Ek[kxi,kyi,:], _ = numpy.linalg.eigh(Hk[kxi,kyi])

print("prepared Ek")

vxk = 0.5*(numpy.roll(Ek[...,2], shift=1, axis=0)-numpy.roll(Ek[...,2], shift=-1, axis=0))/dk

print("prepared vxk")

#get energy window and give some margine

mn = numpy.amin(Ek[:,:,alpha].real)
mx = numpy.amax(Ek[:,:,alpha].real)

center = 0.5*(mn+mx)
width = mx-mn
mx = center + (1.+margine)*width/2
mn = center - (1.+margine)*width/2

print("range: mn,mx:",mn,mx)

#do the histogram

hist, bin_edges = numpy.histogram(
    Ek[:,:,alpha],
    bins=nbins,
    range=(mn,mx),
    weights = vxk**2
)
omegas = numpy.array([bin_edges[0]]+[0.5*(w +  bin_edges[wi+1]) for wi,w in enumerate( bin_edges[:-1])]+[bin_edges[-1]])
domega = omegas[1]-omegas[0]
Phi = numpy.array([0]+list(hist)+[0])/(nk**2*domega)

numpy.savetxt(fnPhi,numpy.concatenate((omegas[:,None], Phi[:,None]), axis=1))
    

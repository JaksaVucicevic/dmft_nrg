import numpy
import scipy
import itertools
import sys
import os
import time

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
        else:
            #print("val:",val) 
            params[key] = type(params[key])(eval(val))

        print("params[%s]: %s"%(key,str(params[key])))
    return params

# read the commandline arguments to set parameters
params = get_params_from_command_line({
  'clean_up': False, #if true, source will be cleaned up first
  'source': 'test', #name of the source folder from which the new calculation will start
  'dest': 'test_new', #name of the destination folder, where the new calculation will be
  'auto_dest': False, #if True, overrides dest, and dest is set as source+key.val.key.val.key.val... for the keys set.
  'from_scratch': False, #if True, deletes all *.dat files in dest
  'keys': [],
  'vals': [],
})

globals().update(params)
if clean_up:
    os.system(f"cd {source}; ../scripts/bin/CLEANUPDMFT;")
    print("cleaned up source")

if auto_dest:
    dest = source+"."+".".join([key+"%.3f"%val for key,val in zip(keys,vals)])

print("dest=",dest)

os.system(f"cp -Tr {source} {dest}")
os.chdir(dest)
os.system("../scripts/bin/CLEANUPDMFT_FOR_RESTART")
if from_scratch: 
    os.system("GLOBIGNORE=DOS.dat; rm *.dat *.err *.out")
    print("will be doing things from scratch")
#os.system(f"cp {source}/DOS.dat {dest}/")

print("cleaned up dest")

for key,val in zip(keys,vals):
  if key in ['T','U']: continue
  os.system('sed -i "/%s = /c\%s = %g" tb_params.py'%(key,key,val))

print("prepared tb_params.py in dest")

if 'mu' in keys:
    os.system("rm param.mu")
    os.system("touch param.mu")
    os.system('echo %g >> param.mu'%(vals[keys.index('mu')]))
    print("prepared param.mu in dest")
if 'epsd' in keys:
    os.system("rm param.eps")
    os.system("touch param.eps")
    os.system('echo %g >> param.eps'%(vals[keys.index('epsd')]))
    print("prepared param.eps in dest")
if 'T' in keys:
    os.system(r'sed -i "/T=/c\T=%g" param.loop'%(vals[keys.index('T')]))
    print("set T in param.loop in dest")
if 'U' in keys:
    os.system(r'sed -i "/U1=/c\U1=%g" param.loop'%(vals[keys.index('U')]))
    print("set U in param.loop in dest")

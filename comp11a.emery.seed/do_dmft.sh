module purge
module load NRGLjubljana/master-foss-2023b

module list

export PATH=$HOME/NRGRUN/TEST/scripts/tools:$PATH
echo $PATH

export EASYBUILD_MODULES_TOOL=EnvironmentModules
export EASYBUILD_MODULE_SYNTAX=Tcl

./DMFT

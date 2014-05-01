qmmm-hydration
==============

QM/MM hydration free energy playground

## Manifest

```
setenvs.tcsh - set up environment variables
setup-hydration-simulations.py - set up hydrated and vacuum simulations for AMBER
ligandtools.py - module used by setup-hydration-simulations.py
```

## To run

```
# Set environment variables
source ../openeye/openeye.csh
source ../amber/amber.csh
source ../anaconda/anaconda.csh
# Clean up molecules directory
rm -rf molecules
# Set up simulations
python setup-hydration-simulations.py
```

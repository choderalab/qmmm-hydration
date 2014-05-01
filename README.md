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
source setenvs.csh # for tcsh/csh
rm -rf molecules ; python setup-hydration-simulations.py
```

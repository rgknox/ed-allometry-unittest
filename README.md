# ed-allometry-unittest
This is a python (v2) framework for unit-testing various allometry schemes in the ED model.  These are offline (not coupled to the CLM-ED, ED2 or other) tests.  The python script drive_allomtests.py calls two F90 shared libraries, EDAllomMod.o and EDAllomUnitWrap.o.  The communication between python and F90 shared libraries are handled via the "ctypes" classes.

##Contents:

drive_allomtests.py: The python driver script

EDAllomMod.F90: The allometry library that is being prototyped, ultimately to make its way into clm-ed.

EDAllomUnitWrap.F90: A library that allocates the pft derived type and handles some classes.

allom_params.xml: An XML file containing the PFT parameters relevant to allometry.

simple_build.sh: A simple bash script that compiles the two libraries.

##Deprecated:

allom_lib_v3.m: The library that is mirrored in moderen fortran inside the model.
drive_allomtests.m:  The user invokes tests through this script.  This is the driver script and contains some control parameters.
gen_param_instance.m:  This is the control card that defines the parameters associated with the different schemes.

The following are various plotting functions that could un-doubtedly be re-written more efficiently into one script with variable controls:

functions/plot_multicase_bag.m
functions/plot_multicase_blmax.m
functions/plot_multicase_blmax_o_dbagdh.m
functions/plot_multicase_blmaxdi.m
functions/plot_multicase_bsap.m
functions/plot_multicase_bsapbag.m
functions/plot_multicase_bsapid.m
functions/plot_multicase_dbhe.m
functions/plot_multicase_h.m
functions/plot_singlecase_cfractions.m

gen_param_instance.m

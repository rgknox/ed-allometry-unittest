# ed-allometry-unittest
This is a matlab unit-testing framework for unit-testing various allometry schemes in the ED model.  These are offline (not coupled to the CLM-ED, ED2 or other) tests.

Contents:

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

import numpy
from ctypes import byref, cdll, c_int, c_double, c_char_p

# edallomunitwrap_MOD_edecophysconalloc

fallomlib = cdll.LoadLibrary('./EDAllomUnitWrap.o')

numpft = c_int(10)

fallomlib.__edallomunitwrap_MOD_edecophysconalloc(byref(numpft))

pname = c_char_p("bl_min")
ipft  = c_int(2)
rval  = c_double(0.04)
ival  = c_int(0)
cval  = c_char_p("p")

fallomlib.__edallomunitwrap_MOD_edecophysconpyset(byref(pname),byref(ipft),byref(rval),byref(ival),byref(cval))



#fadd = cdll.LoadLibrary('./add.so')

#arg1 = c_double(2.0)
#arg2 = c_double(4.0)

#arg3 = c_int(2)
#arg4 = c_int(4)

#fadd.__add_MOD_initpatch(byref(arg4))

#print fadd.__add_MOD_frand(byref(arg4),byref(arg3))

#print fadd.__add_MOD_fiadd(byref(arg3),byref(arg4))

#print fadd.__add_MOD_fdadd(byref(arg1),byref(arg2))


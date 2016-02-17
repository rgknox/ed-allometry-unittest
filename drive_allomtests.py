import numpy
import ctypes
from ctypes import * #byref, cdll, c_int, c_double, c_char_p, c_long
import xml.etree.ElementTree as ET

pft_xml_file = "allom_params.xml"
allom_wrap_object = "./EDAllomUnitWrap.o"
allom_lib_object = "./EDAllomMod.o"


# ==============================================================================
# The following lists define what variables are to be expected for the
# different datatypes (double and int).  It is expected that the XML
# (and netcdf later) and the F90 wrapper have the same entries or more
# ==============================================================================

# These are the expected PFT parameters that are double precision
expt_par_dp = ['c2b','eclim','bl_min','h_max','h_min','slatop', \
               'd_adult','d_sap','f2l_ratio','agb_fraction','latosa_int', \
               'latosa_slp','d2h1_ad','d2h2_ad','d2h3_ad','d2h1_sap', \
               'd2h2_sap','d2bl1_ad','d2bl2_ad','d2bl3_ad','d2bl1_sap', \
               'd2bl2_sap','d2bag1','d2bag2','wood_density']

# These are the expected PFT parameters that are integers
expt_par_int = ['hallom_mode','lallom_mode','fallom_mode','aallom_mode', \
                'callom_mode','sallom_mode']


# ==============================================================================
# Traverse the XML tree and fill the target data structure pftparms
# pftparms is a list of dictionaries
# ==============================================================================

pftroot = ET.parse(pft_xml_file).getroot()
pftparms = []

numpft = 0
for elem in pftroot.iter('pft'):
    plist = {}
    for iv in expt_par_dp:
        pftelem = elem.find(iv)
        if pftelem!=None:
            plist[pftelem.tag] = float(pftelem.text)
        else:
            trypftroot = pftroot.find(iv)
            if trypftroot!=None:
                plist[trypftroot.tag] = float(trypftroot.text)
            else:
                print('Could not find ',iv,' in the XML')
    for iv in expt_par_int:
        pftelem = elem.find(iv)
        if pftelem!=None:
            plist[pftelem.tag] = int(pftelem.text)
        else:
            trypftroot = pftroot.find(iv)
            if trypftroot!=None:
                plist[trypftroot.tag] = int(trypftroot.text)
            else:
                print('Could not find ',iv,' in the XML')
    numpft += 1
    pftparms.append(plist)


# ==============================================================================
# Load the fortran allometry library using python's ctypes library
# ==============================================================================

f90wraplib = ctypes.CDLL(allom_wrap_object,mode=ctypes.RTLD_GLOBAL)
f90funclib = ctypes.CDLL(allom_lib_object,mode=ctypes.RTLD_GLOBAL)

# ==============================================================================
# Allocate fortran PFT arrays
# ==============================================================================

iret=f90wraplib.__edallomunitwrap_MOD_edecophysconalloc(byref(c_int(numpft)))


# ==============================================================================
# Populate the Fortran PFT structure
# ==============================================================================

# First set the arg types
f90wraplib.__edallomunitwrap_MOD_edecophysconpyset.argtypes = \
    [POINTER(c_int),POINTER(c_double),POINTER(c_int),c_char_p,c_long]

for ipft in range(numpft):
    elem=pftparms[ipft]
    for pname in expt_par_dp:
        print 'Sending to F90: {0} = {1}'.format(pname,elem[pname])
        iret=f90wraplib.__edallomunitwrap_MOD_edecophysconpyset(c_int(ipft+1),c_double(elem[pname]),c_int(0),c_char_p(pname),c_long(len(pname)))
    for pname in expt_par_int:
        print 'Sending to F90: {0} = {1}'.format(pname,elem[pname])
        iret=f90wraplib.__edallomunitwrap_MOD_edecophysconpyset(c_int(ipft+1),c_double(-999.9),c_int(elem[pname]),c_char_p(pname),len(pname))



# Some testing constants
ndbh = 2000
maxdbh = 150.0

# =========================================================================
# Generate a vector of diameters that starts at the smallest known diameter
# and extends to 150cm

# =========================================================================
# Initialize Output Arrays

blmaxi  = numpy.zeros((numpft,ndbh))
blmaxd  = numpy.zeros((numpft,ndbh))

bfrmax = numpy.zeros((numpft,ndbh))
hi     = numpy.zeros((numpft,ndbh))
bagi   = numpy.zeros((numpft,ndbh))
bagd   = numpy.zeros((numpft,ndbh))
dbh    = numpy.zeros((numpft,ndbh))
bcr    = numpy.zeros((numpft,ndbh))
bsapi  = numpy.zeros((numpft,ndbh))
bsapd  = numpy.zeros((numpft,ndbh))
bdead  = numpy.zeros((numpft,ndbh))
dbhe   = numpy.zeros((numpft,ndbh))

blmax_o_dbagdh = numpy.zeros((numpft,ndbh))

# Minimum DBH and maximum DBH are diagnosed
# ==============================================================================

f90_h2d = f90funclib.__edallommod_MOD_h2d_allom

for ipft in range(numpft):
    ch_min = c_double(pftparms[ipft]['h_min'])
    cd = c_double(-9.0)
    cdddh = c_double(-9.0)
    cipft = c_int(ipft+1)
    iret=f90_h2d(byref(c_double(pftparms[ipft]['h_min'])),byref(cipft),byref(cd),byref(cdddh))
    pftparms[ipft].update({'d_min':cd.value})
    print 'h_min of {!r} generated d_min of {!r}'.format(pftparms[ipft]['h_min'],pftparms[ipft]['d_min'])

    iret=f90_h2d(byref(c_double(pftparms[ipft]['h_max'])),byref(cipft),byref(cd),byref(cdddh))
    pftparms[ipft].update({'d_max':cd.value})
    print 'h_max of {!r} generated d_max of {!r}'.format(pftparms[ipft]['h_max'],pftparms[ipft]['d_max'])


#    dbh(ip,:)               = linspace(pftcon.dbh_min(ip),maxdbh,ndbh);


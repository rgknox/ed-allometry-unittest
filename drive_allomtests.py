import numpy as np
import matplotlib.pyplot as plt
import ctypes
from ctypes import * #byref, cdll, c_int, c_double, c_char_p, c_long
import xml.etree.ElementTree as ET

pft_xml_file = "allom_params.xml"
allom_wrap_object = "./EDAllomUnitWrap.o"
allom_lib_object = "./EDAllomMod.o"

def wait():
    msvcrt.getch()

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
                print('py: Could not find ',iv,' in the XML')
    for iv in expt_par_int:
        pftelem = elem.find(iv)
        if pftelem!=None:
            plist[pftelem.tag] = int(pftelem.text)
        else:
            trypftroot = pftroot.find(iv)
            if trypftroot!=None:
                plist[trypftroot.tag] = int(trypftroot.text)
            else:
                print('py: Could not find ',iv,' in the XML')
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
        print 'py: sending to F90: {0} = {1}'.format(pname,elem[pname])
        iret=f90wraplib.__edallomunitwrap_MOD_edecophysconpyset(c_int(ipft+1), \
                    c_double(elem[pname]),c_int(0),c_char_p(pname),c_long(len(pname)))
    for pname in expt_par_int:
        print 'py: sending to F90: {0} = {1}'.format(pname,elem[pname])
        iret=f90wraplib.__edallomunitwrap_MOD_edecophysconpyset(c_int(ipft+1), \
                    c_double(-999.9),c_int(elem[pname]),c_char_p(pname),len(pname))



# Some testing constants
ndbh = 200
maxdbh = 150.0

# =========================================================================
# Generate a vector of diameters that starts at the smallest known diameter
# and extends to 150cm

# =========================================================================
# Initialize Output Arrays

blmaxi  = np.zeros((numpft,ndbh))
blmaxd  = np.zeros((numpft,ndbh))

bfrmax = np.zeros((numpft,ndbh))
hi     = np.zeros((numpft,ndbh))
bagi   = np.zeros((numpft,ndbh))
bagd   = np.zeros((numpft,ndbh))
dbh    = np.zeros((numpft,ndbh))
bcr    = np.zeros((numpft,ndbh))
bsapi  = np.zeros((numpft,ndbh))
bsapd  = np.zeros((numpft,ndbh))
bdead  = np.zeros((numpft,ndbh))
dbhe   = np.zeros((numpft,ndbh))

blmax_o_dbagdh = np.zeros((numpft,ndbh))

# Minimum DBH and maximum DBH are diagnosed
# ==============================================================================

f90_h2d   = f90funclib.__edallommod_MOD_h2d_allom  #(h,ipft,d,dddh)
f90_h     = f90funclib.__edallommod_MOD_h_allom    #(d,ipft,h,dhdd)
f90_bag   = f90funclib.__edallommod_MOD_bag_allom  #(d,h,ipft,bag,dbagdd)
f90_blmax = f90funclib.__edallommod_MOD_blmax_allom  #(d,h,ipft,blmax,dblmaxdd)
f90_bsap  = f90funclib.__edallommod_MOD_bsap_allom  #(d,h,blmax,dblmaxdd,dhdd,ipft,bsap,dbsapdd)
f90_bcr   = f90funclib.__edallommod_MOD_bcr_allom  #(d,bag,dbagdd,ipft,bcr,dbcrdd)
f90_bfrmax= f90funclib.__edallommod_MOD_bfrmax_allom  #(d,blmax,dblmaxdd,ipft,bfrmax,dbfrmaxdd)
#(bag,bcr,blmax,bsap,dbagdd,dbcrdd,dblmaxdd,dbsapdd,bdead,dbdeaddd)
f90_bdead = f90funclib.__edallommod_MOD_bdead_allom  
  

for ipft in range(numpft):
    print 'Solving for pft: {}'.format(ipft+1)

    ch_min = c_double(pftparms[ipft]['h_min'])
    cd = c_double(-9.0)
    cdddh = c_double(-9.0)
    cipft = c_int(ipft+1)

    # Calculate the d_min parameter
    iret=f90_h2d(byref(c_double(pftparms[ipft]['h_min'])), \
                 byref(cipft),byref(cd),byref(cdddh))
    pftparms[ipft].update({'d_min':cd.value})
    print 'py: h_min of {!r} generated d_min of {!r}' \
        .format(pftparms[ipft]['h_min'],pftparms[ipft]['d_min'])

    # Calculate the d_max parameter
    iret=f90_h2d(byref(c_double(pftparms[ipft]['h_max'])), \
                 byref(cipft),byref(cd),byref(cdddh))
    pftparms[ipft].update({'d_max':cd.value})
    print 'py: h_max of {!r} generated d_max of {!r}' \
        .format(pftparms[ipft]['h_max'],pftparms[ipft]['d_max'])

    # Generate a vector of diameters (use dbh)
    dbh[ipft,:] = np.linspace(pftparms[ipft]['d_min'],maxdbh,num=ndbh)

    # Initialize various output vectors
    cd = c_double(dbh[ipft,0])

    # Output arguments must be initialized
    ch = c_double(-9.0)
    cdhdd = c_double(-9.0)
    cbag = c_double(-9.0)
    cdbagdd = c_double(-9.0)
    cblmax = c_double(-9.0)
    cdblmaxdd = c_double(-9.0)
    cbfrmax = c_double(-9.0)
    cdbfrmaxdd = c_double(-9.0)
    cbcr = c_double(-9.0)
    cdbcrdd = c_double(-9.0)
    cbsap = c_double(-9.0)
    cdbsapdd = c_double(-9.0)
    cbdead = c_double(-9.0)
    cdbdeaddd = c_double(-9.0)

    # Integrated Height
    iret=f90_h(byref(cd),byref(cipft),byref(ch),byref(cdhdd))
    hi[ipft,0] = ch.value
    print 'py: initialize h[{},0]={}'.format(ipft+1,ch.value)

    # Integrated AGB
    iret=f90_bag(byref(cd),byref(ch),byref(cipft),byref(cbag),byref(cdbagdd))
    bagi[ipft,0] = cbag.value
    print 'py: initialize bagi[{},0]={}'.format(ipft+1,cbag.value)

    # Integrated blmax
    iret=f90_blmax(byref(cd),byref(ch),byref(cipft),byref(cblmax),byref(cdblmaxdd))
    blmaxi[ipft,0] = cblmax.value
    print 'py: initialize blmaxi[{},0]={}'.format(ipft+1,cblmax.value)

    # Deterministic blmax
    blmaxd[ipft,0] = cblmax.value

    # bfrmax (d,blmax,dblmaxdd,ipft,bfrmax,dbfrmaxdd)
    iret=f90_bfrmax(byref(cd),byref(cblmax),byref(cdblmaxdd), \
                    byref(cipft),byref(cbfrmax),byref(cdbfrmaxdd))
    bfrmax[ipft,0] = cbfrmax.value
    print 'py: initialize bfrmax[{},0]={}'.format(ipft+1,cbfrmax.value)
    
    # bcr (d,bag,dbagdd,ipft,bcr,dbcrdd)
    iret=f90_bcr(byref(cd),byref(cbag),byref(cdbagdd), \
                 byref(cipft),byref(cbcr),byref(cdbcrdd))
    bcr[ipft,0] = cbcr.value
    print 'py: initialize bcr[{},0]={}'.format(ipft+1,cbcr.value)

    # integrated bsap (d,h,blmax,dblmaxdd,dhdd,ipft,bsap,dbsapdd)
    iret=f90_bsap(byref(cd),byref(ch),byref(cblmax),byref(cdblmaxdd), \
                  byref(cdhdd),byref(cipft),byref(cbsap),byref(cdbsapdd))
    bsapi[ipft,0] = cbsap.value
    print 'py: initialize bsapi[{},0]={}'.format(ipft+1,cbsap.value)

    # Deterministic bsap
    bsapd[ipft,0] = cbsap.value

    # bdead (bag,bcr,blmax,bsap,dbagdd,dbcrdd,dblmaxdd,dbsapdd,bdead,dbdeaddd)
    iret=f90_bdead(byref(cbag),byref(cbcr),byref(cblmax),byref(cbsap), \
                   byref(cdbagdd),byref(cdbcrdd),byref(cdblmaxdd), \
                   byref(cdbsapdd),byref(cbdead),byref(cdbdeaddd))
    bdead[ipft,0] = cbdead.value
    print 'py: initialize bdead[{},0]={}'.format(ipft+1,cbdead.value)

    # Effective dbh
    dbhe[ipft,0] = dbh[ipft,0]

    # the metric that shan't be spoken
    blmax_o_dbagdh[ipft,0]  = blmaxi[ipft,0]/(cdbagdd.value/cdhdd.value)

    for idi in range(1,ndbh):
        
        dp = dbh[ipft,idi-1]  # previous position
        dc = dbh[ipft,idi]    # current position
        dd = dc-dp
        
        cdp = c_double(dp)
        cdc = c_double(dc)
        cdbhe = c_double(-9.0)
        cddedh = c_double(-9.0)

        # integrate height
        iret=f90_h(byref(cdc),byref(cipft),byref(ch),byref(cdhdd))
        hi[ipft,idi] = hi[ipft,idi-1] + cdhdd.value*dd
        print 'dc = {}, dhdd = {}'.format(dc,cdhdd.value)

        # diagnose effective diameter
        iret=f90_h2d(byref(ch),byref(cipft),byref(cdbhe),byref(cddedh))
        dbhe[ipft,idi] = cdbhe.value

        # diagnose AGB
        iret=f90_bag(byref(cdc),byref(c_double(hi[ipft,idi])),byref(cipft),byref(cbag),byref(cdbagdd))
        bagd[ipft,idi] = cbag.value

        # integrate AGB 
        iret=f90_bag(byref(cdp),byref(c_double(hi[ipft,idi-1])),byref(cipft),byref(cbag),byref(cdbagdd))
        bagi[ipft,idi] = bagi[ipft,idi-1] + cdbagdd.value*dd

        # diagnose blmax
        iret=f90_blmax(byref(cdc),byref(c_double(hi[ipft,idi])),byref(cipft),byref(cblmax),byref(cdblmaxdd))
        blmaxd[ipft,idi] = cblmax.value

        # integrate blmax
        iret=f90_blmax(byref(cdp),byref(c_double(hi[ipft,idi-1])),byref(cipft),byref(cblmax),byref(cdblmaxdd))
        blmaxi[ipft,idi] = blmaxi[ipft,idi-1] + cdblmaxdd.value*dd

        # integrate bfrmax
        iret=f90_bfrmax(byref(cdp),byref(cblmax),byref(cdblmaxdd),byref(cipft),byref(cbfrmax),byref(cdbfrmaxdd))
        bfrmax[ipft,idi] = bfrmax[ipft,idi-1] + cdbfrmaxdd.value*dd
        
        # integrate bcr
        iret=f90_bcr(byref(cdp),byref(cbag),byref(cdbagdd),byref(cipft),byref(cbcr),byref(cdbcrdd))
        bcr[ipft,idi] = bcr[ipft,idi-1] + cdbcrdd.value*dd

        # diagnose bsap
        iret=f90_bsap(byref(cdc),byref(c_double(hi[ipft,idi])), \
                      byref(c_double(blmaxd[ipft,idi])),byref(cdblmaxdd), \
                      byref(cdhdd),byref(cipft),byref(cbsap),byref(cdbsapdd))

        bsapd[ipft,idi] = cbsap.value

        # integrate bsap
        iret=f90_bsap(byref(cdp),byref(c_double(hi[ipft,idi-1])), \
                      byref(c_double(blmaxd[ipft,idi-1])),byref(cdblmaxdd), \
                      byref(cdhdd),byref(cipft),byref(cbsap),byref(cdbsapdd))
        bsapi[ipft,idi] = bsapi[ipft,idi-1] + cdbsapdd.value*dd
        
        # the metric that shan't be spoken 
        # previous t-step derivatives are used for simplicity
        if hi[ipft,idi]>=pftparms[ipft]['h_max']:
            blmax_o_dbagdh[ipft,idi] = 0
        else:
            blmax_o_dbagdh[ipft,idi]  = blmaxi[ipft,idi-1]/(cdbagdd.value/cdhdd.value)

        # Diagnose bdead
        iret=f90_bdead(byref(c_double(bagi[ipft,idi])), \
                       byref(c_double(bcr[ipft,idi])), \
                       byref(c_double(blmaxi[ipft,idi])), \
                       byref(c_double(bsapi[ipft,idi])), \
                       byref(cdbagdd),byref(cdbcrdd),byref(cdblmaxdd), \
                       byref(cdbsapdd),byref(cbdead),byref(cdbdeaddd))
        bdead[ipft,idi] = cbdead.value

fig1 = plt.figure()
for ipft in range(numpft):
    plt.plot(dbh[ipft,:],hi[ipft,:],label="pft{}".format(ipft+1))

plt.legend(loc='lower right')
#plt.plot(np.transpose(dbh),np.transpose(hi))
plt.xlabel('diameter [cm]')
plt.ylabel('height [m]')
plt.title('Integrated Heights')
plt.grid(True)
plt.savefig("hi.png")


fig2=plt.figure()
for ipft in range(numpft):
    plt.plot(dbh[ipft,:],blmaxi[ipft,:],label="pft{}".format(ipft+1))
plt.legend(loc='lower right')
#plt.plot(np.transpose(dbh),np.transpose(hi))
plt.xlabel('diameter [cm]')
plt.ylabel('mass [kgC]')
plt.title('Maximum Leaf Biomass')
plt.grid(True)
plt.savefig("blmaxi.png")
plt.show()


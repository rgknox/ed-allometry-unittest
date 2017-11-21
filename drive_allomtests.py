import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mp
import ctypes
from ctypes import * #byref, cdll, c_int, c_double, c_char_p, c_long
import xml.etree.ElementTree as ET

pft_xml_file = "allom_params.xml"
allom_wrap_object = "./EDAllomUnitWrap.o"
allom_lib_object = "./FatesAllometryMod.o"

def wait():
    msvcrt.getch()

# ==============================================================================
# The following lists define what variables are to be expected for the
# different datatypes (double and int).  It is expected that the XML
# (and netcdf later) and the F90 wrapper have the same entries or more
# ==============================================================================

# These are the expected PFT parameters that are double precision

expt_par_dp = [ 'fates_allom_dbh_maxheight','fates_allom_hmode','fates_allom_amode', \
                'fates_allom_lmode','fates_allom_smode','fates_allom_cmode','fates_allom_fmode', \
                'fates_allom_d2h1','fates_allom_d2h2','fates_allom_d2h3','fates_allom_agb1', \
                'fates_allom_agb2','fates_allom_agb3','fates_allom_agb4','fates_allom_d2bl1', \
                'fates_allom_d2bl2','fates_allom_d2bl3','fates_wood_density','fates_c2b', \
                'fates_allom_latosa_int','fates_allom_latosa_slp','fates_slatop','fates_allom_l2fr', \
                'fates_allom_agb_frac','fates_allom_blca_expnt_diff', \
                'fates_allom_d2ca_coefficient_min','fates_allom_d2ca_coefficient_max' ]

# These are the expected PFT parameters that are integers
expt_par_int = []

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
    plist.update({'name':elem.attrib['tag']})
    pftparms.append(plist)

# ==============================================================================
# Load the fortran allometry library using python's ctypes library
# ==============================================================================

f90wraplib = ctypes.CDLL(allom_wrap_object,mode=ctypes.RTLD_GLOBAL)
f90funclib = ctypes.CDLL(allom_lib_object,mode=ctypes.RTLD_GLOBAL)

# ==============================================================================
# Allocate fortran PFT arrays
# ==============================================================================

iret=f90wraplib.__edpftvarcon_MOD_edpftvarconalloc(byref(c_int(numpft)))


# ==============================================================================
# Populate the Fortran PFT structure
# ==============================================================================

# First set the arg types
f90wraplib.__edpftvarcon_MOD_edpftvarconpyset.argtypes = \
    [POINTER(c_int),POINTER(c_double),POINTER(c_int),c_char_p,c_long]

for ipft in range(numpft):
    elem=pftparms[ipft]
    for pname in expt_par_dp:
        print 'py: sending to F90: {0} = {1}'.format(pname,elem[pname])
        iret=f90wraplib.__edpftvarcon_MOD_edpftvarconpyset(c_int(ipft+1), \
                    c_double(elem[pname]),c_int(0),c_char_p(pname),c_long(len(pname)))
    for pname in expt_par_int:
        print 'py: sending to F90: {0} = {1}'.format(pname,elem[pname])
        iret=f90wraplib.__edpftvarcon_MOD_edpftvarconpyset(c_int(ipft+1), \
                    c_double(-999.9),c_int(elem[pname]),c_char_p(pname),len(pname))



# Some testing constants
ndbh = 2000
maxdbh = 200.0
canopy_trim = 1.0
site_spread = 0.0

# =========================================================================
# Generate a vector of diameters that starts at the smallest known diameter
# and extends to 150cm

# =========================================================================
# Initialize Output Arrays

blmaxi  = np.zeros((numpft,ndbh))
blmaxd  = np.zeros((numpft,ndbh))

bfrmax = np.zeros((numpft,ndbh))
hi     = np.zeros((numpft,ndbh)) # Integrated height
hd     = np.zeros((numpft,ndbh)) # Diagnosed height
bagi   = np.zeros((numpft,ndbh))
bagd   = np.zeros((numpft,ndbh))
dbh    = np.zeros((numpft,ndbh))
bcr    = np.zeros((numpft,ndbh))
bsapi  = np.zeros((numpft,ndbh))
bsapd  = np.zeros((numpft,ndbh))
bdead  = np.zeros((numpft,ndbh))
dbhe   = np.zeros((numpft,ndbh))
camin  = np.zeros((numpft,ndbh))
ldense = np.zeros((numpft,ndbh))

blmax_o_dbagdh = np.zeros((numpft,ndbh))
blmax_o_dbagdd = np.zeros((numpft,ndbh))

# Minimum DBH and maximum DBH are diagnosed
# ==============================================================================

f90_h2d       = f90funclib.__fatesallometrymod_MOD_h2d_allom     #(h,ipft,d,dddh)
f90_h         = f90funclib.__fatesallometrymod_MOD_h_allom       #(d,ipft,h,dhdd)
f90_bag       = f90funclib.__fatesallometrymod_MOD_bag_allom     #(d,h,ipft,bag,dbagdd)
f90_bleaf     = f90funclib.__fatesallometrymod_MOD_bleaf         #(d,h,ipft,canopy_trim,bl,dbldd)
f90_bsap      = f90funclib.__fatesallometrymod_MOD_bsap_allom    #(d,h,ipft,canopy_trim,bsap,dbsapdd)
f90_bcr       = f90funclib.__fatesallometrymod_MOD_bcr_allom     #(d,h,ipft,bcr,dbcrdd)
f90_bfineroot = f90funclib.__fatesallometrymod_MOD_bfineroot     #(d,h,ipft,canopy_trim,bfr,dbfrdd)
f90_carea     = f90funclib.__fatesallometrymod_MOD_carea_allom   #(d,nplant,site_spread,ipft,c_area)

#(bag,bcr,bsap,ipft,bdead,dbagdd,dbcrdd,dbsapdd,dbdeaddd)
f90_bdead = f90funclib.__fatesallometrymod_MOD_bdead_allom
  

for ipft in range(numpft):
    print 'py: Solving for pft: {}'.format(ipft+1)

    ch_min = c_double(1.5)
    cd     = c_double(-9.0)
    cdddh  = c_double(-9.0)
    cipft  = c_int(ipft+1)
    cinit  = c_int(0)

    # Calculate the minimum dbh
    iret=f90_h2d(byref(ch_min),byref(cipft),byref(cd),byref(cdddh),byref(cinit))

    pftparms[ipft].update({'d_min':cd.value})
    
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
    ccamin = c_double(-9.0)


    # Initialize Height   #(d,ipft,h,dhdd)
    iret=f90_h(byref(cd),byref(cipft),byref(ch),byref(cdhdd))
    hi[ipft,0] = ch.value
    hd[ipft,0] = ch.value
    print 'py: initialize h[{},0]={}'.format(ipft+1,ch.value)

    # Initialize AGB      #(d,h,ipft,bag,dbagdd)
    iret=f90_bag(byref(cd),byref(ch_min),byref(cipft),byref(cbag),byref(cdbagdd))
    bagi[ipft,0] = cbag.value
    print 'py: initialize bagi[{},0]={}'.format(ipft+1,cbag.value)

    # Initialize bleaf    #(d,h,ipft,canopy_trim,bl,dbldd)
    iret=f90_bleaf(byref(cd),byref(ch_min),byref(cipft),byref(c_double(1.0)),byref(cblmax),byref(cdblmaxdd))
    blmaxi[ipft,0] = cblmax.value
    blmaxd[ipft,0] = cblmax.value
    print 'py: initialize blmaxi[{},0]={}'.format(ipft+1,cblmax.value)

    # calculate crown area (d,nplant,site_spread,ipft,c_area)  Using nplant = 1, generates units of m2
    #  spread is likely 0.0, which is the value it tends towards when canopies close
    iret= f90_carea(byref(cd),byref(c_double(1.0)),byref(c_double(site_spread)),byref(cipft),byref(ccamin))
    camin[ipft,0]  = ccamin.value
    ldense[ipft,0] = blmaxi[ipft,0]/camin[ipft,0]

    # Initialize fine roots  #(d,h,ipft,canopy_trim,bfr,dbfrdd)
    iret=f90_bfineroot(byref(cd),byref(ch_min),byref(cipft),byref(c_double(1.0)), \
                       byref(cbfrmax),byref(cdbfrmaxdd))
    bfrmax[ipft,0] = cbfrmax.value
    print 'py: initialize bfrmax[{},0]={}'.format(ipft+1,cbfrmax.value)
    
    # Initialize coarse roots #(d,h,ipft,bcr,dbcrdd)
    iret=f90_bcr(byref(cd),byref(ch_min),byref(cipft), \
                 byref(cbcr),byref(cdbcrdd))
    bcr[ipft,0] = cbcr.value
    print 'py: initialize bcr[{},0]={}'.format(ipft+1,cbcr.value)


    # Initialize bsap  #(d,h,ipft,canopy_trim,bsap,dbsapdd)
    iret=f90_bsap(byref(cd),byref(ch_min),byref(cipft),byref(c_double(1.0)),byref(cbsap),byref(cdbsapdd))
    bsapi[ipft,0] = cbsap.value
    bsapd[ipft,0] = cbsap.value
    print 'py: initialize bsapi[{},0]={}'.format(ipft+1,cbsap.value)
    
    # bdead #(bag,bcr,bsap,ipft,bdead,dbagdd,dbcrdd,dbsapdd,dbdeaddd)
    iret=f90_bdead(byref(cbag),byref(cbcr),byref(cbsap),byref(cipft), \
                   byref(cbdead),byref(cdbagdd),byref(cdbcrdd), \
                   byref(cdbsapdd),byref(cdbdeaddd))
    bdead[ipft,0] = cbdead.value
    print 'py: initialize bdead[{},0]={}'.format(ipft+1,cbdead.value)

    # the metric that shan't be spoken
    blmax_o_dbagdh[ipft,0]  = blmaxi[ipft,0]/(cdbagdd.value/cdhdd.value)

    # the metric that shan't be spoken
    blmax_o_dbagdd[ipft,0]  = blmaxi[ipft,0]/(cdbagdd.value)

    for idi in range(1,ndbh):
        
        dp = dbh[ipft,idi-1]  # previous position
        dc = dbh[ipft,idi]    # current position
        dd = dc-dp
        
        cdp = c_double(dp)
        cdc = c_double(dc)
        cdbhe = c_double(-9.0)
        cddedh = c_double(-9.0)

        # integrate height  #(d,ipft,h,dhdd)
        iret=f90_h(byref(cdc),byref(cipft),byref(ch),byref(cdhdd))
        hi[ipft,idi] = hi[ipft,idi-1] + cdhdd.value*dd

        # diagnosed height
        hd[ipft,idi] = ch.value

        # diagnose AGB  #(d,h,ipft,bag,dbagdd)
        iret=f90_bag(byref(cdc),byref(c_double(hi[ipft,idi])),byref(cipft),byref(cbag),byref(cdbagdd))
        bagd[ipft,idi] = cbag.value

        # integrate AGB #(d,h,ipft,bag,dbagdd)
        iret=f90_bag(byref(cdp),byref(c_double(hi[ipft,idi-1])),byref(cipft),byref(cbag),byref(cdbagdd))
        bagi[ipft,idi] = bagi[ipft,idi-1] + cdbagdd.value*dd

        # diagnose bleaf #(d,h,ipft,blmax,dblmaxdd)
        iret=f90_bleaf(byref(cdc),byref(c_double(hi[ipft,idi])),byref(cipft),byref(c_double(1.0)),byref(cblmax),byref(cdblmaxdd))
        blmaxd[ipft,idi] = cblmax.value

        # integrate bleaf #(d,h,ipft,blmax,dblmaxdd)
        iret=f90_bleaf(byref(cdp),byref(c_double(hi[ipft,idi-1])),byref(cipft),byref(c_double(1.0)),byref(cblmax),byref(cdblmaxdd))
        blmaxi[ipft,idi] = blmaxi[ipft,idi-1] + cdblmaxdd.value*dd

        iret=f90_bleaf(byref(cdp),byref(c_double(hi[ipft,idi-1])),byref(cipft),byref(c_double(1.0)),byref(cblmax),byref(cdblmaxdd))

        # calculate crown area (d,nplant,site_spread,ipft,c_area)  Using nplant = 1, generates units of m2
        iret= f90_carea(byref(cdc),byref(c_double(1.0)),byref(c_double(site_spread)),byref(cipft),byref(ccamin))
        camin[ipft,idi]  = ccamin.value

        # leaf mass per square meter of crown
        ldense[ipft,idi] = blmaxd[ipft,idi]/camin[ipft,idi]

        # integrate bfineroot #(d,h,ipft,canopy_trim,bfr,dbfrdd)
        iret=f90_bfineroot(byref(cdp),byref(c_double(hi[ipft,idi-1])),byref(cipft),byref(c_double(1.0)),byref(cbfrmax),byref(cdbfrmaxdd))
        bfrmax[ipft,idi] = bfrmax[ipft,idi-1] + cdbfrmaxdd.value*dd
        
        # integrate bcr #(d,h,ipft,bcr,dbcrdd)
        iret=f90_bcr(byref(cdp),byref(c_double(hi[ipft,idi-1])),byref(cipft),byref(cbcr),byref(cdbcrdd))
        bcr[ipft,idi] = bcr[ipft,idi-1] + cdbcrdd.value*dd

        # diagnose bsap  #(d,h,ipft,canopy_trim,bsap,dbsapdd)
        iret=f90_bsap(byref(cdc),byref(c_double(hi[ipft,idi])),byref(cipft),byref(c_double(1.0)),byref(cbsap),byref(cdbsapdd))
        bsapd[ipft,idi] = cbsap.value

        # integrate bsap
        iret=f90_bsap(byref(cdp),byref(c_double(hi[ipft,idi-1])),byref(cipft),byref(c_double(1.0)),byref(cbsap),byref(cdbsapdd))
        bsapi[ipft,idi] = bsapi[ipft,idi-1] + cdbsapdd.value*dd
        
        # the metric that shan't be spoken 
        # previous t-step derivatives are used for simplicity
        if cdhdd.value<0.000001:
            blmax_o_dbagdh[ipft,idi] = None
        else:
            blmax_o_dbagdh[ipft,idi]  = blmaxi[ipft,idi-1]/(cdbagdd.value/cdhdd.value)

        # the metric that shan't be spoken 
        # previous t-step derivatives are used for simplicity
        blmax_o_dbagdd[ipft,idi]  = blmaxi[ipft,idi-1]/(cdbagdd.value)


        # Diagnose bdead (bag,bcr,bsap,ipft,bdead,dbagdd,dbcrdd,dbsapdd,dbdeaddd)

        iret=f90_bdead(byref(c_double(bagi[ipft,idi])), \
                       byref(c_double(bcr[ipft,idi])), \
                       byref(c_double(bsapi[ipft,idi])), \
                       byref(cipft), byref(cbdead), \
                       byref(cdbagdd),byref(cdbcrdd), \
                       byref(cdbsapdd),byref(cdbdeaddd))
        bdead[ipft,idi] = cbdead.value


linestyles = ['-', '--', '-.', '-', '--', '-.', '-']

my_colors  = ['black','blue','chocolate','darkgrey','darkolivegreen','darkmagenta','lightslategrey']

lwidths    = [1.5, 2.0 , 2.0, 1.5, 2.0, 2.0, 1.5]

#font = {'family' : 'normal',
#        'weight' : 'normal',
#        'size'   : 16}

#mp.rc('font', **font)
mp.rcParams.update({'font.size': 16})
mp.rcParams["savefig.directory"] = ""    #os.chdir(os.path.dirname(__file__))

legfs = 14
lwidth = 2.0

fig1 = plt.figure()
for ipft in range(numpft):
    plt.plot(dbh[ipft,:],hi[ipft,:],linestyle=linestyles[ipft],label="{}".format(pftparms[ipft]['name']),color=my_colors[ipft],linewidth=lwidths[ipft])
plt.legend(loc='lower right',fontsize=legfs)
#plt.plot(np.transpose(dbh),np.transpose(hi))
plt.xlabel('diameter [cm]')
plt.ylabel('height [m]')
plt.title('Integrated Heights')
plt.grid(True)
plt.savefig("plots/hi.png")



fig1_0 = plt.figure()
for ipft in range(numpft):
    plt.plot(dbh[ipft,0:15],hi[ipft,0:15],linestyle=linestyles[ipft],label="{}".format(pftparms[ipft]['name']),color=my_colors[ipft],linewidth=lwidths[ipft])
plt.legend(loc='lower right',fontsize=legfs)
#plt.plot(np.transpose(dbh),np.transpose(hi))
plt.xlabel('diameter [cm]')
plt.ylabel('height [m]')
plt.title('Integrated Heights')
plt.grid(True)
plt.savefig("plots/hi.png")

fig1_1 = plt.figure()
for ipft in range(numpft):
    plt.plot(hd[ipft,:],hi[ipft,:],linestyle=linestyles[ipft],label="{}".format(pftparms[ipft]['name']),color=my_colors[ipft],linewidth=lwidths[ipft])
plt.legend(loc='lower right',fontsize=legfs)
#plt.plot(np.transpose(dbh),np.transpose(hi))
plt.xlabel('height (diagnosed) [m]')
plt.ylabel('height (integrated) [m]')
plt.title('Height')
plt.grid(True)
plt.savefig("plots/hdhi.png")

fig2=plt.figure()
for ipft in range(numpft):
    plt.plot(blmaxd[ipft,:],blmaxi[ipft,:],linestyle=linestyles[ipft],label="{}".format(pftparms[ipft]['name']),color=my_colors[ipft],linewidth=lwidths[ipft])
plt.legend(loc='lower right',fontsize=legfs)
#plt.plot(np.transpose(dbh),np.transpose(hi))
plt.xlabel('diagnosed [kgC]')
plt.ylabel('integrated [kgC]')
plt.title('Maximum Leaf Biomass')
plt.grid(True)
plt.savefig("plots/blmaxdi.png")

fig3=plt.figure()
for ipft in range(numpft):
    plt.plot(dbh[ipft,:],blmaxi[ipft,:],linestyle=linestyles[ipft],label="{}".format(pftparms[ipft]['name']),color=my_colors[ipft],linewidth=lwidths[ipft])
plt.legend(loc='upper left',fontsize=legfs)
#plt.plot(np.transpose(dbh),np.transpose(hi))
plt.xlabel('diameter [cm]')
plt.ylabel('mass [kgC]')
plt.title('Maximum Leaf Biomass')
plt.grid(True)
plt.savefig("plots/blmaxi.png")

fig3_1=plt.figure()
for ipft in range(numpft):
    plt.semilogy(dbh[ipft,1:15],blmaxi[ipft,1:15],linestyle=linestyles[ipft],label="{}".format(pftparms[ipft]['name']),color=my_colors[ipft],linewidth=lwidths[ipft])
plt.legend(loc='upper left',fontsize=legfs)
#plt.ax.set_yscale('log')
#plt.plot(np.transpose(dbh),np.transpose(hi))
plt.xlabel('diameter [cm]')
plt.ylabel('mass [kgC]')
plt.title('Maximum Leaf Biomass')
plt.grid(True)
plt.savefig("plots/blmaxi_small.png")


fig4=plt.figure()
for ipft in range(numpft):
    plt.plot(dbh[ipft,:],camin[ipft,:],linestyle=linestyles[ipft],label="{}".format(pftparms[ipft]['name']),color=my_colors[ipft],linewidth=lwidths[ipft])
plt.legend(loc='upper left',fontsize=legfs)
plt.xlabel('diameter [cm]')
plt.ylabel('[m2] (closed canopy)')
plt.title('Crown Area')
plt.grid(True)
plt.savefig("plots/carea.png")

fig4_1=plt.figure()
for ipft in range(numpft):
    plt.plot(dbh[ipft,:],ldense[ipft,:],linestyle=linestyles[ipft],label="{}".format(pftparms[ipft]['name']),color=my_colors[ipft],linewidth=lwidths[ipft])
plt.legend(loc='upper left',fontsize=legfs)
plt.xlabel('diameter [cm]')
plt.ylabel('[kgC/m2] (closed canopy)')
plt.title('Leaf Mass Per Crown Area')
plt.grid(True)
plt.savefig("plots/ldense.png")



#fig4=plt.figure()
#ax1 = fig4.add_subplot(1,2,1)
#ax2 = fig4.add_subplot(1,2,2)
#ax2.set_yscale('log')
#for ipft in range(numpft):
#    ax1.plot(dbh[ipft,:],bagi[ipft,:],label="{}".format(pftparms[ipft]['name']))
#    ax2.plot(dbh[ipft,:],bagi[ipft,:],label="{}".format(pftparms[ipft]['name']))
#plt.legend(loc='upper left')
#plt.xlabel('diameter [cm]')
#plt.ylabel('mass [kgC]')
#plt.title('Above Ground Biomass')
#plt.grid(True)
#plt.savefig("plots/agbi.png")

fig6=plt.figure()
for ipft in range(numpft):
    plt.plot(dbh[ipft,:],bagi[ipft,:]/1000,linestyle=linestyles[ipft],label="{}".format(pftparms[ipft]['name']),color=my_colors[ipft],linewidth=lwidths[ipft])
plt.legend(loc='upper left',fontsize=legfs)
plt.xlabel('diameter [cm]')
plt.ylabel('AGB [MgC]')
plt.title('Above Ground Biomass')
plt.grid(True)
plt.savefig("plots/agbi.png")

fig5=plt.figure()
for ipft in range(numpft):
    gpmask  = np.isfinite(blmax_o_dbagdh[ipft,:])
    plt.plot(dbh[ipft,gpmask],blmax_o_dbagdh[ipft,gpmask],linestyle=linestyles[ipft],label="{}".format(pftparms[ipft]['name']),color=my_colors[ipft],linewidth=lwidths[ipft])
plt.legend(loc='upper right',fontsize=legfs)
plt.xlabel('diameter [cm]')
plt.ylabel('growth potential: bl/(dAGB/dh) [m]')
plt.title('Height Growth Potential')
plt.grid(True)
plt.savefig("plots/gpot_h.png")

fig6=plt.figure()
for ipft in range(numpft):
    plt.plot(dbh[ipft,:],blmax_o_dbagdd[ipft,:],linestyle=linestyles[ipft],label="{}".format(pftparms[ipft]['name']),color=my_colors[ipft],linewidth=lwidths[ipft])
plt.legend(loc='upper left',fontsize=legfs)
plt.xlabel('diameter [cm]')
plt.ylabel('growth potential: bl/(dAGB/dd) [cm]')
plt.title('Diameter Growth Potential')
plt.grid(True)
plt.savefig("plots/gpot_d.png")


plt.show()


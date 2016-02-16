import numpy
from ctypes import * #byref, cdll, c_int, c_double, c_char_p, c_long
import xml.etree.ElementTree as ET

pft_xml_file = "allom_params.xml"
allom_fortran_archive = "./EDAllomUnitWrap.o"


# ====================================================================
# 1. Load and parse the PFT xml tree
# ====================================================================

pftroot = ET.parse(pft_xml_file).getroot()

# Get the names of the non-pft delineated parameters

# These are the expected PFT parameters that are double precision
expt_par_dp = ['c2b','eclim','llspan','bl_min','h_max','h_min','slatop', \
               'd_adult','d_sap','l2f_ratio','agb_fraction','latosa_int', \
               'latosa_slp','d2h1_ad','d2h2_ad','d2h3_ad','d2h1_sap', \
               'd2h2_sap','d2bl1_ad','d2bl2_ad','d2bl3_ad','d2bl1_sap', \
               'd2bl2_sap','d2bag1','d2bag2','wood_density']


pftparms = []

numpft = 0
for elem in pftroot.iter('pft'):
    print(elem.tag)
    plist = {}
    for iv in expt_par_dp:
        pftelem = elem.find(iv)
        if pftelem!=None:
            plist[pftelem.tag] = float(pftelem.text)
    numpft += 1
    pftparms.append(plist)


print numpft


for elem in pftparms:
    dir(elem)
    for parm in elem:
        print parm,elem[parm]


quit()


# These are the expected PFT parameters that are integers
expt_par_int = ['hallom_mode','lallom_mode','fallom_mode','aallom_mode', \
                'callom_mode','sallom_mode']


pftparms = []

pftparms.append()






quit()
        


#for iv = 1:numel(expt_par_dp)
#    % Check to see if this is a cross-pft parameter designation
#    id = find(strcmp(pnames,expt_par_dp{iv}));
#    if(isempty(id))
#        % Turns out this is probably a pft specific parameter
#        for ip=1:n_pfts
#            id = find(strcmp(expt_par_dp{iv},fieldnames(xmlroot.all.pfts.pft{ip})),1);
#            if(isempty(id))
#                display('Expected parameter cannot be found in XML');
#                display(sprintf('%s, pft %d',expt_par_dp{iv},ip));
#                return;
#            else
#                pftcon.(expt_par_dp{iv})(ip) = str2double(xmlroot.all.pfts.pft{ip}.(expt_par_dp{iv}).Text)#;
#            end
#        end
#    else
#        for ip=1:n_pfts
#            pftcon.(expt_par_dp{iv})(ip) = str2double(xmlroot.all.(pnames{id}).Text);
#        end
#    end
#end




# edallomunitwrap_MOD_edecophysconalloc

fallomlib = cdll.LoadLibrary(allom_fortran_archive)




numpft = c_int(10)

fallomlib.__edallomunitwrap_MOD_edecophysconalloc(byref(numpft))

#pname = c_char_p("bl_min")
#plen  = c_long(80)
ipft  = c_int(2)
rval  = c_double(0.04)
ival  = c_int(10)


pname = "bl_min"
plen = len(pname)

fallomlib.__edallomunitwrap_MOD_edecophysconpyset.argtypes = [POINTER(c_int),POINTER(c_double),POINTER(c_int),c_char_p,c_long]
fallomlib.__edallomunitwrap_MOD_edecophysconpyset(ipft,rval,ival,pname,plen)

#fallomlib.__edallomunitwrap_MOD_edecophysconpyset(byref(pname),byref(ipft),byref(rval),byref(ival))



#fadd = cdll.LoadLibrary('./add.so')

#arg1 = c_double(2.0)
#arg2 = c_double(4.0)

#arg3 = c_int(2)
#arg4 = c_int(4)

#fadd.__add_MOD_initpatch(byref(arg4))

#print fadd.__add_MOD_frand(byref(arg4),byref(arg3))

#print fadd.__add_MOD_fiadd(byref(arg3),byref(arg4))

#print fadd.__add_MOD_fdadd(byref(arg1),byref(arg2))


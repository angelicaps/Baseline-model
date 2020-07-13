

# define:
#       -  # of lightcurves  (2 to start)
#       -  model
#       - jump parameters
#                     - dF
#                     - b
#                     - P
#                     - T0
#                     - a/Rs


# procedure:
#       - set up all modules
#       - read the data
#       - read the input jump parameters and set their initial values
#       - call the mcmc
#       - look at the results

# model:
#       - call a transit lightcurve for each data set
#       - multiply this lightcurve with the baseline
#       - output it as one large array
import itertools
from outputs_v6_bat import *
from jitter_v1 import *
from ecc_om_par import *
from groupweights_v2 import *
from corfac_v8 import *
from plots_v12 import *
import fitfunc_v25_bat2 as ff
from basecoeff_v13 import *
from occultnl import *
from occultquad import *
import mc3
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import sys
import csv

#sys.path.append("/home/bella/Documents/Geneva/CONAN/WASP52/mc3/") 
sys.path.append("/home/bella/Documents/Geneva/CONAN/WASP52/conan/")
sys.path.append("/home/bella/Documents/Geneva/CONAN/WASP52/mandagol/")


# matplotlib.use('TKagg')
matplotlib.use('Agg')


#import MCcubed as mc3

time = [0,1,2,3,4]
am = [0,1,2]
xy = [0,1,2]
fw = [0,1,2]
sk = [0,1,2]

combinations = list(itertools.product(time, am, xy, fw, sk))
rang = len(combinations)
print ('Number of combinations', rang)

outF = open("myBICfile.csv", "w")
for i in range(rang):
    # plt.ion()

    # read the input file
    file = open('/home/bella/Documents/Geneva/CONAN/WASP52/input_v26_BAT.dat')
    dump = file.readline()
    dump = file.readline()
    fpath = file.readline()         # the path where the files are
    fpath = fpath.rstrip()
    dump = file.readline()
    dump = file.readline()

    # now read in the lightcurves one after the other
    names = []                    # array where the LC filenames are supposed to go
    filters = []                  # array where the filter names are supposed to go
    lamdas = []
    bases = []                    # array where the baseline exponents are supposed to go
    groups = []                   # array where the group indices are supposed to go
    grbases = []
    dump = file.readline()        # read the next line
    dump = file.readline()        # read the next line

    while dump[0] != '#':           # if it is not starting with # then
        adump = dump.split()          # split it
        names.append(adump[0])      # append the first field to the name array
        # append the second field to the filters array
        filters.append(adump[1])
        # append the second field to the filters array
        lamdas.append(adump[2])
        # string array of the baseline function exponents
        #strbase = adump[3:10]
        time = combinations[i][0]
        am = combinations[i][1]
        xy = combinations[i][2]
        fw = combinations[i][3]
        sk = combinations[i][4]
        #strbase = [time, '0', '0', '0', '1', '0', '0']
        #base = [int(i) for i in strbase]
        base = [time, am, xy, fw, sk, 0, 0]
        print ('base',base)
        bases.append(base)
        print ('bases',bases)
        group = int(adump[10])
        print ('group',group)
        groups.append(group)
        print ('groups',groups)
        grbase = int(adump[9])
        print ('grbase',grbase)
        grbases.append(grbase)
        print ('grbases',grbases)
        # and read the following line - either a lc or something starting with "#"
        dump = file.readline()

    nphot = len(names)             # the number of photometry input files
    njumpphot = np.zeros(nphot)
    filnames = np.array(list(sorted(set(filters), key=filters.index)))
    ulamdas = np.array(list(sorted(set(lamdas), key=lamdas.index)))
    grnames = np.array(list(sorted(set(groups))))
    nfilt = len(filnames)
    ngroup = len(grnames)
    dump = file.readline()

    RVnames = []
    RVbases = []
    gammas = []
    gamsteps = []
    gampri = []
    gamprilo = []
    gamprihi = []

    dump = file.readline()
    dump = file.readline()
    while dump[0] != '#':           # if it is not starting with # then
        adump = dump.split()
        # append the first field to the RVname array
        RVnames.append(adump[0])
        # string array of the baseline function exponents
        strbase = adump[1:5]
        base = [int(i) for i in strbase]
        RVbases.append(base)
        gammas.append(adump[5])
        gamsteps.append(adump[6])
        gampri.append(adump[8])
        gampriloa = (0. if (adump[7] == 'n' or adump[6] == 0.) else adump[9])
        gamprilo.append(gampriloa)
        gamprihia = (0. if (adump[7] == 'n' or adump[6] == 0.) else adump[10])
        gamprihi.append(gamprihia)
        # and read the following line - either a lc or something starting with "#"
        dump = file.readline()

    nRV = len(RVnames)             # the number of RV input files
    njumpRV = np.zeros(nRV)

    # --> adding RV time series' to groups -->> probably not the best way of doing it
    #grnames_add = np.max(grnames)+np.zeros(nRV)+1
    #grnames = np.append(grnames, grnames_add)
    #groups = np.append(groups,grnames_add)
    # ngroup=ngroup+nRV

    gamma_in = np.zeros((nRV, 7))  # set up array to contain the gamma values

    dump = file.readline()

    dump = file.readline()
    adump = dump.split()
    adump[3] = (0. if adump[1] == 'n' else adump[3])
    adump[7] = (0. if (adump[6] == 'n' or adump[1] ==
                       'n' or float(adump[3]) == 0.) else adump[7])
    adump[8] = (0. if (adump[6] == 'n' or adump[1] ==
                       'n' or float(adump[3]) == 0.) else adump[8])
    adump[9] = (0. if (adump[6] == 'n' or adump[1] ==
                       'n' or float(adump[3]) == 0.) else adump[9])
    rprs_in = [float(adump[2]), float(adump[3]), float(adump[4]), float(
        adump[5]), float(adump[7]), float(adump[8]), float(adump[9])]
    rprs0 = np.copy(rprs_in[0])
    if rprs_in[1] != 0.:
        njumpphot = njumpphot+1

    # note: erprs0 is set to zero here! because we don't actually need it any more
    erprs0 = np.copy(0)
    dump = file.readline()
    adump = dump.split()
    adump[3] = (0. if adump[1] == 'n' else adump[3])
    adump[7] = (0. if (adump[6] == 'n' or adump[1] ==
                       'n' or float(adump[3]) == 0.) else adump[7])
    adump[8] = (0. if (adump[6] == 'n' or adump[1] ==
                       'n' or float(adump[3]) == 0.) else adump[8])
    adump[9] = (0. if (adump[6] == 'n' or adump[1] ==
                       'n' or float(adump[3]) == 0.) else adump[9])
    inc_in = [float(adump[2]), float(adump[3]), float(adump[4]), float(
        adump[5]), float(adump[7]), float(adump[8]), float(adump[9])]

    if inc_in[1] != 0.:
        njumpphot = njumpphot+1

    dump = file.readline()
    adump = dump.split()
    adump[3] = (0. if adump[1] == 'n' else adump[3])
    adump[7] = (0. if (adump[6] == 'n' or adump[1] ==
                       'n' or float(adump[3]) == 0.) else adump[7])
    adump[8] = (0. if (adump[6] == 'n' or adump[1] ==
                       'n' or float(adump[3]) == 0.) else adump[8])
    adump[9] = (0. if (adump[6] == 'n' or adump[1] ==
                       'n' or float(adump[3]) == 0.) else adump[9])
    aRs_in = [float(adump[2]), float(adump[3]), float(adump[4]), float(
        adump[5]), float(adump[7]), float(adump[8]), float(adump[9])]
    if aRs_in[1] != 0.:
        njumpphot = njumpphot+1

    dump = file.readline()
    adump = dump.split()
    adump[3] = (0. if adump[1] == 'n' else adump[3])
    adump[7] = (0. if (adump[6] == 'n' or adump[1] ==
                       'n' or float(adump[3]) == 0.) else adump[7])
    adump[8] = (0. if (adump[6] == 'n' or adump[1] ==
                       'n' or float(adump[3]) == 0.) else adump[8])
    adump[9] = (0. if (adump[6] == 'n' or adump[1] ==
                       'n' or float(adump[3]) == 0.) else adump[9])
    T0_in = [float(adump[2]), float(adump[3]), float(adump[4]), float(
        adump[5]), float(adump[7]), float(adump[8]), float(adump[9])]
    if T0_in[1] != 0.:
        njumpRV = njumpRV+1
        njumpphot = njumpphot+1

    dump = file.readline()
    adump = dump.split()
    adump[3] = (0. if adump[1] == 'n' else adump[3])
    adump[7] = (0. if (adump[6] == 'n' or adump[1] ==
                       'n' or float(adump[3]) == 0.) else adump[7])
    adump[8] = (0. if (adump[6] == 'n' or adump[1] ==
                       'n' or float(adump[3]) == 0.) else adump[8])
    adump[9] = (0. if (adump[6] == 'n' or adump[1] ==
                       'n' or float(adump[3]) == 0.) else adump[9])
    per_in = [float(adump[2]), float(adump[3]), float(adump[4]), float(
        adump[5]), float(adump[7]), float(adump[8]), float(adump[9])]
    if per_in[1] != 0.:
        njumpRV = njumpRV+1
        njumpphot = njumpphot+1

    dump = file.readline()
    adump = dump.split()
    adump[3] = (0. if adump[1] == 'n' else adump[3])
    adump[7] = (0. if (adump[6] == 'n' or adump[1] ==
                       'n' or float(adump[3]) == 0.) else adump[7])
    adump[8] = (0. if (adump[6] == 'n' or adump[1] ==
                       'n' or float(adump[3]) == 0.) else adump[8])
    adump[9] = (0. if (adump[6] == 'n' or adump[1] ==
                       'n' or float(adump[3]) == 0.) else adump[9])
    eccpri = np.copy(adump[6])
    ecc_in = [float(adump[2]), float(adump[3]), float(adump[4]), float(
        adump[5]), float(adump[7]), float(adump[8]), float(adump[9])]
    if ecc_in[1] != 0.:
        njumpRV = njumpRV+1

    dump = file.readline()
    adump = dump.split()
    adump[3] = (0. if adump[1] == 'n' else adump[3])
    adump[7] = (0. if (adump[6] == 'n' or adump[1] ==
                       'n' or float(adump[3]) == 0.) else adump[7])
    adump[8] = (0. if (adump[6] == 'n' or adump[1] ==
                       'n' or float(adump[3]) == 0.) else adump[8])
    adump[9] = (0. if (adump[6] == 'n' or adump[1] ==
                       'n' or float(adump[3]) == 0.) else adump[9])
    opri = np.copy(adump[6])
    omega_in = [float(adump[2]), float(adump[3]), float(adump[4]), float(
        adump[5]), float(adump[7]), float(adump[8]), float(adump[9])]
    omega_in = np.multiply(omega_in, np.pi)/180.
    if omega_in[1] != 0.:
        njumpRV = njumpRV+1

    dump = file.readline()
    adump = dump.split()
    adump[3] = (0. if adump[1] == 'n' else adump[3])
    adump[7] = (0. if (adump[6] == 'n' or adump[1] ==
                       'n' or float(adump[3]) == 0.) else adump[7])
    adump[8] = (0. if (adump[6] == 'n' or adump[1] ==
                       'n' or float(adump[3]) == 0.) else adump[8])
    adump[9] = (0. if (adump[6] == 'n' or adump[1] ==
                       'n' or float(adump[3]) == 0.) else adump[9])
    Kpri = np.copy(adump[6])
    K_in = [float(adump[2]), float(adump[3]), float(adump[4]), float(
        adump[5]), float(adump[7]), float(adump[8]), float(adump[9])]
    K_in = np.divide(K_in, 1000.)  # convert to km/s
    if K_in[1] != 0.:
        njumpRV = njumpRV+1

    for i in range(nRV):
        gamma_in[i, :] = [float(gammas[i]), float(
            gamsteps[i]), -1000, 1000, float(gampri[i]), float(gamprilo[i]), float(gamprihi[i])]
        if (gamma_in[i, 1] != 0.):
            njumpRV[i] = njumpRV[i]+1

    # adapt the eccentricity and omega jump parameters sqrt(e)*sin(o), sqrt(e)*cos(o)

    if ((eccpri == 'y' and opri == 'n') or (eccpri == 'n' and opri == 'y')):
        print('priors on eccentricity and omega: either both on or both off')
        print(nothing)

    eos_in, eoc_in = ecc_om_par(ecc_in, omega_in)

    dump = file.readline()
    dump = file.readline()
    dump = file.readline()
    adump = dump.split()

    ddfYN = adump[0]      # (y/n) fit ddFs?

    drprs_op = [0., float(adump[1]), float(adump[2]), float(
        adump[3]), 0., float(adump[5]), float(adump[6])]  # the dRpRs options
    drprs_op[5] = (0 if (adump[4] == 'n' or ddfYN ==
                         'n' or float(adump[2]) == 0.) else adump[5])
    drprs_op[6] = (0 if (adump[4] == 'n' or ddfYN ==
                         'n' or float(adump[2]) == 0.) else adump[6])
    divwhite = adump[7]      # (y/n) divide-white?
    dump = file.readline()

    grprs = np.zeros(ngroup)   # the group rprs values
    egrprs = np.zeros(ngroup)  # the uncertainties of the group rprs values
    dwfiles = []             # the filenames of the white residuals

    for i in range(ngroup):
        dump = file.readline()
        adump = dump.split()
        grprs[i] = np.copy(float(adump[1]))
        egrprs[i] = np.copy(float(adump[2]))
        dwfiles.append(adump[3])

    dwCNMarr = np.array([])      # initializing array with all the dwCNM values
    # initializing array with the indices of each group's dwCNM values
    dwCNMind = []
    dwind = np.array([])
    if (divwhite == 'y'):           # do we do a divide-white?
        for i in range(ngroup):   # read fixed dwCNMs for each group
            tdwCNM, dwCNM = np.loadtxt(
                fpath+dwfiles[i], usecols=(0, 1), unpack=True)
            dwCNMarr = np.concatenate((dwCNMarr, dwCNM), axis=0)
            dwind = np.concatenate(
                (dwind, np.zeros(len(dwCNM), dtype=np.int)+i), axis=0)
            indices = np.where(dwind == i)
            dwCNMind.append(indices)
        if (ddfYN == 'n'):
            print('you can not do divide-white and not fit ddfs!')
            print(nothing)

        for i in range(nphot):
            if (bases[i][6] > 0):
                print('you can not have CNMs active and do divide-white')
                print(nothing)

    if (ddfYN == 'n' and np.max(grbases) > 0):
        print('no dDFs but groups? Not a good idea!')
        print(base)
        print(nothing)

    dump = file.readline()
    dump = file.readline()
    # setup the arrays for the LD coefficients
    c1_in = np.zeros((nfilt, 7))
    c2_in = np.zeros((nfilt, 7))
    c3_in = np.zeros((nfilt, 7))
    c4_in = np.zeros((nfilt, 7))

    for i in range(nfilt):
        dump = file.readline()
        adump = dump.split()
        # make sure the sequence in this array is the same as in the "filnames" array
        j = np.where(filnames == adump[0])
        k = np.where(np.array(filters) == adump[0])

        c1_in[j, :] = [float(adump[2]), float(adump[3]), -3., 3., float(adump[2]),
                       float(adump[4]), float(adump[5])]  # the limits are -3 and 3 => very safe
        c1_in[j, 5] = (0. if (adump[1] == 'n' or float(
            adump[3]) == 0.) else c1_in[j, 5])
        c1_in[j, 6] = (0. if (adump[1] == 'n' or float(
            adump[3]) == 0.) else c1_in[j, 6])
        if c1_in[j, 1] != 0.:
            njumpphot[k] = njumpphot[k]+1

        c2_in[j, :] = [float(adump[6]), float(
            adump[7]), -3., 3., float(adump[6]), float(adump[8]), float(adump[9])]
        c2_in[j, 5] = (0. if (adump[1] == 'n' or float(
            adump[7]) == 0.) else c2_in[j, 5])
        c2_in[j, 6] = (0. if (adump[1] == 'n' or float(
            adump[7]) == 0.) else c2_in[j, 6])
        if c2_in[j, 1] != 0.:
            njumpphot[k] = njumpphot[k]+1

        c3_in[j, :] = [float(adump[10]), float(
            adump[11]), -3., 3., float(adump[10]), float(adump[12]), float(adump[13])]
        c3_in[j, 5] = (0. if (adump[1] == 'n' or float(
            adump[11]) == 0.) else c3_in[j, 5])
        c3_in[j, 6] = (0. if (adump[1] == 'n' or float(
            adump[11]) == 0.) else c3_in[j, 6])
        if c3_in[j, 1] != 0.:
            njumpphot[k] = njumpphot[k]+1

        c4_in[j, :] = [float(adump[14]), float(
            adump[15]), -3., 3., float(adump[14]), float(adump[16]), float(adump[17])]
        c4_in[j, 5] = (0. if (adump[1] == 'n' or float(
            adump[15]) == 0.) else c4_in[j, 5])
        c4_in[j, 6] = (0. if (adump[1] == 'n' or float(
            adump[15]) == 0.) else c4_in[j, 6])
        if c4_in[j, 1] != 0.:
            njumpphot[k] = njumpphot[k]+1

        # convert the input u1, u2 to c1, c2 if the limb-darkening law is quadratic
        if (c3_in[j, 0] == 0. and c4_in[j, 0] == 0 and c3_in[j, 1] == 0. and c4_in[j, 1] == 0.):
            print('Limb-darkening law: quadratic')
            v1 = 2.*c1_in[j, 0]+c2_in[j, 0]
            v2 = c1_in[j, 0]-c2_in[j, 0]
            ev1 = np.sqrt(4.*c1_in[j, 1]**2+c2_in[j, 1]**2)
            ev2 = np.sqrt(c1_in[j, 1]**2+c2_in[j, 1]**2)
            lov1 = np.sqrt(4.*c1_in[j, 5]**2+c2_in[j, 5]**2)
            lov2 = np.sqrt(c1_in[j, 5]**2+c2_in[j, 5]**2)
            hiv1 = np.sqrt(4.*c1_in[j, 6]**2+c2_in[j, 6]**2)
            hiv2 = np.sqrt(c1_in[j, 6]**2+c2_in[j, 6]**2)
            c1_in[j, 0] = np.copy(v1)
            c2_in[j, 0] = np.copy(v2)
            c1_in[j, 4] = np.copy(v1)
            c2_in[j, 4] = np.copy(v2)
            c1_in[j, 1] = np.copy(ev1)
            c2_in[j, 1] = np.copy(ev2)
            if (adump[1] == 'y'):    # prior on LDs
                c1_in[j, 5] = np.copy(lov1)
                c1_in[j, 6] = np.copy(hiv1)
                c2_in[j, 5] = np.copy(lov2)
                c2_in[j, 6] = np.copy(hiv2)

    dump = file.readline()
    dump = file.readline()
    # setup the arrays for the contamination factors
    cont = np.zeros((nfilt, 2))   # contains for each filter [value, error]

    # read the contamination factors
    for i in range(nfilt):
        dump = file.readline()
        adump = dump.split()
        # make sure the sequence in this array is the same as in the "filnames" array
        j = np.where(filnames == adump[0])
        cont[j, :] = [float(adump[1]), float(adump[2])]

    # read the stellar input properties
    dump = file.readline()
    dump = file.readline()
    dump = file.readline()
    adump = dump.split()
    Rs_in = float(adump[1])
    sRs_lo = float(adump[2])
    sRs_hi = float(adump[3])
    dump = file.readline()
    adump = dump.split()
    Ms_in = float(adump[1])
    sMs_lo = float(adump[2])
    sMs_hi = float(adump[3])
    dump = file.readline()
    adump = dump.split()
    howstellar = adump[1]
    dump = file.readline()
    adump = dump.split()
    sec_mass = adump[1]

    # read the MCMC setup
    dump = file.readline()
    dump = file.readline()
    adump = dump.split()
    nsamples = int(adump[1])   # total number of integrations
    dump = file.readline()
    adump = dump.split()
    nchains = int(adump[1])  # number of chains
    dump = file.readline()
    adump = dump.split()
    nproc = int(adump[1])  # number of processes
    dump = file.readline()
    adump = dump.split()
    burnin = int(adump[1])    # Length of bun-in
    dump = file.readline()
    adump = dump.split()
    walk = adump[1]            # Differential Evolution?
    dump = file.readline()
    adump = dump.split()
    grtest = True if adump[1] == 'y' else False  # GRtest done?
    dump = file.readline()
    adump = dump.split()
    plots = True if adump[1] == 'y' else False  # Make plots done
    dump = file.readline()
    adump = dump.split()
    leastsq = True if adump[1] == 'y' else False  # Do least-square?
    dump = file.readline()
    adump = dump.split()
    savefile = adump[1]   # Filename of save file
    dump = file.readline()
    adump = dump.split()
    savemodel = adump[1]   # Filename of model save file
    dump = file.readline()
    adump = dump.split()
    adaptBL = adump[1]   # Adapt baseline coefficent
    dump = file.readline()
    adump = dump.split()
    paraCNM = adump[1]   # remove parametric model for CNM computation
    dump = file.readline()
    adump = dump.split()
    # do a leas-square minimization for the baseline (not jump parameters)
    baseLSQ = adump[1]
    dump = file.readline()
    adump = dump.split()
    # use Levenberg-Marquardt algorithm for minimizer?
    lm = True if adump[1] == 'y' else False
    dump = file.readline()
    adump = dump.split()
    cf_apply = adump[1]  # which CF to apply
    dump = file.readline()
    adump = dump.split()
    jit_apply = adump[1]  # apply jitter

    inc_in = np.multiply(inc_in, np.pi)/180.
    # cos_in=[np.cos(inc_in[0]),np.multiply(np.sin(inc_in[0]),inc_in[1]),np.cos(inc_in[3]),np.cos(inc_in[2])]

    tarr = np.array([])  # initializing array with all timestamps
    farr = np.array([])  # initializing array with all flux values
    earr = np.array([])  # initializing array with all error values
    xarr = np.array([])  # initializing array with all x_shift values
    yarr = np.array([])  # initializing array with all y_shift values
    aarr = np.array([])  # initializing array with all airmass values
    warr = np.array([])  # initializing array with all fwhm values
    sarr = np.array([])  # initializing array with all sky values
    lind = np.array([])  # initializing array with the lightcurve indices
    barr = np.array([])  # initializing array with all bisector values
    carr = np.array([])  # initializing array with all contrast values

    indlist = []    # the list of the array indices
    bvars = []   # a list that will contain lists of [0, 1] for each of the baseline parameters, for each of the LCs. 0 means it's fixed. 1 means it's variable
    bvarsRV = []   # a list that will contain lists of [0, 1] for each of the baseline parameters, for each of the RV curves. 0 means it's fixed. 1 means it's variable

    # if ddFs are fit: set the Rp/Rs to the value specified at the jump parameters, and fix it.
    if ddfYN == 'y':
        rprs_in = [rprs_in[0], 0, 0, 1, 0, 0, 0]
        nddf = nfilt
    else:
        nddf = 0

    # set up the parameters
    params = np.array([T0_in[0], rprs_in[0], inc_in[0], aRs_in[0],
                       per_in[0], eos_in[0], eoc_in[0], K_in[0]])  # initial guess params
    stepsize = np.array([T0_in[1], rprs_in[1], inc_in[1], aRs_in[1],
                         per_in[1], eos_in[1], eoc_in[1], K_in[1]])  # stepsizes
    pmin = np.array([T0_in[2], rprs_in[2], inc_in[2], aRs_in[2],
                     per_in[2], eos_in[2], eoc_in[2], K_in[2]])  # Boundaries (min)
    pmax = np.array([T0_in[3], rprs_in[3], inc_in[3], aRs_in[3],
                     per_in[3], eos_in[3], eoc_in[3], K_in[3]])  # Boundaries (max)
    prior = np.array([T0_in[4], rprs_in[4], inc_in[4], aRs_in[4],
                      per_in[4], eos_in[4], eoc_in[4], K_in[4]])  # Prior centers
    priorlow = np.array([T0_in[5], rprs_in[5], inc_in[5], aRs_in[5],
                         per_in[5], eos_in[5], eoc_in[5], K_in[5]])  # Prior sigma low side
    priorup = np.array([T0_in[6], rprs_in[6], inc_in[6], aRs_in[6],
                        per_in[6], eos_in[6], eoc_in[6], K_in[6]])  # Prior sigma high side
    pnames = np.array(['T_0', 'RpRs', 'inc_[d]', 'aRs', 'Period_[d]',
                       'esin(w)', 'ecos(w)', 'K'])  # Parameter names

    if (divwhite == 'y'):           # do we do a divide-white? If yes, then fix all the transit shape parameters
        stepsize[0:6] = 0
        prior[0:6] = 0

    # if ddFs are fit: set the Rp/Rs to the specified value, and fix it.
    if ddfYN == 'y':
        drprs_in = np.zeros((nfilt, 7))
        njumpphot = njumpphot+1   # each LC has another jump pm

        # and make an array with the drprs inputs | drprs_op=[0.,float(adump[3]),float(adump[4]),float(adump[5])]  # the dRpRs options
        for i in range(nfilt):
            drprs_in[i, :] = drprs_op
            # add them to the parameter arrays
            params = np.concatenate((params, [drprs_in[i, 0]]))
            stepsize = np.concatenate((stepsize, [drprs_in[i, 1]]))
            pmin = np.concatenate((pmin, [drprs_in[i, 2]]))
            pmax = np.concatenate((pmax, [drprs_in[i, 3]]))
            prior = np.concatenate((prior, [drprs_in[i, 4]]))
            priorlow = np.concatenate((priorlow, [drprs_in[i, 5]]))
            priorup = np.concatenate((priorup, [drprs_in[i, 6]]))
            pnames = np.concatenate((pnames, [filnames[i]+'_dRpRs']))

    for i in range(nfilt):  # add the LD coefficients for the filters to the parameters
        params = np.concatenate(
            (params, [c1_in[i, 0], c2_in[i, 0], c3_in[i, 0], c4_in[i, 0]]))
        stepsize = np.concatenate(
            (stepsize, [c1_in[i, 1], c2_in[i, 1], c3_in[i, 1], c4_in[i, 1]]))
        pmin = np.concatenate(
            (pmin, [c1_in[i, 2], c2_in[i, 2], c3_in[i, 2], c4_in[i, 2]]))
        pmax = np.concatenate(
            (pmax, [c1_in[i, 3], c2_in[i, 3], c3_in[i, 3], c4_in[i, 3]]))
        prior = np.concatenate(
            (prior, [c1_in[i, 4], c2_in[i, 4], c3_in[i, 4], c4_in[i, 4]]))
        priorlow = np.concatenate(
            (priorlow, [c1_in[i, 5], c2_in[i, 5], c3_in[i, 5], c4_in[i, 5]]))
        priorup = np.concatenate(
            (priorup, [c1_in[i, 6], c2_in[i, 6], c3_in[i, 6], c4_in[i, 6]]))
        pnames = np.concatenate(
            (pnames, [filnames[i]+'_c1', filnames[i]+'_c2', filnames[i]+'_c3', filnames[i]+'_c4']))

    for i in range(nRV):
        params = np.concatenate((params, [gamma_in[i, 0]]), axis=0)
        stepsize = np.concatenate((stepsize, [gamma_in[i, 1]]), axis=0)
        pmin = np.concatenate((pmin, [gamma_in[i, 2]]), axis=0)
        pmax = np.concatenate((pmax, [gamma_in[i, 3]]), axis=0)
        prior = np.concatenate((prior, [gamma_in[i, 4]]), axis=0)
        priorlow = np.concatenate((priorlow, [gamma_in[i, 5]]), axis=0)
        priorup = np.concatenate((priorup, [gamma_in[i, 6]]), axis=0)
        pnames = np.concatenate((pnames, [RVnames[i]+'_gamma']), axis=0)

    # total number of baseline coefficients let to vary (leastsq OR jumping)
    nbc_tot = np.copy(0)

    for i in range(nphot):
        print(names[i])
        t, flux, err, xshift, yshift, airm, fwhm, sky, eti = np.loadtxt(
            fpath+names[i], usecols=(0, 1, 2, 3, 4, 5, 6, 7, 8), unpack=True)  # reading in the data
        if (divwhite == 'y'):  # if the divide - white is activated, divide the lcs by the white noise model before proceeding
            dwCNM = np.copy(dwCNMarr[dwCNMind[groups[i]-1]])
            flux = np.copy(flux/dwCNM)

        sky = sky-np.mean(sky)
        tarr = np.concatenate((tarr, t), axis=0)
        farr = np.concatenate((farr, flux), axis=0)
        earr = np.concatenate((earr, err), axis=0)
        xarr = np.concatenate((xarr, xshift), axis=0)
        yarr = np.concatenate((yarr, yshift), axis=0)
        aarr = np.concatenate((aarr, airm), axis=0)
        warr = np.concatenate((warr, fwhm), axis=0)
        sarr = np.concatenate((sarr, sky), axis=0)
        # bisector array: filled with 0s
        barr = np.concatenate((barr, np.zeros(len(t), dtype=np.int)), axis=0)
        # contrast array: filled with 0s
        carr = np.concatenate((carr, np.zeros(len(t), dtype=np.int)), axis=0)
        lind = np.concatenate((lind, np.zeros(len(t), dtype=np.int)+i), axis=0)
        indices = np.where(lind == i)
        indlist.append(indices)

        # the baseline coefficients for this lightcurve; each is a 2D array
        A_in, B_in, C_in, D_in, E_in, G_in, H_in, nbc = basecoeff(bases[i])

        nbc_tot = nbc_tot+nbc  # add up the number of jumping baseline coeff
        # if the least-square fitting for the baseline is turned on (baseLSQ = 'y'), then set the stepsize of the jump parameter to 0
        if (baseLSQ == "y"):
            abvar = np.concatenate(
                ([A_in[1, :], B_in[1, :], C_in[1, :], D_in[1, :], E_in[1, :], G_in[1, :], H_in[1, :]]))
            abind = np.where(abvar != 0.)
            bvars.append(abind)
            # the step sizes are set to 0 so that they are not interpreted as MCMC JUMP parameters
            A_in[1, :] = B_in[1, :] = C_in[1, :] = D_in[1,
                                                        :] = E_in[1, :] = G_in[1, :] = H_in[1, :] = 0

        # append these to the respective mcmc input arrays
        params = np.concatenate(
            (params, A_in[0, :], B_in[0, :], C_in[0, :], D_in[0, :], E_in[0, :], G_in[0, :], H_in[0, :]))
        stepsize = np.concatenate(
            (stepsize, A_in[1, :], B_in[1, :], C_in[1, :], D_in[1, :], E_in[1, :], G_in[1, :], H_in[1, :]))
        pmin = np.concatenate(
            (pmin, A_in[2, :], B_in[2, :], C_in[2, :], D_in[2, :], E_in[2, :], G_in[2, :], H_in[2, :]))
        pmax = np.concatenate(
            (pmax, A_in[3, :], B_in[3, :], C_in[3, :], D_in[3, :], E_in[3, :], G_in[3, :], H_in[3, :]))
        prior = np.concatenate((prior, np.zeros(len(A_in[0, :])+len(B_in[0, :])+len(
            C_in[0, :])+len(D_in[0, :])+len(E_in[0, :])+len(G_in[0, :])+len(H_in[0, :]))))
        priorlow = np.concatenate((priorlow, np.zeros(len(A_in[0, :])+len(B_in[0, :])+len(
            C_in[0, :])+len(D_in[0, :])+len(E_in[0, :])+len(G_in[0, :])+len(H_in[0, :]))))
        priorup = np.concatenate((priorup, np.zeros(len(A_in[0, :])+len(B_in[0, :])+len(
            C_in[0, :])+len(D_in[0, :])+len(E_in[0, :])+len(G_in[0, :])+len(H_in[0, :]))))
        pnames = np.concatenate((pnames, [names[i]+'_A0', names[i]+'_A1', names[i]+'_A2', names[i]+'_A3', names[i]+'_A4', names[i]+'_B1', names[i]+'_B2', names[i]+'_C1', names[i]+'_C2', names[i] +
                                          '_C3', names[i]+'_C4', names[i]+'_C5', names[i]+'_D1', names[i]+'_D2', names[i]+'_E1', names[i]+'_E2', names[i]+'_G1', names[i]+'_G2', names[i]+'_G3', names[i]+'_H1', names[i]+'_H2']))
        # note currently we have the following parameters in these arrays:
        #   [T0,RpRs,b,inc,per,eos, eoc, ddf_1, ..., ddf_n, c1_f1,c2_f1,c3_f1,c4_f1, c1_f2, .... , c4fn, A0_lc1,A1_lc1,A2_lc0,A3_lc0,A4_lc0, B0_lc1,B1_lc1,C0_lc1,C1_lc1,C2_lc1,C3_lc1,C4_lc1,D0_lc1,D1_lc1,E0_lc1,E1_lc1,G0_lc1,G1_lc1,H0_lc1,H1_lc1,H2_lc1,A0_lc2, ...]
        #    0  1    2  3   4  5   6    |  7 - 4+nddf     |          [7+nddf -- 6+nddf+4*n_filt]       |                           7+nddf+4*n_filt -- 7+nddf+4*n_filt + 15                                                |    7+nddf+4*n_filt + 16
        #    p a r a m e t e r s   |  d d F s        | Limb Darkening                             | const.  time  (5)     |      AM (2)  |     coordinate shifts   (5)      |     FWHM  (2) |   sky  (2)   | SIN (3)  | CNM (2) |
        #    each lightcurve has 21 baseline jump parameters, starting with index  8+nddf+4*n_filt+nRV

    for i in range(nRV):
        t, rv, err, bis, fwhm, contrast = np.loadtxt(
            fpath+RVnames[i], usecols=(0, 1, 2, 3, 4, 5), unpack=True)  # reading in the data

        tarr = np.concatenate((tarr, t), axis=0)
        # ! add the RVs to the "flux" array !
        farr = np.concatenate((farr, rv), axis=0)
        # ! add the RV errors to the "earr" array !
        earr = np.concatenate((earr, err), axis=0)
        # xshift array: filled with 0s
        xarr = np.concatenate((xarr, np.zeros(len(t), dtype=np.int)), axis=0)
        # yshift array: filled with 0s
        yarr = np.concatenate((yarr, np.zeros(len(t), dtype=np.int)), axis=0)
        # airmass array: filled with 0s
        aarr = np.concatenate((aarr, np.zeros(len(t), dtype=np.int)), axis=0)
        warr = np.concatenate((warr, fwhm), axis=0)
        # sky array: filled with 0s
        sarr = np.concatenate((sarr, np.zeros(len(t), dtype=np.int)), axis=0)
        barr = np.concatenate((barr, bis), axis=0)  # bisector array
        carr = np.concatenate((carr, contrast), axis=0)  # contrast array
        lind = np.concatenate(
            (lind, np.zeros(len(t), dtype=np.int)+i+nphot), axis=0)   # indices
        indices = np.where(lind == i+nphot)
        indlist.append(indices)

        # the baseline coefficients for this lightcurve; each is a 2D array
        W_in, V_in, U_in, S_in, nbcRV = basecoeffRV(RVbases[i])
        nbc_tot = nbc_tot+nbcRV  # add up the number of jumping baseline coeff
        abvar = np.concatenate(
            ([W_in[1, :], V_in[1, :], U_in[1, :], S_in[1, :]]))
        abind = np.where(abvar != 0.)
        njumpRV[i] = njumpRV[i]+len(abind)

        if (baseLSQ == "y"):
            bvarsRV.append(abind)
            # the step sizes are set to 0 so that they are not interpreted as MCMC JUMP parameters
            W_in[1, :] = V_in[1, :] = U_in[1, :] = S_in[1, :] = 0
        # append these to the respective mcmc input arrays
        params = np.concatenate(
            (params, W_in[0, :], V_in[0, :], U_in[0, :], S_in[0, :]))
        stepsize = np.concatenate(
            (stepsize, W_in[1, :], V_in[1, :], U_in[1, :], S_in[1, :]))
        pmin = np.concatenate(
            (pmin, W_in[2, :], V_in[2, :], U_in[2, :], S_in[2, :]))
        pmax = np.concatenate(
            (pmax, W_in[3, :], V_in[3, :], U_in[3, :], S_in[3, :]))
        prior = np.concatenate((prior, np.zeros(
            len(W_in[0, :])+len(V_in[0, :])+len(U_in[0, :])+len(S_in[0, :]))))
        priorlow = np.concatenate((priorlow, np.zeros(
            len(W_in[0, :])+len(V_in[0, :])+len(U_in[0, :])+len(S_in[0, :]))))
        priorup = np.concatenate((priorup, np.zeros(
            len(W_in[0, :])+len(V_in[0, :])+len(U_in[0, :])+len(S_in[0, :]))))
        pnames = np.concatenate((pnames, [RVnames[i]+'_W1', RVnames[i]+'_W2', RVnames[i]+'_V1',
                                          RVnames[i]+'_V2', RVnames[i]+'_U1', RVnames[i]+'_U2', RVnames[i]+'_S1', RVnames[i]+'_S2']))
        # note currently we have the following parameters in these arrays:
        #   [T0,RpRs,inc, aRs,per,eos, eoc, ddf_1, ..., ddf_n, c1_f1,c2_f1,c3_f1,c4_f1, c1_f2, .... , c4fn, A0_lc1,A1_lc1,A2_lc0,A3_lc0,A4_lc0, B0_lc1,B1_lc1,C0_lc1,C1_lc1,C2_lc1,C3_lc1,C4_lc1,D0_lc1,D1_lc1,E0_lc1,E1_lc1,G0_lc1,G1_lc1,H0_lc1,H1_lc1,H2_lc1,A0_lc2, ...]
        #    0  1    2  3   4  5   6    |  7 - 4+nddf     |          [7+nddf -- 6+nddf+4*n_filt]       |                           7+nddf+4*n_filt -- 7+nddf+4*n_filt + 15                                                |    7+nddf+4*n_filt + 16
        #    p a r a m e t e r s   |  d d F s        | Limb Darkening                             | const.  time  (5)     |      AM (2)  |     coordinate shifts   (5)      |     FWHM  (2) |   sky  (2)   | SIN (3)  | CNM (2) | time_rv (2) | bisector (2) | fwhm_rv (2) | contrast (2)
        #    each rv curve has 8 baseline jump parameters, starting with index  8+nddf+4*n_filt+nRV + nphot* 21

    # ============  CONTINUE HERE =================
    # calculate the weights for the lightcurves to be used for the CNM calculation later: do this in a function!
    ewarr = grweights(earr, indlist, grnames, groups, ngroup, nphot)

    for i in range(len(params)):
        print(pnames[i], params[i], stepsize[i],
              pmin[i], pmax[i], priorup[i], priorlow[i])

    inmcmc = 'n'
    indparams = [tarr, farr, xarr, yarr, warr, aarr, sarr, barr, carr, nphot, nRV, indlist, filters, nfilt, filnames, nddf, rprs0, erprs0, grprs, egrprs, grnames, groups,
                 ngroup, ewarr, inmcmc, paraCNM, baseLSQ, bvars, bvarsRV, cont, names, RVnames, earr, divwhite, dwCNMarr, dwCNMind, Rs_in, sRs_lo, sRs_hi, Ms_in, sMs_lo, sMs_hi]

    func = ff.fitfunc

    # ========= plot the initial lightcurves to see where we start=================
    yval0 = func(params, tarr, farr, xarr, yarr, warr, aarr, sarr, barr, carr, nphot, nRV, indlist, filters, nfilt, filnames, nddf, rprs0, erprs0, grprs, egrprs, grnames,
                 groups, ngroup, ewarr, inmcmc, paraCNM, baseLSQ, bvars, bvarsRV, cont, names, RVnames, earr, divwhite, dwCNMarr, dwCNMind, Rs_in, sRs_lo, sRs_hi, Ms_in, sMs_lo, sMs_hi)

    mcmc_plots(yval0, tarr, farr, earr, xarr, yarr, warr, aarr, sarr, barr,
               carr, lind, nphot, nRV, indlist, filters, names, RVnames, 'init_', params)

    # ================== calculate the CF factors! ===========================
    # first a chisq minimization

    inparams = np.copy(params)
    insteps = np.copy(stepsize)

    # for i in range(len(params)):
    #    print (params[i], stepsize[i])

    print('start least-square fit')

    # chisq, chibp, chim, lsfit = mc3.fit.modelfit(data=farr, uncert=earr, func=func, indparams=indparams,
    # params=params, pmin=pmin, pmax=pmax, stepsize=stepsize, lm=lm)

    output_opt = mc3.fit(farr, earr, func, params, indparams=indparams,
                         pmin=pmin, pmax=pmax, pstep=stepsize, leastsq='lm')

    chibp = output_opt['bestp']
    chim = output_opt['best_model']  # the model of the chi2 fit

    ######################GET BIC FROM LS###################    
    ijnames = np.where(stepsize != 0.)    
    jnames = pnames[[ijnames][0]]  # jnames are the names of the jump parameters
    nijnames = np.where(stepsize == 0.)
    njnames = pnames[[nijnames][0]]  # njnames are the names of the fixed parameters
    npar=len(jnames)

    ndat = len(tarr)
    if (baseLSQ == "y"):
          npar = npar + nbc_tot   # add the baseline coefficients if they are done by leastsq


    chisq = np.sum((chim-farr)**2/earr**2)    
    bic = chisq + npar * np.log(ndat)
    #bic = get_BIC(npar, ndat,chim)
    
    outF.write(str(time))
    outF.write(',')
    outF.write(str(am))
    outF.write(',')
    outF.write(str(xy))
    outF.write(',')
    outF.write(str(fw))
    outF.write(',')
    outF.write(str(sk))
    outF.write(',')
    outF.write(str(bic))
    outF.write(',')
    outF.write(str(npar))
    outF.write("\n")
outF.close()

#finds best model depending on the bayes factor

bic = []
ndat = []
new_bic = []
new_ndat = []
fin_bic = []
fin_ndat = []

time = []
am = []
xy = []
fw = []
sk = []

new_time = []
new_am = []
new_xy = []
new_fw = []
new_sk = []

fin_time = []
fin_am = []
fin_xy = []
fin_fw = []
fin_sk = []

best_time = []
best_am = []
best_xy = []
best_fw = []
best_sk = []

with open('myBICfile.csv') as csvFile:
    reader = csv.reader(csvFile)
    for row in reader:
        time.append(float(row[0]))    
        am.append(float(row[1]))    
        xy.append(float(row[2]))    
        fw.append(float(row[3]))    
        sk.append(float(row[4]))
        bic.append(float(row[5]))
        ndat.append(float(row[6]))
        min_bic = (min(bic))
    for i in range(len(bic)):
        if bic[i]-min_bic < 10:  #Bf=exp(deltaBIC/2)>150
            new_bic.append(bic[i])
            new_ndat.append(ndat[i])
            new_time.append(time[i])            
            new_am.append(am[i])            
            new_xy.append(xy[i])            
            new_fw.append(fw[i])            
            new_sk.append(sk[i])  
                      
print ('Lowest BIC',min_bic)                            
best_ndat = min(new_ndat)
for i in range(len(new_bic)): 
	if (new_ndat[i])==best_ndat:
	        fin_bic.append(new_bic[i])
	        fin_time.append(new_time[i])
	        fin_am.append(new_am[i])
	        fin_xy.append(new_xy[i])	
	        fin_fw.append(new_fw[i])	
	        fin_sk.append(new_sk[i])     
best_bic = min(fin_bic)

for i in range(len(new_bic)): 
	if (new_bic[i])==best_bic:
	        best_time.append(new_time[i])	
	        best_am.append(new_am[i])	
	        best_xy.append(new_xy[i])		        
	        best_fw.append(new_fw[i])		        
	        best_sk.append(new_sk[i])		        	                

print ('Best Model Parameters:')
print ('Number of free parameters', best_ndat)
print ('BIC', best_bic)
print ('Time', best_time)
print ('am', best_am)
print ('xy', best_xy)
print ('FWHM', best_fw)
print ('sky', best_sk)


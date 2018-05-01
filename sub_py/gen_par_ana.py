import numpy as np
import time
from math import sqrt as sqrt
from math import pi as pi
import os, sys
from subprocess import Popen, PIPE, STDOUT

cwd = os.getcwd()

from isotope import *
from ranges import *
from maxwellian import *

# g(enerate)p(arse)a(nalysis) for a provided isotope, energy, and name to call the output freya data. 
def gpa(Z, A, Energy, output_file, **kwargs):
    e = kwargs.get('e')
    x = kwargs.get('x')
    c = kwargs.get('c')
    T = kwargs.get('T')
    d = kwargs.get('d')
    generate_number = kwargs.get('generate_number')

    # if parameters are prescribed, rewrite the parameter file
    if e is not None:
        print("Begin Rewriting Parameter File...")
        if Energy == '-1' or -1:
            reac_type = 'sf'
        else:
            reac_type = '(n,f)'
        i = isotope(Z,A)[1]

        os.chdir(cwd+'/../fission_v2.0.3/data_freya/')
        infile = open("inputparameters.dat","r+") 
        
        content = infile.readlines() #reads line by line and outputs a list of each line
        content[i] = str(Z)+"  "+str(A)+"   '"+str(reac_type)+"'      "+str(e)+"     "+str(x)+"  "+str(c)+" "+str(T)+" 0.150  -                "+str(d)+"\n" 
        #replaces content of the (i+1)th line with chosen parameter values
        infile.close()
        infile = open("inputparameters.dat", 'w') #clears content of file. 
        infile.close
        infile = open("inputparameters.dat", 'r+')
        for item in content: #rewrites file content from list 
            infile.write("%s" % item)
        infile.close()
        os.chdir(cwd)
        print('Successful.')

    #  generate event files using the recipe executable for freya
    print('Begin Generating Events...')

    generate_begin = time.time()

    if generate_number is None:
        generate_number = 1000000

    generate_number = int(generate_number)

    if generate_number != 1000000:
        print('WARNING: 1,000,000 events is recommended')

    recipe_input = str(Z)+'\n'+str(A)+'\n'+str(Energy)+'\n'+ str(generate_number) + '\n'+str(output_file)+'\n'
    brecipe_input = recipe_input.encode()

    os.chdir(cwd )

    p = Popen("cd ../fission_v2.0.3/ \n ./recipe", shell = True, stdout=PIPE, stdin=PIPE, stderr=STDOUT)
    #  p = Popen("ls", shell = True, stdout=PIPE, stdin=PIPE, stderr=STDOUT)

    grep_stdout = p.communicate(input=b''+brecipe_input)[0]

    if grep_stdout is None:
        print('ERROR: Events not generated.')
    else:
        #  print(grep_stdout)
        print('Successful.')

    generate_end = time.time()
    gen_time = generate_end - generate_begin
    print('Time: ' + str(gen_time) + " sec")

    #  open output data made by recipe
    file = open(cwd+'/../fission_v2.0.3/sample_codes/events/build/' + output_file,'r')

    #  parse out first line
    initial_line = file.readline()
    initial_words = initial_line.split()

    freya_output = {}
    times = {}

    lines = file.readlines()
    lines = lines[0:len(lines)]
    file.close()

    #Build structure of observables
    #  as many rows as we have bins
    #  3 columns for bins, values, and variance
    #  no additional depth (in data_parse.py you will notice there is depth, this is to allow for the experimental uncertainty that is obviously not present here)

    Fragment_A = np.zeros((ranges_x['A'][1]- ranges_x['A'][0] + 1 , 3))
    Fragment_A[:,0] = range(ranges_x['A'][0],ranges_x['A'][1]+1)
    
    Product_A = np.copy(Fragment_A)

    n_Af = np.copy(Fragment_A)

    m_Af = np.copy(Fragment_A)

    TKE_A = np.copy(Fragment_A)

    #  egrid = np.logspace(-3,2,bin_number['mannhart'])
    egrid = np.zeros((bin_number['mannhart']+1))
    egrid[1:] = np.array(mannhart_bins)
    degrid = egrid[1:]-egrid[:-1]

    mannhart = np.zeros((bin_number['mannhart'], 3))
    mannhart[:,0]= mannhart_bins

    n_mult = np.zeros((ranges_x['n_mult'][1] - ranges_x['n_mult'][0] , 3))
    n_mult[:,0] = range(ranges_x['n_mult'][0],ranges_x['n_mult'][1])

    m_mult = np.zeros((ranges_x['m_mult'][1] - ranges_x['m_mult'][0] , 3))
    m_mult[:,0] = range(ranges_x['m_mult'][0],ranges_x['m_mult'][1] )

    energy_per_photon = np.zeros((ranges_x['A'][1]- ranges_x['A'][0]  , 4))
    energy_per_photon[:,0] = range(ranges_x['A'][0],ranges_x['A'][1] )

    total_photon_energy = np.copy(Fragment_A)

    #  the neutron and gamma spectrums, and the angular correlation will be calculated from the list which we will fill event by event. 
    n_spec = []
    rest_n_spec = []

    n_ang = []

    m_spec = []

    light_neutrons = []
    heavy_neutrons = []
    neutrons = []
    photons = []
    Freya_Ah = []
    Freya_Al = []
    freya_TKE = []

    #  the n_A_TKE bin columns need to be 'repeated' and 'tiled' so as to get every combinatorial pair
    n_A_TKE = np.zeros((bin_number['n_A_TKE'] * (ranges_x['A'][1] - ranges_x['A'][0]) , 4))
    n_A_TKE_energy_bins = np.arange(ranges_y['n_A_TKE'][0], ranges_y['n_A_TKE'][1])
    n_A_TKE_mass_bins = np.arange(ranges_x['A'][0],ranges_x['A'][1] )
    n_A_TKE[:,0] = np.repeat( n_A_TKE_mass_bins , int( len(n_A_TKE)  / len(n_A_TKE_mass_bins ) ) )
    n_A_TKE[:,1] = np.tile(n_A_TKE_energy_bins,int( len(n_A_TKE) / len(n_A_TKE_energy_bins) ) )

    n_TKE = np.zeros((bin_number['n_TKE'],4))
    n_TKE[:,0] = np.arange(ranges_x['n_TKE'][0], ranges_x['n_TKE'][1], bin_width['n_TKE'])


    import_begin = time.time()

    print('Importing Events...')

    n = 1
    events = []
    event = []
    for i in range(0,len(lines)):
        line = lines[i]
        if( i == (len(lines)-1) ):
            events.append(event[:])
        elif ( not(line.split()[0].strip() == str(n+1)) ):
            event.append(line.strip())
        else:
            events.append(event[:])
            event = [line.strip()]
            n += 1

    lines = [] #clear this portion of memory

    print(initial_words[3] + ' events imported.')
    import_end = time.time()
    import_time = import_end - import_begin
    print('Time: '+ str(import_time) + " sec")

    begin_parse = time.time()
    print('Parsing events...')

    def rest_frame_boost(old_energy , dot_product , fragment_energy , fragment_mass):
        new_energy = (np.sqrt(old_energy*2) + (dot_product * np.sqrt(fragment_energy * 2 * fragment_mass)))**2 / 2
        #  new_energy = old_energy + (dot_product * np.sqrt(fragment_energy) * fragment_mass)
        return new_energy

    for event in events:
        initial = []
        light = []
        heavy = []
        momenta_n = []
        momenta_m = []

        #  assign nn to be the number of neutrons for this event
        nn = int( event[0].split()[4] )
        if(nn > 0):
            nn_tf = 1
        else:
            nn_tf = 0
        initial = event[0:nn_tf+2]
        nucleus_direction = [float(initial[1].split()[1]),float(initial[1].split()[2]),float(initial[1].split()[3])]
        nucleus_TKE = float(initial[1].split()[0])
        nucleus_A = float(initial[0].split()[2]) - 1

        m = int( event[0].split()[5] )
        if (m>0):
            m_tf = 1
        else:
            m_tf = 0

        
        #  nnl is the integer number of neutrons coming from the light fragment
        nnl = int( event[nn_tf+2].split()[4] )

        #  ngl is the integer number of gammas coming from the light fragment
        ngl = int( event[nn_tf+2].split()[5] )
        #  if we have a positive number for nnl, set the nnl_tf to be 1 for true
        if(nnl > 0):
            nnl_tf = 1
        #  if we don't have any neutrons for the light fragment, set nnl_tf to 0 for false
        else:
            nnl_tf = 0
        if(ngl > 0):
            ngl_tf = 1
        else:
            ngl_tf = 0

        #  if we have (resp. don't have) neutrons for the fragments for an event we have (resp. don't have) additional lines to include in what we separate out as belonging to this event
        light = event[nn_tf+2:nn_tf+2+nnl_tf+ngl_tf+2]
        
        #  now we do the same thing for the heavy fragment
        nnh = int( event[nn_tf+2+nnl_tf+ngl_tf+2].split()[4] )
        ngh = int( event[nn_tf+2+nnl_tf+ngl_tf+2].split()[5] )

        if(nnh > 0):
            nnh_tf = 1
        else:
            nnh_tf = 0
        if(ngh > 0):
            ngh_tf = 1
        else:
            ngh_tf = 0
        heavy = event[nn_tf+2+nnl_tf+ngl_tf+2:nn_tf+2+nnl_tf+ngl_tf+2+nnh_tf+ngh_tf+2]

        Afrag_l = int( light[0].split()[2] )
        Afrag_h = int( heavy[0].split()[2] )

        #  write the appropriate values into the appropriate arrays
        #  the contents of each line for each events is detailed in the document: event_record.pdf


        Fragment_A[Afrag_l - ranges_x['A'][0] ][1] += 1
        Fragment_A[Afrag_h - ranges_x['A'][0] ][1] += 1
        
        n_Af[Afrag_l - ranges_x['A'][0] ][1] += nnl
        n_Af[Afrag_l - ranges_x['A'][0] ][2] += nnl**2

        n_Af[Afrag_h - ranges_x['A'][0] ][1] += nnh
        n_Af[Afrag_h - ranges_x['A'][0] ][2] += nnh**2

        m_Af[Afrag_l - ranges_x['A'][0] ][1] += ngl
        m_Af[Afrag_l - ranges_x['A'][0] ][2] += ngl**2

        m_Af[Afrag_h - ranges_x['A'][0] ][1] += ngh
        m_Af[Afrag_h - ranges_x['A'][0] ][2] += ngh**2

        Aprod_l = Afrag_l - nnl
        Aprod_h = Afrag_h - nnh

        light_neutrons.append(nnl)
        heavy_neutrons.append(nnh)
        Freya_Ah.append(Afrag_h)
        Freya_Al.append(Afrag_l)


        Product_A[Aprod_l - ranges_x['A'][0]][1] += 1
        Product_A[Aprod_h - ranges_x['A'][0]][1] += 1
        
        ntot = nnl + nnh + nn
        n_mult[ntot][1] += 1
        neutrons.append(ntot)
        
        mtot = ngl + ngh + m
        m_mult[mtot][1] += 1
        photons.append(mtot)

        TKE = float( light[1].split()[0] ) + float( heavy[1].split()[0] )

        light_TKE = float( light[1].split()[0] )
        heavy_TKE = float( heavy[1].split()[0] )
        light_direction = [float(light[1].split()[1]),float(light[1].split()[2]),float(light[1].split()[3])]
        heavy_direction = [float(heavy[1].split()[1]),float(heavy[1].split()[2]),float(heavy[1].split()[3])]

        freya_TKE.append(TKE)

        TKE_A[Afrag_l - ranges_x['A'][0]][1] += TKE
        TKE_A[Afrag_l - ranges_x['A'][0]][2] += TKE**2

        TKE_A[Afrag_h - ranges_x['A'][0]][1] += TKE
        TKE_A[Afrag_h - ranges_x['A'][0]][2] += TKE**2

        energy_bin =  np.searchsorted(n_A_TKE_energy_bins,TKE)
        energy_bin_l = (Afrag_l - ranges_x['A'][0])*len(n_A_TKE_energy_bins) + energy_bin
        n_A_TKE[energy_bin_l,2] += nnl
        n_A_TKE[energy_bin_l,3] += nnl**2
        energy_bin_h = (Afrag_h - ranges_x['A'][0]) * len(n_A_TKE_energy_bins) + energy_bin
        n_A_TKE[energy_bin_h,2] += nnh
        n_A_TKE[energy_bin_h,3] += nnh**2
        #  energy_bin_n = int((A - ranges_x['A'][0]) * len(n_A_TKE_energy_bins) + energy_bin)
        #  n_A_TKE[energy_bin_n,2] += nn
        #  n_A_TKE[energy_bin_n,3] += nn**2

        n_TKE_bin = np.searchsorted(n_TKE[:,0], TKE)
        n_TKE[n_TKE_bin,1] += ntot
        n_TKE[n_TKE_bin,2] += ntot ** 2
        n_TKE[n_TKE_bin,3] += 1

        #  treat each case of neutrons and gamma separately (as determined by the 1 or 0 value of the indicating boolean nnl_tf for ex.)
        #  for each of these cases we write the appropriate values into the appropriate lists and arrays

        if(nn_tf > 0):
            parts = initial[2].split()
            for j in range( 0,int(len(parts)/4) ):
                En = float(parts[j*4])
                momentum = []
                momentum.append( float(parts[j*4 +1]) )
                momentum.append( float(parts[j*4 +2]) )
                momentum.append( float(parts[j*4 +3]) )
                momenta_n.append( momentum )
                n_spec.append(En)
                dot_product = np.sum(np.array(momentum) * np.array(nucleus_direction))
                rest_n_spec.append(rest_frame_boost(En,dot_product,nucleus_TKE,nucleus_A))
        if( (nnl_tf > 0) and (ngl > 0) ):
            parts = light[2].split()
            for j in range( 0,int(len(parts)/4) ):
                En = float(parts[j*4])
                momentum = []
                momentum.append( float(parts[j*4 +1]) )
                momentum.append( float(parts[j*4 +2]) )
                momentum.append( float(parts[j*4 +3]) )
                momenta_n.append( momentum )
                n_spec.append(En)
                dot_product = np.sum(np.array(momentum) * np.array(light_direction))
                rest_En = rest_frame_boost(En,dot_product,light_TKE,Afrag_l)
                rest_n_spec.append(rest_En)
            parts = light[3].split()
            for j in range( 0,int(len(parts)/4) ):
                Eg = float(parts[j*4])
                momentum = []
                momentum.append( float(parts[j*4 +1]) )
                momentum.append( float(parts[j*4 +2]) )
                momentum.append( float(parts[j*4 +3]) )
                momenta_m.append( momentum )
                m_spec.append(Eg)

                energy_per_photon[Afrag_l - ranges_x['A'][0] + 1][1] += Eg
                energy_per_photon[Afrag_l - ranges_x['A'][0] + 1][2] += Eg**2
                energy_per_photon[Afrag_l - ranges_x['A'][0] + 1][3] += 1

                total_photon_energy[Afrag_l - ranges_x['A'][0] + 1][1] += Eg
                total_photon_energy[Afrag_l - ranges_x['A'][0] + 1][2] += Eg**2
        elif( nnl_tf > 0 ):
            parts = light[2].split()
            for j in range( 0,int(len(parts)/4) ):
                En = float(parts[j*4])
                momentum = []
                momentum.append( float(parts[j*4 +1]) )
                momentum.append( float(parts[j*4 +2]) )
                momentum.append( float(parts[j*4 +3]) )
                momenta_n.append( momentum )
                n_spec.append(En)
                dot_product = np.sum(np.array(momentum) * np.array(light_direction))
                rest_En = rest_frame_boost(En,dot_product , light_TKE,Afrag_l)
                rest_n_spec.append(rest_En)
        elif( ngl > 0 ):
            parts = light[2].split()
            for j in range( 0,int(len(parts)/4) ):
                Eg = float(parts[j*4])
                momentum = []
                momentum.append( float(parts[j*4 +1]) )
                momentum.append( float(parts[j*4 +2]) )
                momentum.append( float(parts[j*4 +3]) )
                momenta_m.append( momentum )
                m_spec.append(Eg)

                energy_per_photon[Afrag_l - ranges_x['A'][0] ][1] += Eg
                energy_per_photon[Afrag_l - ranges_x['A'][0] ][2] += Eg**2
                energy_per_photon[Afrag_l - ranges_x['A'][0] ][3] += 1

                total_photon_energy[Afrag_l - ranges_x['A'][0] ][1] += Eg
                total_photon_energy[Afrag_l - ranges_x['A'][0] ][2] += Eg**2

        if( (nnh_tf > 0) and (ngh > 0) ):
            parts = heavy[2].split()
            for j in range( 0,int(len(parts)/4) ):
                En = float(parts[j*4])
                momentum = []
                momentum.append( float(parts[j*4 +1]) )
                momentum.append( float(parts[j*4 +2]) )
                momentum.append( float(parts[j*4 +3]) )
                momenta_n.append( momentum )
                n_spec.append(En)
                dot_product = np.sum(np.array(momentum) * np.array(heavy_direction))
                rest_En = rest_frame_boost(En , dot_product , heavy_TKE,Afrag_h)
                rest_n_spec.append(rest_En)
            parts = heavy[3].split()
            for j in range( 0,int(len(parts)/4) ):
                Eg = float(parts[j*4])
                momentum = []
                momentum.append( float(parts[j*4 +1]) )
                momentum.append( float(parts[j*4 +2]) )
                momentum.append( float(parts[j*4 +3]) )
                momenta_m.append( momentum )
                m_spec.append(Eg)

                energy_per_photon[Afrag_h - ranges_x['A'][0] ][1] += Eg
                energy_per_photon[Afrag_h - ranges_x['A'][0] ][2] += Eg**2
                energy_per_photon[Afrag_h - ranges_x['A'][0] ][3] += 1

                total_photon_energy[Afrag_h - ranges_x['A'][0] ][1] += Eg
                total_photon_energy[Afrag_h - ranges_x['A'][0] ][2] += Eg**2
        elif( nnh_tf > 0 ):
            parts = heavy[2].split()
            for j in range( 0,int(len(parts)/4) ):
                En = float(parts[j*4])
                momentum = []
                momentum.append( float(parts[j*4 +1]) )
                momentum.append( float(parts[j*4 +2]) )
                momentum.append( float(parts[j*4 +3]) )
                momenta_n.append( momentum )
                n_spec.append(En)
                dot_product = np.sum(np.array(momentum) * np.array(heavy_direction))
                rest_En = rest_frame_boost(En , dot_product , heavy_TKE,Afrag_h)
                rest_n_spec.append(rest_En)
        elif( ngh > 0 ):
            parts = heavy[2].split()
            for j in range( 0,int(len(parts)/4) ):
                Eg = float(parts[j*4])
                momentum = []
                momentum.append( float(parts[j*4 +1]) )
                momentum.append( float(parts[j*4 +2]) )
                momentum.append( float(parts[j*4 +3]) )
                momenta_m.append( momentum )
                m_spec.append(Eg)

                energy_per_photon[Afrag_h - ranges_x['A'][0] ][1] += Eg
                energy_per_photon[Afrag_h - ranges_x['A'][0] ][2] += Eg**2
                energy_per_photon[Afrag_h - ranges_x['A'][0] ][3] += 1

                total_photon_energy[Afrag_h - ranges_x['A'][0] ][1] += Eg
                total_photon_energy[Afrag_h - ranges_x['A'][0] ][2] += Eg**2
        
        for p in range(0,len(momenta_n)-1):
            for q in range(p+1,len(momenta_n)):
                dotprod = momenta_n[p][0] * momenta_n[q][0] + momenta_n[p][1] * momenta_n[q][1] + momenta_n[p][2] * momenta_n[q][2]
                n_ang.append(dotprod)

    print('All events parsed.')
    end_parse = time.time()
    parse_time = end_parse - begin_parse
    print('Time: '+ str(parse_time) + ' sec')

    begin_calc = time.time()
    print('Calculating Statistics... ')

    #  calculate variances, averages for all the appropriate observables
    #  (using numpy functions for efficiency)
    with np.errstate(divide='ignore', invalid='ignore'):

        long_frag_counts = np.repeat(Fragment_A[:,1][:-1],int(len(n_A_TKE)/(ranges_x['A'][1] - ranges_x['A'][0])))
        n_A_TKE[:,2] = np.divide(n_A_TKE[:,2], long_frag_counts)
        n_A_TKE[:,2]= 100 * n_A_TKE[:,2]
        #  n_A_TKE[:,3] = np.divide(n_A_TKE[:,3], long_frag_counts)
        n_A_TKE[:,3] = np.subtract(n_A_TKE[:,3],np.square(n_A_TKE[:,2]))
        n_A_TKE[:,3] = np.sqrt(n_A_TKE[:,3])
        n_A_TKE[:,2] = np.array( n_A_TKE[:,2], dtype = np.float)
        n_A_TKE[n_A_TKE == 0] = 'nan'

        total_photon_energy[:,1] = np.divide(total_photon_energy[:,1], Fragment_A[:,1])
        total_photon_energy[:,2] = np.divide(total_photon_energy[:,2], Fragment_A[:,1])
        total_photon_energy_squared = np.divide(total_photon_energy[:,1], Fragment_A[:,1])
        total_photon_energy[:,2] = np.subtract(total_photon_energy[:,2], total_photon_energy_squared)
        total_photon_energy[:,2] = np.sqrt(total_photon_energy[:,2])
        
        energy_per_photon[:,1] = np.divide(energy_per_photon[:,1], energy_per_photon[:,3])
        energy_per_photon[:,2] = np.divide(energy_per_photon[:,2], energy_per_photon[:,3])
        energy_per_photon[:,2] = np.subtract(energy_per_photon[:,2], np.square(energy_per_photon[:,1]))
        energy_per_photon[:,2] = np.sqrt(energy_per_photon[:,2])

        average_photon_energy = np.average(m_spec)
        ape_sigma = np.sqrt(np.mean(np.square(m_spec)) - average_photon_energy**2)

        n_Af[:,1] = np.divide(n_Af[:,1] , Fragment_A[:,1])
        n_Af[:,2] = np.divide(n_Af[:,2] , Fragment_A[:,1])
        n_Af[:,2] = np.subtract(n_Af[:,2],np.square(n_Af[:,1]))
        n_Af[:,2] = np.sqrt(n_Af[:,2])

        TKE_A[:,1] = np.divide(TKE_A[:,1] , Fragment_A[:,1])
        TKE_A[:,2] = np.divide(TKE_A[:,2] , Fragment_A[:,1])
        TKE_A[:,2] = np.subtract(TKE_A[:,2],np.square(TKE_A[:,1]))
        TKE_A[:,2] = np.sqrt(TKE_A[:,2])

        m_Af[:,1] = np.divide(m_Af[:,1],Fragment_A[:,1])
        m_Af[:,2] = np.divide(m_Af[:,2],Fragment_A[:,1])
        m_Af[:,2] = np.subtract(m_Af[:,2],np.square(m_Af[:,1]))
        m_Af[:,2] = np.sqrt(m_Af[:,2])


        hlab_freya , binEdges_freya = np.histogram(n_spec,bins = egrid)
        #  binCenters_freya = 0.5 * (binEdges_freya[1:] + binEdges_freya[:-1])
        n = len(n_spec)

        mannhart[:,0] = binEdges_freya[:-1]
        mannhart[:,1] = hlab_freya/degrid/n
###
        #  mannhart[:,2] = 10/(np.sqrt(hlab_freya))
        mannhart[:,2] = 2/(np.sqrt(hlab_freya))
###

        n_TKE[:,1] = np.divide(n_TKE[:,1], n_TKE[:,3])
        n_TKE[:,2] = np.divide(n_TKE[:,2], n_TKE[:,3])
        n_TKE[:,2] = np.subtract(n_TKE[:,2], np.square(n_TKE[:,1]))
        n_TKE[:,2] = np.sqrt(n_TKE[:,2])

        nubar = float(len(n_spec))/float(len(events))
        nubar_sigma = np.sqrt(np.mean(np.square(np.array(neutrons))) - nubar**2)
        gammabar = float(len(m_spec)) / float(len(events))
        gammabar_sigma = np.sqrt(np.mean(np.square(np.array(photons))) - gammabar**2)

        n_mult[:,2] = np.sqrt(1/n_mult[:,1])
        n_mult[:,1] = np.divide(n_mult[:,1] , np.sum(n_mult[:,1]) )

###
        m_mult[:,2] = np.sqrt(1/m_mult[:,1])/2
###
        m_mult_smudge = np.copy(m_mult)
        #  m_mult_smudge[1:-1,1] = m_mult[1:-1,1] + m_mult[:-2,1] + m_mult[2:,1]
        m_mult_smudge[1:-1,1] = m_mult[:-2,1] + m_mult[2:,1]

        m_mult_smudge[:,1] = np.divide(m_mult_smudge[:,1] , np.sum(m_mult_smudge[:,1]) )
        m_mult[:,1] = np.divide(m_mult[:,1] , np.sum(m_mult[:,1]) )

        nu1 = 0.0
        nu1_sigma = 0.0
        nu2 = 0.0
        nu2_sigma = 0.0
        nu3 = 0.0
        nu3_sigma = 0.0
        nu4 = 0.0
        nu4_sigma = 0.0

        for i in n_mult[:,0]:
            i = int(i)
            if sum(n_mult[:,1]) != 0:
                nu1 += (i) * n_mult[i,1]
                nu1_sigma += (i)**2 * n_mult[i,1]
                nu2 += (i) * (i-1) * n_mult[i,1]
                nu2_sigma += (i)**2 * (i-1)**2 * n_mult[i,1]
                nu3 += (i) * (i-1) * (i-2) *  n_mult[i,1]
                nu3_sigma += (i)**2 * (i-1)**2 * (i-2)**2 *  n_mult[i,1]
                nu4 += (i) * (i-1) * (i-2) * (i-3) * n_mult[i,1]
                nu4_sigma += (i)**2 * (i-1)**2 * (i-2)**2 * (i-3)**2 * n_mult[i,1]
        nu1_sigma = np.sqrt(nu1_sigma - nu1**2)
        nu2_sigma = np.sqrt(nu2_sigma - nu2**2)
        nu3_sigma = np.sqrt(nu3_sigma - nu3**2)
        nu4_sigma = np.sqrt(nu4_sigma - nu4**2)

        end_calc = time.time()

    #  now for the histogram calculations (angular, spectrums)
    #  take the lists written to above, and calculate the appropriate statistics

    number_of_bins = int(np.ceil(bin_number['n_angular']))
    hist, bin_edges = np.histogram(n_ang, bins = number_of_bins, range = (ranges_x['n_angular'][0],ranges_x['n_angular'][1]), normed = True)
    n_angular = np.zeros((number_of_bins,3))
    n_angular[0] = [ -1 , 1 , bin_width['n_angular']]
    n_angular[:,0] = np.add(bin_edges, bin_width['n_angular']/2)[:-1]
    n_angular[:,1] = hist
    n_angular[:,2] = np.nan

    n_spectrum_bin_widths = bin_width['n_spectrum']
    number_of_bins = int(np.ceil(bin_number['n_spectrum']))
    hist, bin_edges = np.histogram(n_spec, bins = number_of_bins, range = (ranges_x['n_spectrum'][0],ranges_x['n_spectrum'][1]))
    hist = hist.astype('float')
    hist[hist == 0] = None
    normed_hist, bin_edges = np.histogram(n_spec, bins = number_of_bins, range = (ranges_x['n_spectrum'][0],ranges_x['n_spectrum'][1]), normed = True)
    n_spectrum = np.zeros((number_of_bins ,3))
    n_spectrum[0] = [0,ranges_x['n_spectrum'][1],n_spectrum_bin_widths]
    n_spectrum[:,0] = bin_edges[:-1]
    n_spectrum[:,1] = normed_hist
    #  n_spectrum[:,2] = n_spectrum_bin_widths * 1/np.sqrt(hist)
    n_spectrum[:,2] = np.nan

    n_spectrum_bin_widths = bin_width['n_spectrum']
    number_of_bins = int(np.ceil(bin_number['n_spectrum']))
    hist, bin_edges = np.histogram(rest_n_spec, bins = number_of_bins, range = (ranges_x['n_spectrum'][0],ranges_x['n_spectrum'][1]))
    hist = hist.astype('float')
    hist[hist == 0] = None
    normed_hist, bin_edges = np.histogram(rest_n_spec, bins = number_of_bins, range = (ranges_x['n_spectrum'][0],ranges_x['n_spectrum'][1]), normed = True)
    rest_n_spectrum = np.zeros((number_of_bins ,3))
    rest_n_spectrum[0] = [0,ranges_x['n_spectrum'][1],n_spectrum_bin_widths]
    rest_n_spectrum[:,0] = bin_edges[:-1]
    rest_n_spectrum[:,1] = normed_hist
    #  rest_n_spectrum[:,2] = rest_n_spectrum_bin_widths * 1/np.sqrt(hist)
    rest_n_spectrum[:,2] = np.nan

    m_spectrum_bin_widths = bin_width['m_spectrum']
    number_of_bins = int(np.ceil(bin_number['m_spectrum']))
    hist, bin_edges = np.histogram(m_spec, bins = number_of_bins, range = (ranges_x['m_spectrum'][0] - m_spectrum_bin_widths, ranges_x['m_spectrum'][1]), normed = True)
    m_spectrum = np.zeros((number_of_bins,3))
    m_spectrum[0] = [0,ranges_x['m_spectrum'][1],m_spectrum_bin_widths]
    m_spectrum[:,0] = np.add(bin_edges,m_spectrum_bin_widths/2)[:-1]
    m_spectrum[:,1] = hist
    m_spectrum[:,2] = np.nan

    Fragment_A[:,1] = np.divide(Fragment_A[:,1], 0.5 * np.average(Fragment_A[1:,1]))

    Product_A[:,1] = np.divide(Product_A[:,1], 0.5 * np.average(Product_A[1:,1]))

    print('Time: ' + str(end_calc - begin_calc) + ' sec')
    if initial_words[1] is 253 or '253' or "253":
        initial_words[1] = 252

    # assign elements of dictionary to appropriate freya outputs
    freya_output['Z00'] = initial_words[0]
    freya_output['A00'] = initial_words[1]
    freya_output['isotope'] = 'Z = ' + str(freya_output['Z00']) + ' A = ' + str(freya_output['A00'])
    freya_output['number_of_events'] = initial_words[3]
    freya_output['nubar_moments'] = [np.array([[1,nu1,nu1_sigma],[2,nu2,nu2_sigma],[3,nu3,nu3_sigma]]),
                    "Moments",
                    "Neutron multiplicity moments for: "+freya_output['isotope'],
                    "Moment Number",
                    "Value"]
    freya_output['nubar'] = [np.array([[[1],[nubar],[nubar_sigma]]]),"nubar","nubar","N/A","nubar"]
    freya_output['nu1'] = np.array([nu1,nu1_sigma])
    freya_output['nu2'] = np.array([nu2,nu2_sigma])
    freya_output['nu3'] = np.array([nu3,nu3_sigma])
    freya_output['nu4'] = np.array([nu4,nu4_sigma])
    freya_output['gammabar'] = [np.array([[[1],[gammabar],[gammabar_sigma]]]),"gammabar","gammabar","N/A","gammabar"]
    freya_output['Fragment_A'] = Fragment_A
    freya_output['Product_A'] = [Product_A,
                    'Product Mass (A)',
                    'Mass of Fission Products for: ' + freya_output['isotope'],
                    'Fragment Mass (A)',
                    'Probability']
    freya_output['n_A_TKE'] = [n_A_TKE,
                    'nu(A,TKE)', 
                    'Neutron Multiplicity vs Mass Number A and Total Kinetic Energy for: ' + freya_output['isotope'],
                    'Mass Number(A)',
                    'Total Kinetic Energy',
                    'Neutron Multiplicity'
                    ]
    freya_output['mannhart'] = [mannhart,
                    r'$ \bar\nu (E) $',
                    'Nubar vs. Total Kinetic Energy for: '+freya_output['isotope'],
                    'Total Kinetic Energy of Fission Fragments (MeV)', 'Mean Neutrons Multiplicity (Count/Bin Width)']
    freya_output['n_angular'] = [n_angular,
                    'Total Angular Correlation',
                    r'$n-n \ Angular\ Correlation\ for: \ $'+freya_output['isotope'],
                    r'$ \cos ( \theta )$',
                    'Correlation']
    freya_output['n_mult'] = [n_mult,
                    'Neutron Multiplicity Distribution',
                    r'$Neutron \ Multiplicity\ Distribution\ for: \ $'+freya_output['isotope'],
                    'Neutron Multiplicity', 
                    'Probability'] 
    freya_output['m_mult'] = [m_mult,
                    'Gamma Multiplicity',
                    r'$Gamma \ Multiplicity\ Distribution\ for: \ $'+freya_output['isotope'],
                    'Gamma Multiplicity',
                    'Probability']
    freya_output['m_mult_smudge'] = [m_mult_smudge,
                    'Gamma Multiplicity',
                    r'$Gamma \ Multiplicity\ Distribution\ for: \ $'+freya_output['isotope'],
                    'Gamma Multiplicity',
                    'Probability']
    freya_output['n_spectrum'] = [n_spectrum,
                    'Neutron Spectral Shape',
                    r'$Fission \ Neutron \ Spectrum\ for: \ $'+freya_output['isotope'], 
                    'Neutron Energy (MeV)', 
                    'Probability'] 
    freya_output['rest_n_spectrum'] = [rest_n_spectrum,
                    'Neutron Spectral Shape',
                    r'$Fission \ Neutron \ Spectrum\ for: \ $'+freya_output['isotope'], 
                    'Neutron Energy (MeV)', 
                    'Probability'] 
    freya_output['m_spectrum'] = [m_spectrum,
                    'Photon Spectral Shape',
                    r'$Fission \ Photon \ Spectrum\ for: \ $'+freya_output['isotope'], 
                    'Photon Energy(MeV)', 
                    'Probability'] 
    freya_output['total_photon_energy'] = [total_photon_energy,
                    'Total Photon Energy',
                    r'$Total \ Photon \ Energy\ for: \ $'+freya_output['isotope'],
                    'Fragment Mass A', 
                    'Total Photon Energy (MeV)'] 
    freya_output['energy_per_photon'] = [energy_per_photon,
                    'Energy Per Photon',
                    r'$Energy \ per \ Photon\ for: \ $'+freya_output['isotope'], 
                    'Fragment Mass A', 'Energy per Photon (MeV)']
    freya_output['average_photon_energy'] = [np.array([[[1],[average_photon_energy],[ape_sigma]]]),"average_photon_energy","average_photon_energy","N/A","average_photon_energy"]
    freya_output['n_Af'] = [n_Af,
                    r'$ \bar \nu(A)$',
                    r'$Neutron \ Multiplicity \ vs. \ Fragment \ Mass \ for: \ $'+freya_output['isotope'], 
                    'Fragment Mass A', 
                    r'$\bar \nu(A)$']
    freya_output['TKE_A'] = [TKE_A,
                    r'$  TKE(A)$',
                    r'$Total \ Kinetic \ Energy vs. \ Fragment \ Mass \ for: \ $'+freya_output['isotope'], 
                    'Fragment Mass A', 
                    r'$ TKE (A)$']
    freya_output['m_Af'] = [m_Af,
                    'Gamma Multiplicity',
                    r'$Gamma \ Multiplicity \ vs. \ Fragment \ Mass \ for: \ $'+freya_output['isotope'], 
                    'Fragment Mass A', 
                    'Gamma Multiplicity'] 
    freya_output['n_TKE'] = [n_TKE,
                    r'$ \bar \nu(TKE)$',
                    r'$Neutron \ Multiplicity \ vs. \ Total \ Kinetic \ Energy \ for: \ $'+freya_output['isotope'], 
                    'TKE (MeV)', 
                    r'$\bar \nu(TKE)$' ]

    freya_output['light_neutrons'] = np.array(light_neutrons)
    freya_output['heavy_neutrons'] = np.array(heavy_neutrons)
    freya_output['Ah'] = np.array(Freya_Ah)
    freya_output['Al'] = np.array(Freya_Al)
    freya_output['TKE'] = np.array(freya_TKE)
    freya_output['neutrons'] = np.array(neutrons)
    freya_output['n_TKE_alt'] = freya_output['n_TKE']


    times['gen_time'] = gen_time
    times['parse_time'] = parse_time + import_time

    return freya_output, times


import pandas as pd
import numpy
import rebound

G=6.67428e-11 # m3 kg-1 s-2
M_sun=1.98855e30 # kg
AU=1.496e11 # m

#-------------------------------------------------------------------------------
def load_dat_file(f):
    '''Load dat.txt file and return a time and a dataframe containing all particle data

    Parameters
    ----------
    f
        path to data file, formatted with a time stamp on the first line, and then x,y,z,vx,vy,vz,m,r; one line per particle
    '''
    arrays = [numpy.array(map(float,line.split())) for line in open(f)] #reads in line by line
    t=float(arrays[0]) #first line is time
    dat=numpy.array(arrays[1:])

    '''N=len(dat[:,0])
    pos=dat[:,:3]
    vel=dat[:,3:6]
    m=dat[:,6]
    r=dat[:,7]
    return pos'''

    df_dat=pd.DataFrame(dat,columns=['x(m)','y(m)','z(m)','vx(ms^-1)','vy(ms^-1)','vz(ms^-1)','m(kg)','r(m)'])
    return t,df_dat
#-------------------------------------------------------------------------------
def file_list_no_restart(files):
    ''' take a list of files dat0000000_0.txt, dat0000001_0.txt, dat0000001_1.txt, ...
    and account for an restarts
    input should be sorted'''

    files.sort()
    final_files=[files[0]]
    for j in range(1,len(files)):
        fnum=int(files[j][3:10]) # number of current file
        fnum2=int(files[j-1][3:10]) # number of previous file
        if fnum!=fnum2: # if they are not equal, we can safely append
            final_files.append(files[j])
        else:
            final_files[-1]=files[j] # otherwise we must reset the value of the last entry
    return final_files
#-------------------------------------------------------------------------------
def load_run_params(file_path):
    '''
    Load the run_params file as a pandas dataframe.
    Generally expects a run that has an identifying number somewhere in the name.
    ----------
    file_path
        path to run_params txt file
    '''

    params_columns=['run','run_dir','N_tot','a_orb(m)','R_eq(m)','rho(kgm-3)',
    'M_tot(kg)','R_c(m)','OM_orb(s-1)','X','OM_circ(s-1)','f','dt(s)','t_run(s)',
    'n_cores','initial']
    df_params= pd.DataFrame([],columns=params_columns)
    # params=numpy.genfromtxt(file_path,dtype=None) # this lines causes warnings
    with open(file_path) as f:
        content = f.readlines()
    params = [x.strip() for x in content]

    d=None
    run_name_split=params[0].split('/')[-1].split('_')
    for i in range(len(run_name_split))[::-1]:
        try:
            d=int(run_name_split[i]) #load a run with run number at front of dirname
        except:
            continue

    params=numpy.insert(params,0,0) #insert placeholder for run number
    df_params_add=pd.DataFrame([params],columns=params_columns)
    df_params_add['run']=d
    df_params=df_params.append(df_params_add,ignore_index=True)
    print df_params
    return df_params
#-------------------------------------------------------------------------------
def rotating_to_heliocentric_array(r,v,a,t):
    '''Function to transform from rotating frame (circular orbit of radius a, at R(t=0)=(a,0,0)) to the heliocentric frame'''
    Om_k=numpy.sqrt(G*(M_sun)/(a**3.0)) #angular velocity at semimajor axis a_k
    a_vec=a*numpy.array([numpy.cos(Om_k*t),numpy.sin(Om_k*t),0.0]) #time gives location of frame
    Om_vec=numpy.array([0.0,0.0,Om_k]) # define angular velocity vector
    theta=Om_k*t # angle of rotationa s a function of time
    rot_mat=numpy.array([numpy.cos(theta),-numpy.sin(theta),0,numpy.sin(theta),numpy.cos(theta),0,0,0,1]).reshape((3,3)) # z axis rotation vector
    Om_skew=numpy.array([0,-Om_k,0,Om_k,0,0,0,0,0]).reshape((3,3)) #define skew symmetric matrix for angular velocity (0,0,Om_k)
    rot_mat_dot=numpy.dot(Om_skew,rot_mat) #time derivative of the rotation matrix

    if len(r.shape)>1:
        N=len(r)
        # R=numpy.zeros((N,3))
        # V=numpy.zeros((N,3))
        for i in range(N):
            # R[i,:]=a_vec+numpy.dot(rot_mat,r[i,:]) # transform position
            # V[i,:]=numpy.cross(Om_vec,a_vec)+numpy.dot(rot_mat_dot,r[i,:])+numpy.dot(rot_mat,v[i,:]) # transform velocity
            R=a_vec+numpy.dot(rot_mat,r.T).T # transform position
            V=numpy.cross(Om_vec,a_vec)+numpy.dot(rot_mat_dot,r.T).T+numpy.dot(rot_mat,v.T).T # transform velocity
    else:
        R=a_vec+numpy.dot(rot_mat,r) # transform position
        V=numpy.cross(Om_vec,a_vec)+numpy.dot(rot_mat_dot,r)+numpy.dot(rot_mat,v) # transform velocity
    return R,V
#-------------------------------------------------------------------------------
def load_orb_file(f):
    '''Function to load orbits file as a dataframe'''
    arrays = [numpy.array(map(float,line.split())) for line in open(f)] #reads in line by line
    orb_dat=numpy.array(arrays)

    try:
        df_orb=pd.DataFrame(orb_dat,columns=['t(s)','i_o','j','file','a(m)','e','I(rad)',
        'omega(rad)','OMEGA(rad)','f_true(rad)','E_J(J)','m1(kg)','m2(kg)','a_hel(m)',
        'e_hel','I_hel(rad)','omega_hel(rad)','OMEGA_hel(rad)','f_true_hel(rad)'])
    except:
        df_orb=pd.DataFrame(orb_dat,columns=['t(s)','file','i','j','a(m)','e','I(rad)',
        'omega(rad)','OMEGA(rad)','f_true(rad)','m1(kg)','m2(kg)'])

    return df_orb
#-------------------------------------------------------------------------------
def orbit_func_faster(f,df_params,m_lim_1=0,m_lim_2=2,N_lim=100,coord_trans=1):
    '''
    Search a file for orbits, returns a list of orbits (if any)

    f - file name including path
    df_params - pandas dataframe containing the run parameters
    m_lim_1 - minimum mass ratio of interest. NB that this is the mass ratio to the largest body, could smaller body binaries be hidden?
    m_lim_2 - minimum factor that particle mass must increase from the minimum
    N_lim - max number of particles to search for orbits
    coord_trans=1 #set to 1 if we want to transform from rotating to helio

    '''

    print "search {}".format(f)

    #create rebound sim once at start
    sim = rebound.Simulation()
    sim.G=G

    orb_list=[] # list will remain empty if no orbits are found

    fnum=int(f.split("/")[-1][3:10])

    # retrieve from run params
    m_min=float(df_params.iloc[0]['M_tot(kg)'])/float(df_params.iloc[0]['N_tot'])

    # load the data file
    t,df_dat=load_dat_file(f)

    # Sort from highest to lowest mass and keep only the first N_lim particles
    df_dat=df_dat.sort_values(by=['m(kg)'],ascending=False)
    m_cut_1=numpy.amax(df_dat['m(kg)'])*m_lim_1 # Only consider pairs of particles above a certain mass ratio
    m_cut_2=m_min*m_lim_2 # only consider particles that have accreted to a mass m_cut_2
    df_dat=df_dat[(df_dat['m(kg)']>=m_cut_1) & (df_dat['m(kg)']>m_cut_2)]
    N=len(df_dat)
    if N<=1: # deal with case of no viable particles, and therefore no orbits
        return []
    if N>N_lim: # Have an additional cut, never search between > N_lim particles
        df_dat=df_dat.iloc[:N_lim]
        N=len(df_dat)
    df_dat['i']=range(len(df_dat)) # add an index to rank by mass, 0 is most massive
    # print f,N

    #Do coordinate transform from the simulation rotating frame, to the heliocentric frame
    if coord_trans==1:
        print "coord transform"
        r=numpy.array(df_dat.loc[:,['x(m)','y(m)','z(m)']])
        v=numpy.array(df_dat.loc[:,['vx(ms^-1)','vy(ms^-1)','vz(ms^-1)']])
        R,V=rotating_to_heliocentric_array(r,v,float(df_params.iloc[0]['a_orb(m)']),t)
        # # OPTIONAL: transform to relative frame
        # R=R-R[0,:]
        # V=V-V[0,:]
        df_dat.loc[:,['x(m)','y(m)','z(m)']]=R
        df_dat.loc[:,['vx(ms^-1)','vy(ms^-1)','vz(ms^-1)']]=V
    else:
        print "no coord transform"
    # Add particles to rebound sim
    for i in df_dat['i']:
        pi=df_dat.iloc[i] #primary properties
        sim.add(x=pi['x(m)'],y=pi['y(m)'],z=pi['z(m)'],
        vx=pi['vx(ms^-1)'],vy=pi['vy(ms^-1)'],vz=pi['vz(ms^-1)'],
        m=pi['m(kg)'])

    # Search for orbits
    for i in df_dat['i'][:-1]: # the index means we don't search the last particle, it will have already been searched (symmetry!)
        #print "search particle {} out of {}".format(i,N-1)
        pi=df_dat.iloc[i] #primary properties
        for j in df_dat['i']:
            if i==j or j<i: # Do not search for orbit with self, or with a particle of higher mass (it's already been checked)
                continue
            pj=df_dat.iloc[j] # secondary properties

            # use rebound to calculate the orbit, better handling for certain cases, e.g. circular orbits
            orbit = sim.particles[j].calculate_orbit(sim.particles[i])
            # print "number of rebound particles = {}".format(len(sim.particles))
            # Only save closed orbits
            if orbit.e>=0 and orbit.e<1.0:
                # Save the orbit as: ['t(s)','file','i','j','a(m)','e','I(rad)','omega(rad)','OMEGA(rad)','f_true(rad)','m1(kg)','m2(kg)'] where i,j is now the actual particle index in the file
                # add to array
                orb_list.append([t,fnum,int(pi.name),int(pj.name),orbit.a,orbit.e,orbit.inc,orbit.omega,orbit.Omega,orbit.f,pi['m(kg)'],pj['m(kg)']])

    return orb_list # returns a list of orbits
#-------------------------------------------------------------------------------
def planet_orbit(orb,n):
    '''function to find the xyz orbit data relative to the reference point'''
    #orb is the 6 orbital elements: a,e,I,OM,om,f and n is number of points to plot with
    a=orb[0]
    e=orb[1]
    I=orb[2]
    OMEGA=orb[3]
    omega=orb[4]
    #specify theta from 0 to 2pi radians: i.e. number of points on orbit, THE TRUE ANOMALY
    theta = numpy.linspace(0,2*numpy.pi, n)
    #-------------------------------------------------------------------------------
    r = numpy.zeros(len(theta))
    for i in range(len(theta)):
        r[i]=r_elliptical(a,e,theta[i])
    #-------------------------------------------------------------------------------
    #determine orbital basis vectors
    e_p=find_e_p(OMEGA,omega,I)
    e_Q=find_e_Q(OMEGA,omega,I)
    #-------------------------------------------------------------------------------
    #define the r(x,y,z) position array
    pos=numpy.zeros((len(theta),3))
    for n in range(len(pos[:,0])):
        pos[n,:]=r[n]*((numpy.cos(theta[n])*e_p)+(numpy.sin(theta[n])*e_Q))#PSS eq 11.36a
        pos[n,:]=pos[n,:]
    return pos
#-------------------------------------------------------------------------------
def r_elliptical(a,e,theta):
    '''
    Function to find distance, r, at true anomaly, theta, around a keplerian orbit

    Parameters
    ----------
    a
        semi major axis (m)
    e
        orbital eccentricity
    theta
        true anomaly (rad)
    '''
    r = (a*(1-e**2))/(1+e*numpy.cos(theta))
    return r
#-------------------------------------------------------------------------------
def find_e_p(OMEGA,omega,I):
    '''Function to find the normalised Lenz vector along the apsidal line, uses PSS eq11.39a'''
    e_p=numpy.zeros(3)
    e_p[0]=(numpy.cos(OMEGA)*numpy.cos(omega))-(numpy.cos(I)*numpy.sin(OMEGA)*numpy.sin(omega)) #x
    e_p[1]=(numpy.sin(OMEGA)*numpy.cos(omega))+(numpy.cos(I)*numpy.cos(OMEGA)*numpy.sin(omega)) #y
    e_p[2]=numpy.sin(I)*numpy.sin(omega) #z
    return e_p
#-------------------------------------------------------------------------------
def find_e_Q(OMEGA,omega,I):
    '''Function to find the normalised e_Q vector which is perp to h and e_a and
    lies in the orbital plane, uses PSS eq11.39b'''
    e_Q=numpy.zeros(3)
    e_Q[0]=-(numpy.cos(OMEGA)*numpy.sin(omega))-(numpy.cos(I)*numpy.sin(OMEGA)*numpy.cos(omega)) #x
    e_Q[1]=-(numpy.sin(OMEGA)*numpy.sin(omega))+(numpy.cos(I)*numpy.cos(OMEGA)*numpy.cos(omega)) #y
    e_Q[2]=numpy.sin(I)*numpy.cos(omega) #z
    return e_Q
#-------------------------------------------------------------------------------

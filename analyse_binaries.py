import binorb_func as bf
import pandas as pd
import matplotlib.pyplot as pyplot
import matplotlib.gridspec as gridspec
import numpy

# change these parameters to load a particular run
run_path="test_data" # location of the directory with the run data
run_dir="test_run_0" # name of the directory with the run data
orb_file="test_run_0_orbit_search_faster_hel.txt" # name of the orbit file
run_params_file="{}/{}/run_params_0.txt".format(run_path,run_dir)

# These parameters control some plotting options
orbit_plot_limit=5 # maximum number of orbits to show on plot
plot_lim_factor=2.0 # The plot limits are set to be (-lim, +lim) where lim = plot_lim_factor * semimajor axis of the first orbit

# load the run params
df_rp=bf.load_run_params(run_params_file)
print df_rp

# load the orbit file
df_orb=bf.load_orb_file(orb_file)

# save the orbit dataframe to a file if you want to keep it
# df_orb.to_csv("df_{}".format(orb_file.split("/")[-1]))
# df_orb=pd.read_csv("df_{}".format(orb_file.split("/")[-1]),index_col=0)
# print df_orb

# retrieve the particle indices that have bound orbits
pri_list=numpy.array(df_orb['i']).astype(int)
sec_list=numpy.array(df_orb['j']).astype(int)
particle_list=numpy.unique(numpy.append(pri_list,sec_list))
# print particle_list

# load the data file and keep only the particles that are in particle_list
n_file=int(numpy.unique(df_orb['file'])[0])
file_iter=0 # some dat file might be of the format dat_..._0.txt, dat_..._1.txt, dat_..._2.txt etc
while file_iter<10:
    dat_file_name="{}/{}/dat{:07d}_{}.txt".format(run_path,run_dir,n_file,file_iter)
    try:
        t,df_dat=bf.load_dat_file(dat_file_name)
        break
    except:
        file_iter+=1
df_cut=df_dat.iloc[particle_list]

# N.B. if the dat file contains >1e5 particles it will be slow to load and cut.
# You may want to save the cut dataframe once and then load it instead
# df_cut.to_csv("cut_dat{:07d}_0.txt".format(n_file))
# t=bf.read_dat_file_time(dat_file_name)
# df_cut=pd.read_csv("cut_dat{:07d}_0.txt".format(n_file),index_col=0)

df_dat=df_cut

# transform from rotating to helio coords
pos_hel,vel_hel=bf.rotating_to_heliocentric_array(numpy.array(df_dat[['x(m)','y(m)','z(m)']]),numpy.array(df_dat[['vx(ms^-1)','vy(ms^-1)','vz(ms^-1)']]),float(df_rp['a_orb(m)']),t)
df_dat[['x(m)','y(m)','z(m)']]=pos_hel
df_dat[['vx(ms^-1)','vy(ms^-1)','vz(ms^-1)']]=vel_hel

# transform to coords relative to the position of the primary
pri=int(df_orb.iloc[0]['i'])
# print pri,df_dat.loc[pri][['x(m)','y(m)','z(m)']]
df_dat[['x(m)','y(m)','z(m)']]=df_dat[['x(m)','y(m)','z(m)']]-df_dat.loc[pri][['x(m)','y(m)','z(m)']]
# print df_dat

fig = pyplot.figure()
gs = gridspec.GridSpec(1,1)
ax1 = pyplot.subplot(gs[0,0])
ax1.set_aspect("equal")

# plot all orbits relative to the primary of that orbit
for k in range(len(df_orb)):
    if k>orbit_plot_limit: # only plot some of the orbits
        break
    df=df_orb.iloc[k]
    pri=int(df['i'])
    pos_pri=df_dat.loc[pri][['x(m)','y(m)','z(m)']]
    orb=numpy.array(df[['a(m)','e','I(rad)','OMEGA(rad)','omega(rad)']])
    pos_orb=bf.planet_orbit(orb,100)
    ax1.plot((pos_orb[:,0]+pos_pri[0]),(pos_orb[:,1]+pos_pri[1]),zorder=2)

# plot all particles (in particle_list), where colour/size scales with mass
df_dat=df_dat.sort_values('m(kg)')
# print df_dat
s1=ax1.scatter(df_dat['x(m)'],df_dat['y(m)'],c=numpy.log10(df_dat['m(kg)']),zorder=3)

lim=plot_lim_factor*float(df_orb.iloc[0]['a(m)'])
ax1.set_xlim(-lim,lim)
ax1.set_ylim(-lim,lim)
cbar=pyplot.colorbar(s1)
cbar.set_label("log10(m(kg))")
ax1.set_xlabel("x(m)")
ax1.set_ylabel("y(m)")

pyplot.show()

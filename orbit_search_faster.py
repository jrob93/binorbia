'''
python script that launches the multiprocess orbit search.

Accepts as input the path and dir to search.
-p path to where the run directories live
-d name of the run directory to be searched

e.g.:
python orbit_search_faster.py -p test_data -d test_run_0

Extra options:
choose data files to search with:
-i starting file number
-j end file number
-n file step

Note that if you pass -ve numbers for i and j, i and j will be used as list indices.
e.g. "-i -2 -j -1" chooses the last data file
'''
# coding: utf-8

from multiprocessing import Pool
import time
import numpy
import binorb_func as bf
import pandas as pd
import rebound
import subprocess
import os
from optparse import OptionParser

#input options
parser = OptionParser()
parser.add_option( "-p", "--path", dest="path", help="path", metavar="path")
parser.add_option( "-d", "--dir", dest="dir", help="dir", metavar="dir" )
parser.add_option( "-t", "--threads", dest="threads", help="threads", metavar="threads" )
parser.add_option( "-s", "--save_path", dest="save_path", help="save_path", metavar="save_path" )
parser.add_option( "-m", "--m_lim_1", dest="m_lim_1", help="m_lim_1", metavar="m_lim_1" )
parser.add_option( "-M", "--m_lim_2", dest="m_lim_2", help="m_lim_2", metavar="m_lim_2" )
parser.add_option( "-N", "--N_lim", dest="N_lim", help="N_lim", metavar="N_lim" )
parser.add_option( "-i", "--file_low", dest="file_low", help="file_low", metavar="file_low" )
parser.add_option( "-j", "--file_high", dest="file_high", help="file_high", metavar="file_high" )
parser.add_option( "-n", "--file_step", dest="file_step", help="file_step", metavar="file_step" )
parser.add_option( "-o", "--orb_f_base", dest="orb_f_base", help="orb_f_base", metavar="orb_f_base" )
parser.add_option( "-T", "--coord_trans", dest="coord_trans", help="coord_trans", metavar="coord_trans" )

(options,args)=parser.parse_args()

if options.path:
    path=options.path
else:
    cmap="."
if options.dir:
    dir=options.dir
else:
    dirname="."
if options.threads:
    threads=int(options.threads)
else:
    threads=1
if options.save_path:
    save_path=options.save_path
else:
    save_path="."
if options.m_lim_1:
    m_lim_1=float(options.m_lim_1)
else:
    m_lim_1=0
if options.m_lim_2:
    m_lim_2=float(options.m_lim_2)
else:
    m_lim_2=2
if options.N_lim:
    N_lim=int(options.N_lim)
else:
    N_lim=100
if options.file_step:
    file_step=int(options.file_step)
else:
    file_step=1
if options.orb_f_base:
    orb_f_base=options.orb_f_base
else:
    orb_f_base="orbit_search_faster_hel"
if options.coord_trans:
    coord_trans=int(options.coord_trans)
else:
    coord_trans=1

print "path:",path
print "dir:",dir
print "threads:",threads
print "save_path:",save_path
print "m_lim_1:",m_lim_1
print "m_lim_2:",m_lim_2
print "N_lim:",N_lim
print "coord_trans:",coord_trans

# path="restart_dirs"
# dir="226_cloud_order_kelvin_fix_dt_2_cloud_runs_fix_dt_2_56381"

dir_path="{}/{}".format(path,dir)
orb_fname="{}/{}_{}.txt".format(save_path,dir,orb_f_base)
print "dir path: {}".format(dir_path)
print "save file: {}".format(orb_fname)

# find the files in the dir to search
files=next(os.walk("{}".format(dir_path)))[2] #retrieve the files in the run directory
files.sort() #ensure that the files are always sorted the same?
d_files = [ fi for fi in files if fi.endswith(".txt") and fi.startswith("dat") ] #keep only dat*.txt files
d_files = bf.file_list_no_restart(d_files) #account for restarted runs

# set up the dat file list
if options.file_low:
    if int(options.file_low)<0:
        file_low=int(d_files[int(options.file_low)][3:10])
    else:
        file_low=int(options.file_low)
else:
    file_low=int(d_files[0][3:10])
if options.file_high:
    if int(options.file_high)<0:
        file_high=int(d_files[int(options.file_high)][3:10])
    else:
        file_high=int(options.file_high)
else:
    file_high=int(d_files[-1][3:10])
print "file low:",file_low
print "file high:",file_high
print "file_step:",file_step
d_files = [ fi for fi in d_files if int(fi[3:10])>=file_low and int(fi[3:10])<=file_high ] #keep only dat*.txt files in the desired ranges
d_end=d_files[-1] # always keep the last file
d_files=d_files[::file_step] # make the required steps
if d_files[-1]!=d_end: # add the last file if neccessary
    d_files.append(d_end)
d_files=d_files[::-1] # reverse the order of the files (search from end to start)

# load the run params file
rp_files = [ fi for fi in files if fi.endswith(".txt") and fi.startswith("run_params") ]
print rp_files[0]
df_params=bf.load_run_params("{}/{}".format(dir_path,rp_files[0]))

# set up pool and find orbits
t_start=time.time()
pool = Pool(threads)
multiple_results = [pool.apply_async(bf.orbit_func_faster, args=("{}/{}".format(dir_path,f),df_params,m_lim_1,m_lim_2,N_lim,coord_trans)) for f in d_files]
pool.close()
pool.join()

# recombine all the orbits arrays
orbits_all=[]
count=0
for i,f in enumerate(d_files):
    # print i,f
    orbits=multiple_results[i].get()
    if len(orbits)==0:
         print f,numpy.shape(orbits),numpy.shape(orbits_all)
         continue
    else:
        if count==0:
            orbits_all=orbits
            count+=1
        else:
            orbits_all=numpy.vstack((orbits_all,orbits))
        print f,numpy.shape(orbits),numpy.shape(orbits_all)

# save the arrays as a dataframe
df_orb=pd.DataFrame(orbits_all)
df_orb.to_csv(orb_fname,sep="\t",header=False,index=False)
print orb_fname
t_end=time.time()

print "time = {}s".format(t_end-t_start)

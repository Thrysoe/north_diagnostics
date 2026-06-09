"""
Only the imput parameter section should be modified.

They are: - [shot] int
              the number of the shot
          - [bias_type] str
              the way to bias the probes in the given shot
              Note that it can only be a temperature or ion_saturation_current bias
              Other bias type will trigger an error message
          - [data_type] str
              the type of data plotted
              Can be raw data or fluctuations (data minus mean value)
          - [data_origin] str
              the way data were collected. 
              Can be experiment or simulation.
          - [current_value] float
              the value of the current in the vertical coils, only relevant when [data_origin] is simulation
          - [path] str
              the name of the folder in which all the code architecture is 
          - [fps] float
              the frame rate of the saved movie 
          - [length] int
              the time lenght in µs of the movie/image sequence
          - [sus_probes] numpy array
              the discarded probes considering their current or voltage plots
          - [k] int
              the extent of the colorbar in terms of standard variation of the signal

Take the data from the DDAQ file or computed by the extract_all_data program to plot images of NORTH and then 
generate a video. The obtained plots are stored in a sub folder the Figure folder created if not existing.
Can be raw or fluctuating data.
"""

#Find the path to local libraries (modify the PYTHON PATH)
import sys
sys.path.append('C:\\Users\\tulla\\Perso\\north_diagnostics\\diagnostics')
sys.path.append('C:\\Users\\tulla\\Perso\\north_diagnostics\\utils')
sys.path.append('C:\\Users\\tulla\\Perso\\north_diagnostics')

#Import all useful libraries
import numpy as np
import matplotlib.pyplot as plt
import os
import utils.dau



"""
Input parameters
"""
shot = 10918
bias_type = 'sweep'
data_type = 'raw'
data_origin = 'experiment'
current_value = 0
path = 'C:/Users/tulla/Perso/north_diagnostics'
fps = 15
length = 100000
sus_probes = []
k = 2




"""
Main program: from now, all shall remain untouched.
"""
#Extract and cleaning data from .txt file
#Activated probes for that type of data
studied_probes = np.array([1,2,3,5,6,7,8,13,14,15,20,21,23,24,25, 
                           26,28,29,30,31,32,33,38,39,41,42,44,45,46,47,48,49])

#data = np.genfromtxt(f"{path}/Data/probe_data{shot}.txt", delimiter=';', skip_header=1)
if data_origin == 'experiment':
      start = 300000
      all_data = utils.dau.read_probe_data(shot, path, bias_type, 1/100, 
                                         studied_probes-1, 0.7, 0.8, False)
      time = np.mean(all_data[:,0,0])
      data = np.mean(all_data, axis=0)
      data[0,:] = time
      print(f"Data shot {shot} loaded with success at time {time} s")
    
      #Create a folder to store data if it doesn't exist
      if not os.path.exists(f"{path}/Figures/{shot}_{bias_type}_{data_type}/"):
          os.makedirs(f"{path}/Figures/{shot}_{bias_type}_{data_type}/")
          print("A new directory for storing data was created")
        
elif data_origin == 'simulation': 
      #Activated probes for that type of data
      studied_probes = np.array([1,2,3,4,5,6,7,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,
                                 26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50])
      
      all_data = np.genfromtxt(f"{path}/Data/probeIsat_Bvertsweep.txt", skip_header=5)
      print("Data probe loaded with success")
      
      #Take only the data for the given current value
      n_sweep = int(current_value//2.5)
      data = all_data[n_sweep*1001+100:(n_sweep+1)*1001-100, 1:]
      data[:, 0] = data[:, 0]*1E-3
      
      #Create a folder to store data if it doesn't exist
      if not os.path.exists(f"{path}/Figures/{int(current_value)}_{bias_type}_{data_type}/"):
        os.makedirs(f"{path}/Figures/{int(current_value)}_{bias_type}_{data_type}/")
        print("A new directory for storing data was created") 
        
else:
    print('WARNING: the data origin is not recognized')
    sys.exit()
    
#Boundary condition and name of each map
bc_all = [-2E-5, 17000, -10, 0, 0, 0, -10, 2E16]
sweep_all = ['Ion_saturation_current', 'Temperature', 'Floating Potential',
             'Ion_saturation_current relative uncertainty', 'Temperature relative uncertainty', 
             'Floating Potential relative uncertainty', 'Plasma potential', 'Plasma density']

print(f"Data ({data_type}) loaded with success")

#Testing activated probes
activated_probes = []
for i in range(50):
    if i+1 in studied_probes and i+1 not in sus_probes:
        activated_probes.append(True)
    else:
        activated_probes.append(False)
        


for j in range(len(bc_all)):
    #Create figure 
    fig = plt.figure(figsize=(10,10))
    bc = bc_all[j]
    
    #Data plotting
    if j==1:
        vmin = 0
    if j==len(bc_all)-2:
        Vp = data[:, 2]+1.38E-23*data[:, 1]/(2*1.6E-19)*np.log(np.sqrt(6.68E-27/(2*np.pi*9.11E-31*0.61**2)))
        Vp[0] = time
        mean, std = np.mean(Vp[studied_probes]), np.std(Vp[studied_probes]) #[:, studied_probes]
        vmin, vmax = mean-k*std, mean+k*std
        output = utils.dau.plot_2D_data(Vp,shot, path, bias_type, data_type, j, activated_probes, vmin, vmax, 
                                        fig, bc)
    elif j==len(bc_all)-1:
        dens = -data[:, 0]/(0.61*1E-6*1.6E-19*np.sqrt(1.38E-23*(data[:, 1]+1)/6.68E-27))
        dens[0] = time
        mean, std = np.mean(dens[studied_probes]), np.std(dens[studied_probes]) #[:, studied_probes]
        vmin, vmax = mean-k*std, mean+k*std 
        output = utils.dau.plot_2D_data(dens, shot, path, bias_type, data_type, j, activated_probes, vmin, vmax, 
                                        fig, bc)
    else:       
        mean, std = np.mean(data[studied_probes, j]), np.std(data[studied_probes, j]) #[:, studied_probes]
        vmin, vmax = mean-k*std, mean+k*std
        output = utils.dau.plot_2D_data(data[:, j], shot, path, bias_type, data_type, j, 
                                    activated_probes, vmin, vmax, fig, bc)
    plt.legend()
    std_format = "{:0"+str(int(np.log(len(data)+1)/np.log(10))+1)+"d}"
    std_format = std_format.format(i)
    if data_origin=='experiment':
        plt.suptitle(f"Shot {shot} {data_type} {bias_type} "+sweep_all[j])
        plt.savefig(f"{path}/Figures/{shot}_{bias_type}_{data_type}/{sweep_all[j]}.jpg")
    else:
        plt.suptitle(f"Simulation {data_type} {bias_type} at time {data[i,0]:.6} s")
        plt.savefig(f"{path}/Figures/{int(current_value)}_{bias_type}_{data_type}/{std_format}.jpg")
   
    print(output)

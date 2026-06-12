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
          - [toroidal_current_value] float
              the value of the current in the toroidal coils, only relevant when [data_origin] is simulation
          - [vertical_current_value] float
              the value of the current in the vertical coils, only relevant when [data_origin] is simulation
          - [path] str
              the name of the folder in which all the code architecture is 
          - [fps] float
              the frame rate of the saved movie 
          - [start] int
              the time in µs at which the movie will begin
          - [length] int
              the time lenght in µs of the movie/image sequence
          - [sus_probes] numpy array
              the discarded probes considering their current or voltage plots
          - [k] int
              the extent of the colorbar in terms of standard variation of the signal

Take the data from the DDAQ file or computed by the simulation to plot images of a poloidal cross-section
of NORTH and then generate a video. The CRIO file is needed for the heating zone.
The obtained plots are stored in a sub folder the Figure folder created if not existing.
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
shot = 11152
bias_type = 'ion_saturation_current'
data_type = 'fluctuations' #'fluctuations'
data_origin = 'simulation'
toroidal_current_value = 1000
vertical_current_value = 30
path = 'C:/Users/tulla/Perso/north_diagnostics'
fps = 15
start = 400000
length = 1000
sus_probes = [22]
k = 2




"""
Main program: from now, all shall remain untouched.
"""
#Extract and cleaning data from .txt file
#Activated probes for that type of data
studied_probes = np.array([1,2,3,5,6,7,8,9,10,11,13,14,15,17,18,19,20,21,23,24,25,
                           26,28,29,30,31,32,33,38,39,41,42,43,44,45,46,47,48,49])

#data = np.genfromtxt(f"{path}/Data/probe_data{shot}.txt", delimiter=';', skip_header=1)
if data_origin == 'experiment':
      all_data = utils.dau.read_probe_data(shot, path, bias_type, 1/75, 
                                         studied_probes-1, 0, 1, False)
      data = all_data[start:start+length, :]
      print(f"Data shot {shot} loaded with success")
    
      #Create a folder to store data if it doesn't exist
      if not os.path.exists(f"{path}/Figures/{shot}_{bias_type}_{data_type}/"):
          os.makedirs(f"{path}/Figures/{shot}_{bias_type}_{data_type}/")
          print("A new directory for storing data was created")
        
elif data_origin == 'simulation': 
      #Activated probes for that type of data
      studied_probes = np.array([1,2,3,4,5,6,7,9,10,11,12,13,14,15,16,17,18,20,21,22,23,24,25,
                                 26,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50])
      
      all_data = np.genfromtxt(f"{path}/Data/probeIsat_{int(vertical_current_value)}_{toroidal_current_value}.txt", skip_header=5)
      print("Data probe loaded with success")
      
      #Take only the data for the given current value
      start = 3000
      n_sweep = 0 #int((current_value)//20)
      data = - all_data[n_sweep*5001+start:n_sweep*5001+start+length, :]*1E3
      data[:, 0] = - data[:, 0]*1E-6
      
      #Create a folder to store data if it doesn't exist
      if not os.path.exists(f"{path}/Figures/{int(vertical_current_value)}_{toroidal_current_value}_{bias_type}_{data_type}/"):
        os.makedirs(f"{path}/Figures/{int(vertical_current_value)}_{toroidal_current_value}_{bias_type}_{data_type}/")
        print("A new directory for storing data was created") 
        
else:
    print('WARNING: the data origin is not recognized')
    sys.exit()
    
#Add a boundary condition for the interpolation method
bc = - 2.26E-5 #in A

#If activated, compute the fluctuations of the data and adapt boundary condition
if data_type=='fluctuations':
  data[:,1:] = data[:,1:] - np.mean(data[:,1:], axis=0, keepdims=True)
  bc = 0
print(f"Data ({data_type}) loaded with success")

#Testing activated probes
activated_probes = []
for i in range(50):
    if i+1 in studied_probes and i+1 not in sus_probes:
        activated_probes.append(True)
    else:
        activated_probes.append(False)

#Parameters for plotting a single colorbar
mean, std = np.mean(data[:, studied_probes]), np.std(data[:, studied_probes])
vmin, vmax = mean - k*std, mean + k*std

#Create figure 
fig = plt.figure(figsize=(10,10))
  
#Plot all images of the vessel + the probes + the data for each time
for i in range(len(data[:,0])): 
    data[i, 1:] = data[i, 1:]*1E3
    output = utils.dau.plot_2D_data(data[i, :], shot, path, bias_type, data_type, data_origin, 
                             toroidal_current_value, i, activated_probes, vmin*1E3, vmax*1E3, fig, bc*1E3)
    plt.legend(fontsize=18)
    
    #Necessary to name images with a fixed format to avoid confusion when the video is read
    std_format = "{:0"+str(int(np.log(len(data[:,0])+1)/np.log(10))+1)+"d}"
    std_format = std_format.format(i)
    
    #Adapt the layout
    if data_origin=='experiment':
        plt.suptitle(f"Shot {shot} {data_type} {bias_type} at time {data[i,0]:.6} s", fontsize=18)
        plt.savefig(f"{path}/Figures/{shot}_{bias_type}_{data_type}/{std_format}.jpg")
    else:
        plt.suptitle(f"Simulation {data_type} {bias_type} at time {data[i,0]:.6} s", fontsize=18)
        plt.savefig(f"{path}/Figures/{int(vertical_current_value)}_{toroidal_current_value}_{bias_type}_{data_type}/{std_format}.jpg")
    plt.clf()
    
    #Display a confirmation message that the process happened
    print(output)
plt.close()

#Save those images in a movie (.avi format) => to be found in the figure folder
output = utils.dau.video_2D(data, data_origin, shot, toroidal_current_value, vertical_current_value, path, 
                            bias_type, data_type, fps, int(start/1E3))
print(output)

"""
Only the imput parameter section should be modified.

They are: - [shot] int
              the number of the shot
          - [bias_type] str
              the way to bias the probes in the given shot
              Note that it can only be a temperature or ion_saturation_current bias
              Other bias type will trigger an error message
          - [T_sweep] float
              the period of the sweep of the bias voltage. 
              Only relevant for temperature mesurements
          - [studied_channels] numpy array
              the number of the probe collecting data in the given shot
              Array with all probes: np.array([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,
                                       26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50])
          - [path] str
              the name of the folder in which all the code architecture is 
          - [t_start] float
              the start time of the data collection. 
              Must be between 0 and 1 second
          - [t_end] float
              the end time of the data collection. 
              Must be between 0 and 1 second and greater than t_start


Extract all data from the DDAQ and CRIO files and plot it. Both .tdms file must be saved in the data folder.
Save all useful data in a .txt file with a short header description, plot machine parameters and I-t and
U-t curves to see if all of the data collected aren't anomalous.
"""

#Find the path to local libraries (modify the PYTHON PATH)
import sys
sys.path.append('C:\\Users\\tulla\\Perso\\north_diagnostics\\diagnostics')
sys.path.append('C:\\Users\\tulla\\Perso\\north_diagnostics\\utils')
sys.path.append('C:\\Users\\tulla\\Perso\\north_diagnostics')

#Import all useful libraries
import matplotlib.pyplot as plt
import numpy as np
import os
import utils.dau



"""
Input parameters
"""
shot = 10919
bias_type = 'ion_saturation_current'
T_sweep = 50 # in s, one other the sweeping voltage frequency
studied_probes = np.array([1,2,3,5,6,7,8,9,10,13,14,15,17,18,19,20,21,23,24,25,
                           26,28,29,30,31,32,33,38,39,41,42,44,45,46,47,48,49])
path = 'C:/Users/tulla/Perso/north_diagnostics'
t_start = 300E-3 # in s
t_end = 310E-3 # in s





"""
Main program: from now, all shall remain untouched.
"""
#Create a folder to store images if they don't exist
if not os.path.exists(f"{path}/Figures/IandVplots_{shot}/"):
  os.makedirs(f"{path}/Figures/IandVplots_{shot}/")
  print("A new directory for storing data was created")

#Generate the data and plot the curves
machine_data = utils.dau.read_machine_data(shot, path)
probe_data = utils.dau.read_probe_data(shot, path, bias_type, T_sweep, studied_probes-1, 
                                       t_start, t_end, True)
plt.close()

#Saving all data in the Data folder
#Machine data
head = 'Time; Light sensor; Coil current; Pressure sensor; LFS power; HFS power in SI units'
np.savetxt(f"{path}/Data/machine_data{shot}.txt", machine_data, delimiter=';', header=head)

#Probes data
#head = 'Time; probes in the numerical order in SI units'
#np.savetxt(f"{path}/Data/probe_data{shot}.txt", probe_data, delimiter=';', header=head)

#Plot and save machine parameters
fig = plt.figure(figsize=(10,10))

#Light sensor
plt.subplot(2,2,1)
plt.plot(machine_data[:,0], machine_data[:,1])
plt.xlabel('time (s)')
plt.ylabel('light sensor signal (U.A.)')

#Toroidal current
plt.subplot(2,2,2)
plt.plot(machine_data[:,0], machine_data[:,2])
plt.xlabel('time (s)')
plt.ylabel('Coil current (A)')

#Neutral gas pressure
plt.subplot(2,2,3)
plt.plot(machine_data[:,0], machine_data[:,3])
plt.xlabel('time (s)')
plt.ylabel('Pressure (Pa)')

#Heating power delivered by the micro-waves systems
plt.subplot(2,2,4)
plt.plot(machine_data[:,0], machine_data[:,4], label='LFS Power')
plt.plot(machine_data[:,0], machine_data[:,5], label='HFS Power')
plt.xlabel('time (s)')
plt.ylabel('Heating Power (W)')

#Add legend and save
plt.legend()
plt.show()
plt.savefig(f"{path}/Figures/machine_data/machine_data{shot}")

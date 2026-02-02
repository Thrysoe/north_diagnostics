"""
This program was designed to run with python 3 in a Spyder environnement, with usual packages plus npTDMS
and local libraries found in the folder diagnostics and utils of the environment.
Only the imput parameter section should be modified.

They are: - [shot] the number of the shot
          - [bias_type] the way to bias the probes in the given shot
              Note that it can only be a temperature or ion_saturation_current bias
              Other bias type will trigger an error message
          - [T_sweep] the period of the sweep of the bias voltage. 
              Only relevant for temperature mesurements
          - [studied_channels] the number of the probe collecting data in the given shot
          - [path_to_data] the name of the data folder with the full path to acceed it
          - [path_to_figure] the name of the figure folder with the full path to acceed it
          - [t_start] the start time of the data collection. 
              Must be between 0 and 1 second
          - [t_end] the end time of the data collection. 
              Must be between 0 and 1 second and greater than t_start


Extract all data from the DDAQ and CRIO files and plot it. Both .tdms file must be saved in the data folder,
but there is a possibility to change the path file at the beginning of the main program.
Save all useful data in a .txt file with a short header description, plot machine parameters and I-t and
U-t curves to see if all of the data collected aren't anomalous.
"""

#Find the path to local libraries (modify the PYTHON PATH)
import sys
sys.path.append('C:\\Users\\Saïd\\Downloads\\north_diagnostics\\diagnostics')
sys.path.append('C:\\Users\\Saïd\\Downloads\\north_diagnostics\\utils')
sys.path.append('C:\\Users\\Saïd\\Downloads\\north_diagnostics')

#Import all useful libraries
import matplotlib.pyplot as plt
import numpy as np
import os
import utils.dau

#Input parameters
shot = 10164
bias_type = 'ion_saturation_current'
T_sweep = 1/75
studied_probes = np.array([1,2,3,4,5,6,7,8,9,10,11,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,38,39,40,41,42,43,44,45,46,47,48,49])
path_to_data = 'E:/north_diagnostics/Data'
path_to_figure = 'E:/north_diagnostics/Figures'
t_start = 300E-3
t_end = 400E-3





"""
Main program: from now, all shall remain untouched.
"""
#Create a folder to store images if they don't exist
if not os.path.exists(f"{path_to_figure}/IandVplots_{shot}/"):
  os.makedirs(f"{path_to_figure}/IandVplots_{shot}/")
  print("A new directory for storing data was created")

#Generate the data and plot the curves
machine_data = utils.dau.read_machine_data(shot, path_to_data)
probe_data = utils.dau.read_probe_data(shot, path_to_data, bias_type, T_sweep, studied_probes-1, t_start, t_end, path_to_figure)
plt.close()

#Saving all data in the Data folder
head = 'Time; Light sensor; Coil current; Pressure sensor; LFS power; HFS power in SI units'
np.savetxt(f"{path_to_data}/machine_data{shot}.txt", machine_data, delimiter=';', header=head)
head = 'Time; probes in the numerical order in SI units'
np.savetxt(f"{path_to_data}/probe_data{shot}.txt", probe_data, delimiter=';', header=head)

#Plot and save machine parameters
fig = plt.figure(figsize=(10,10))
plt.subplot(2,2,1)
plt.plot(machine_data[:,0]*1E3, machine_data[:,1])
plt.xlabel('time (s)')
plt.ylabel('light sensor signal (U.A.)')

plt.subplot(2,2,2)
plt.plot(machine_data[:,0]*1E3, machine_data[:,2])
plt.xlabel('time (s)')
plt.ylabel('Coil current (A)')

plt.subplot(2,2,3)
plt.plot(machine_data[:,0]*1E3, machine_data[:,3])
plt.xlabel('time (s)')
plt.ylabel('Pressure (Pa)')

plt.subplot(2,2,4)
plt.plot(machine_data[:,0]*1E3, machine_data[:,4], label='LFS Power')
plt.plot(machine_data[:,0]*1E3, machine_data[:,5], label='HFS Power')
plt.xlabel('time (s)')
plt.ylabel('Heating Power (W)')

plt.legend()
plt.show()
plt.savefig(f"{path_to_figure}/machine_data{shot}")
plt.close(fig)

"""
Only the imput parameter section should be modified.

They are: - [shot] int
              the number of the shot
          - [path] str
              the name of the folder in which all the code architecture is 
          - [current_value] float
              the value of the current in the vertical coils, only relevant when [data_origin] is simulation
          - [data_origin] str
              the way data were collected. Can be experiment or simulation.
          - [NFFT] int
              the window on which the signal is sampled to compute the FFT
              With great NFFT comes great frequency resolution at the cost of time resolution
              Recommended values: 2048 experiment
              
Take the data from the data from the DDAQ file or simulation file.
Plot it into a PSD spectrogram of each probes to analyze the presence or abscence of probes.
Not adapted to the simulation data due to the too short duration.
"""

#Find the path to local libraries
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
shot = 11342
path = 'C:/Users/tulla/Perso/north_diagnostics'
current_value = 0
data_origin = 'experiment'
NFFT = 2048




"""
Main program: from now, all shall remain untouched.
"""
#Extract experiment or simulation data
if data_origin == 'experiment':
      #Activated probes for that type of data
      studied_probes = np.array([1,2,3,4,5,6,7,8,9,10,11,13,14,15,16,17,18,19,20,21,22,23,24,25,
                                 26,27,28,29,30,31,32,33,38,39,40,41,42,43,44,45,46,47,48,49])
      data = utils.dau.read_probe_data(shot, path, 'ion_saturation_current', 1/75, studied_probes-1, 
                                     0, 1, False)
      print(f"Data shot {shot} loaded with success")
    
      #Create a folder to store data if it doesn't exist
      if not os.path.exists(f"{path}/Figures/{shot}_spectrogram/"):
          os.makedirs(f"{path}/Figures/{shot}_spectrogram/")
          print("A new directory for storing data was created")
    
      #The time step between two recorded points
      time_step = 1E-6  
    
elif data_origin == 'simulation': 
      #Activated probes for that type of data
      studied_probes = np.array([1,2,3,4,5,6,7,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,
                                 26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50])
      
      all_data = np.genfromtxt(f"{path}/Data/probeIsat_fullmodel.txt", skip_header=5)
      print("Data probe loaded with success")
      
      #Take only the data for the given current value
      n_sweep = int((current_value)//20)
      data = all_data[n_sweep*5001:(n_sweep+1)*5001, :]
      
      #Create a folder to store data if it doesn't exist
      if not os.path.exists(f"{path}/Figures/simu_spectrogram_{int(current_value)}/"):
        os.makedirs(f"{path}/Figures/simu_spectrogram_{int(current_value)}/")
        print("A new directory for storing data was created")  
    
      #The time step between two recorded points
      time_step = 1E-6 
        
#Send an error message if their is a mispelling in the data_origin variable  
else:
    print('WARNING: the data origin is not recognized')
    sys.exit()

#Create figure 
fig = plt.figure(figsize=(10,10))
  
#Plot all images of the vessel + the probes + the data for each time
for i in studied_probes:
    output = utils.dau.plot_spectrogram_fft(data[:, i], shot, current_value, path, i, 
                                            data_origin, time_step, NFFT)
    print(output)
    
plt.close()

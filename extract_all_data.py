"""
Extract all data from the DDAQ and CRIO files and plot it. Both .tdms file must be saved in the data folder, but there is a possibility to change the path file 
at the beginning of the main program.
Usefull for further data treatment and save all usefull data in a .txt file.
Define shot number just after the main programm begins.
"""

#Import all useful libraries
import nptdms
import scipy.optimize import curve_fit
from north_diagnostics.diagnostics import Probe, Diagnostic
import matplotlib.pyplot as plt
import numpy as np

# Define all useful functions
def current_fit(U, I_isat, k_BT_e, U_f):
  """Function fitted by scipy in the transition function of the I-V curve for electron temperature measurements"""
  e0 = 1.60E-19 # in C
  return I_isat*(np.exp(e0*(U-U_f)/k_BT_e)-1)
  
def read_machine_data(shot, path_to_data):
  """
  Read the channels:  1 (Light sensor)
                      2 (Coil currents)
                      3 (Pressure sensor)
                      6 (LFS power)
                      7 (HFS power)
  to control machine parameters. Plotted in the main program and saved in a .txt files in the Data folder.
  """
  #Read data file
  file = nptdms.TdmsFile.read(f"{path_to_data}/CRIO{shot}.tdms")

  #Define the time interval of the study (the mask variable)
  t_start = 0
  t_end = machine_file['Data']['Time'][-1]
  mask = (machine_file['Data']['Time'][:] >= t_start) & (machine_file['Data']['Time'][:] <= t_end)

  #Extracting all machine parameters
  data[:,0] = t
  data = np.zeros((len(t), 6))
  head = 'Time; Light sensor; Coil current; Pressure sensor; LFS power; HFS power'
  data[:,0] = t*1E-3 # in s
  data[:,1] = machine_file['Data']['Light'][mask] # in ??; whatever, not for a quantitative analysis
  data[:,2] = machine_file['Data']['I_TF'][mask] # in A
  data[:,3] = machine_file['Data']['Pressure'][mask]*1E2 # in Pa
  data[:,4] = machine_file['Data']['LFSset'][mask]*450/3000 # in W
  data[:,5] = machine_file['Data']['HFSset'][mask]*450/3000 # in W
  return data

def read_probe_data(shot, path_to_data, m_i, A, T_sweep, k_B, e):
  """
  Read all the probes channels which are n_channel=n_probe+14. 
  Saved in a .txt files in the Data folder.
  """
  #Useful variables
  probe = {}
  data = np.zeros((len(t), 51))

  #Read data file
  for i in range(Probe.TOTAL_PROBES):
    probe[i] = Probe(path = path_to_data, shot = shot, number = i+1, caching = True)
    t = probe[i].time
    U = probe[i].bias_voltage
    I = probe[i].current
    if bias_type == 'density':
      if i==0:
        data = np.zeros((len(t), 51))
        data[:,0] = t
      data[:,i+1] = I/(0.61*e*np.sqrt(k_B*T_e_estim/m_i)*A)
    elif bias_type == 'temperature':
      if i==0:
        N = int(len(t)/T_sweep) #Number of temperature measurements possible, take a full sweep to be unbothered by hysterisis effects
        data = np.zeros(N, 51))
        data[:,0] = range(N)*T_sweep + Tsweep/2
      for j in range(N):
        start, end = get_time_indices(j*T_sweep, (j+1)*T_sweep)
        guess = [0.1, 1E-19, 10]
        popt, pcov = curve_fit(current_fit, U[start, end], I[start, end], guess)
        data[j,i+1] = popt[1]/k_B
    else:
      print('WARNING: the bias type is not recognized')
      exit()
  return data

#Main program: plot all machine and some probe data to verify that the shot "looks fine"
if __name__=="__main__":
  #Input parameters
  shot = 9974
  A_gas = 4 #Helium 4
  T_sweep = 13.33E-3 #period of the sweep for electronic temperature measurements in s
  bias_type = 'temperature' #Probes can be biased to measure 'density' or 'temperature' (the same bias is applied on every probe)
  path_to_data = './north_diagnostics/Data/'
  path_to_figure = './north_diagnostics/Figures/'

  #Physical constants
  e = 1.602E-19 # in C
  k_B = 1.38E-23 # in J/K
  m_i = A_gas*1.67E-27 # ion mass in kg

  #Experiment parameters
  A = 1E-6 # probe surface in m^2
  T_e_estim = 10*e/k_B # estimate of temperature in K (the numerical value is in eV) to measure electron density

  #Generate the data
  machine_data = read_machine_data(shot, path_to_data)
  probe_data = read_probe_data(shot, path_to_data, m_i, A, T_sweep, k_B, e)
  
  #Saving all data in the Data folder
  head = 'Time; Light sensor; Coil current; Pressure sensor; LFS power; HFS power in SI units'
  np.savetxt(f"{path_to_data}machine_data{shot}.txt", machine_data, delimiter=';', header=head)
  head = 'Time; probes in the numerical order in SI units'
  np.savetxt(f"{path_to_data}probe_data{shot}.txt", probe_data, delimiter=';', header=head)
  
  #Plot and save figures
  plt.subplot(2,2,1)
  plt.plot(machine_data[:,0]*1E3, machine_data[:,1])
  plt.xlabel('time (ms)')
  plt.ylabel('light sensor signal (U.A.)')
  
  plt.subplot(2,2,2)
  plt.plot(machine_data[:,0]*1E3, machine_data[:,2])
  plt.xlabel('time (ms)')
  plt.ylabel('Coil current (A)')
  
  plt.subplot(2,2,3)
  plt.plot(machine_data[:,0]*1E3, machine_data[:,3])
  plt.xlabel('time (ms)')
  plt.ylabel('Pressure (Pa)')
  
  plt.subplot(2,2,4)
  plt.plot(machine_data[:,0]*1E3, machine_data[:,4], label='LFS Power')
  plt.plot(machine_data[:,0]*1E3, machine_data[:,5], label='HFS Power')
  plt.xlabel('time (ms)')
  plt.ylabel('Heating Power (W)')
  
  plt.legend()
  plt.show()
  plt.savefig(f"{path_to_figure}machine_data{shot})

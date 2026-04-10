"""
Contains the functions relevant for data analysis in the NORTH tokamak. 
For specific details, see the documentation of each function.
"""

#Import libraries
import nptdms
import scipy.optimize as scopt
import diagnostics.probe, diagnostics.diagnostic
import matplotlib.pyplot as plt
import numpy as np
import scipy.interpolate as sci
import os
import cv2
import diagnostics.probe
import sys
import scipy.stats as scst
import scipy.signal as scs



#Definition of the functions
def current_fit(U, I_isat, k_BT_e, U_f):
  """
  Function fitted by scipy in the transition function of the I-V curve for electron 
  temperature measurements.
  Only k_BT_e is the electron temperature in J, other quantities must be given in SI.
  
  Input:  - [U] int
              the bias voltage
          - [I_isat] int
              the estimated ion saturation current
          - [k_BT_e] int
              the estimated electron temperature
          - [U_f] int
              the estimated floating potential
          
  Output: The value of the probe current I for a given bias volatge I
  """
  e0 = 1.602E-19 # in C
  return I_isat*(np.exp(e0*(U-U_f)/k_BT_e)-1)

def rolling_average(data, av_size):
    """
    Computes the rolling average of a 1D array. Can be used to smooth noisy data.
    
    Input:  - [data] numpy array
                the one dimensional data array to smooth
            - [av_size] int
                the length of the average segment
            
    Output: - [data_av_roll] numpy array
                the rolling average of the data
    """
    #Define the boundary of the domain on which rolling average can be done
    start_av = (av_size-1)//2
    end_av = len(data)+(av_size-1)//2
    
    #Compute the rolling average
    data_roll_av = np.array([np.mean(data[i-start_av:i+start_av]) for i in range(start_av, end_av)])
    return data_roll_av

def median_filter(data, med_size):
    """
    Computes the median of the points of a 1D array on a length med_size. Can be used to smooth noisy data.
    
    Input:  - [data] numpy array
                the one dimensional data array to smooth
            - [av_size] int
                the length of the segment on which the median will be calculated
            
    Output: - [data_av_roll] numpy array
                the median filter of the data
    """
    #Define the boundary of the domain on which rolling average can be done
    start_av = (med_size-1)//2
    end_av = len(data)+(med_size-1)//2
    
    #Compute the rolling average
    data_med_av = np.array([np.median(data[i-start_av:i+start_av]) for i in range(start_av, end_av)])
    return data_med_av

def read_machine_data(shot, path):
  """
  Read the channels:  1 (Light sensor)
                      2 (Coil currents)
                      3 (Pressure sensor)
                      6 (LFS power)
                      7 (HFS power)
  to control machine parameters. Plotted in the main program and saved in a .txt files in the Data folder.
  
  Input:  - [shot] int
              the number of the studied shot
          - [path] str
              the name of the folder in which all the code architecture is 
          
  Output: - [data] numpy array
              The data are stored into a 6 column table:
                  - 0: time of the experiment
                  - 1: the signal collected by the light sensor
                  - 2: the toroidal coil current
                  - 3: the gas pressure 
                  - 4: the heating power of the LFS system
                  - 5: the heating power of the HFS system
  """
  #Read data file
  file = nptdms.TdmsFile.read(f"{path}/Data/CRIO{shot}.tdms")

  #Define the time interval of the study (the mask variable)
  t_start = 0
  t_end = file['Data']['Time'][-1]
  mask = (file['Data']['Time'][:] >= t_start) & (file['Data']['Time'][:] <= t_end)

  #Extracting all machine parameters
  t = file['Data']['Time'][mask]
  data = np.zeros((len(t), 6))
  data[:,0] = t*1E-3 # in s
  data[:,1] = file['Data']['Light'][mask] # in ??; whatever, not for a quantitative analysis
  data[:,2] = file['Data']['I_TF'][mask] # in A
  data[:,3] = file['Data']['Pressure'][mask]*1E2 # in Pa
  data[:,4] = file['Data']['LFSset'][mask]*450/3000 # in W
  data[:,5] = file['Data']['HFSset'][mask]*450/3000 # in W
  return data

def read_probe_data(shot, path, bias_type, T_sweep, studied_probes, t_start, t_end, plot_curves):
  """
  Read all the activated probe channels which correspondance is given in the mappings file in the utils 
  folder. Saved in a .txt files in the Data folder.
  While data are saved, it plots all of the I-t and V-t curves.
  
  Input:  - [shot] int
              the number of the studied shot
          - [path] str
              the name of the folder in which all the code architecture is 
          - [bias_type] str
              the way used to bias the probes, which changes the quantity measured by it. Note 
              that it can only be a temperature or ion_saturation_current bias. Other bias type will 
              trigger an error message.
          - [T_sweep] float
              the period of the bias voltage sweep. Only relevant for temperature measurement.
          - [studied_channels] numpy array
              the list of the channels activated for the analysis
          - [t_start] float
              the start time of the data collection. Must be between 0 and 1 second.
          - [t_end] float
              the end time of the data collection. Must be between 0 and 1 second and greater 
              than t_start.
          - [plot_curves] boolean
              the boolean variable to control if we plot test curves or not
              
  Output: - [data] numpy array
              The data are stored into a len(studied_channels) + 1 column table:
                  - 0: time of the measurement
                  - i: the signal collected by the i Langmuir probe according to the selected bias type.
  """
  #Initialisation of probe variable
  probe = {}
  k_B = 1.38E-23 #the Boltzmann constant in SI

  for i in studied_probes: 
    #Read data file
    probe[i] = diagnostics.Probe(path = f"{path}/Data", shot = shot, number = i+1, caching = True)
    
    #Extracting all data from probe i in the DDAQ file
    if i==0:
      ind_start, ind_end = probe[i].get_time_indices(t_start, t_end)
      ind_Pon, ind_Poff = probe[i].get_time_indices(90E-3, 910E-3)
    t = probe[i].time[ind_start:ind_end] 
    U = probe[i].bias_voltage[ind_start:ind_end] 
    I = probe[i].current[ind_start:ind_end] - np.mean(probe[i].current[:ind_Pon])
    
    #Collect data
    if bias_type == 'ion_saturation_current':
      if i==0:
        data = np.zeros((len(t), 51))
        data[:,0] = t
      data[:,i+1] = I 
      
    elif bias_type == 'temperature':
      if i==0:
        #Number of temperature measurements possible
        N = int((t_end-t_start)/T_sweep)
        data = np.zeros((N, 51))
        data[:,0] = np.array(range(N))*T_sweep + T_sweep/2 + t_start
      #Fit the IV-curve to measure temperature
      for j in range(N):
        start, end = probe[i].get_time_indices(j*T_sweep, (j+1)*T_sweep)
        guess = [0.0001, 1E-18, -10]
        popt, pcov = scopt.curve_fit(current_fit, U[start:end], I[start:end], guess)
        data[j,i+1] = popt[1]/k_B
        
    elif bias_type == 'potential':
      if i==0:
        #Number of potential measurements possible
        N = int((t_end-t_start)/T_sweep)
        data = np.zeros((N, 51))
        data[:,0] = np.array(range(N))*T_sweep + T_sweep/2 + t_start
      #Fit the IV-curve to measure potential
      for j in range(N):
        start, end = probe[i].get_time_indices(j*T_sweep, (j+1)*T_sweep)
        guess = [0.0001, 1E-18, -10]
        popt, pcov = scopt.curve_fit(current_fit, U[start:end], I[start:end], guess)
        data[j,i+1] = popt[2]
        
    elif bias_type == 'sweep':
      if i==0:
        #Number of sweep measurements possible
        N = int((t_end-t_start)/T_sweep)
        data = np.zeros((N, 51, 6))
        data[:,0,0] = np.array(range(N))*T_sweep + T_sweep/2 + t_start
        data[:,0,1] = np.array(range(N))*T_sweep + T_sweep/2 + t_start
        data[:,0,2] = np.array(range(N))*T_sweep + T_sweep/2 + t_start
      #Smoothing of the data
      I_smooth = scs.savgol_filter(median_filter(I, 10), 100, 6)
      #Fit the IV-curve to measure Iisat, temperature and potential
      for j in range(N):
        #get the time indices for a given period
        start, end = probe[i].get_time_indices(j*T_sweep, (j+1)*T_sweep)
        
        #Take only the data before the electron saturation knee
        uj = U[start:end]
        ij = I_smooth[start:end]
        
        #Fit the curve
        guess = [-1E-4, 1E-18, -10]
        ind = [i for i in range(len(uj)) if uj[i] < 12]
        popt, pcov = scopt.curve_fit(current_fit, uj[ind], ij[ind], guess)
        
        #Save the data
        data[j, i+1, 0] = - popt[0]
        data[j, i+1, 1] = popt[1]/k_B
        data[j, i+1, 2] = popt[2]
        data[j, i+1, 3] = 2*np.sqrt(np.abs(pcov[0,0]/popt[0]))
        data[j, i+1, 4] = 2*np.sqrt(np.abs(pcov[1,1]/popt[1]))
        data[j, i+1, 5] = 2*np.sqrt(np.abs(pcov[2,2]/popt[2]))
        
    #Send an error message if there is a misspellin g the input parameters
    else:
      print('WARNING: the bias type is not recognized')
      sys.exit()
    
    #Plot the I-t and V-t curves if the plot_curves variable is True
    if plot_curves == True:
        if bias_type == 'ion_saturation_current':
            #Voltage bias curve
            plt.subplot(2,1,1)
            plt.plot(probe[i].time*1e3, probe[i].bias_voltage, color='k', linewidth=2.5)
            plt.xlim(0, 1000)
            plt.ylabel("$V_{\\rm bias}$ [V]", fontsize=14)
        
            plt.title(f"Shot {probe[i].shot} : Probe {probe[i].number}", fontsize=18)
        
            #Probe current curve
            plt.subplot(2,1,2)
            plt.plot(probe[i].time*1e3, (probe[i].current - np.mean(probe[i].current[:ind_Pon]))*1e3, color='k', linewidth=2.5)
            plt.xlim(0, 1000) 
            plt.xlabel("$t$ [ms]", fontsize=14)
            plt.ylabel("$I_{\\rm probe}$ [mA]", fontsize=14)
        
            plt.savefig(f"{path}/Figures/IandVplots_{shot}/probe{i+1}.png", dpi=300)
            plt.clf()
        else:
            #Plot the IV curve
            plt.plot(probe[i].bias_voltage, (probe[i].current - np.mean(probe[i].current[:ind_Pon]))*1e3, color='k', linewidth=2.5)
            plt.xlabel("$V_{\\rm bias}$ [V]", fontsize=14)
            plt.ylabel("$I$ [mA]", fontsize=14)
        
            plt.title(f"Shot {probe[i].shot} : Probe {probe[i].number}", fontsize=18)
            plt.savefig(f"{path}/Figures/IandVplots_{shot}/probe{i+1}.png", dpi=300)
            plt.clf()
  return data

def plot_2D_data(data, shot, path, bias_type, data_type, time, activated_probes, vmin, vmax, fig, bc, layout=True):
  """
  Takes the data from the txt files and probe position from the .json file in utils to build the vessel, 
  the probes and the image at a given time of raw or fluctuating data.
  
  Input:  - [data] numpy array
              the data array at time [time]
          - [shot] int
              the number of the studied shot
          - [path] str
              the name of the folder in which all the code architecture is 
          - [bias_type] str
              the way used to bias the probes, which changes the quantity measured by it. Note 
              that it can only be a temperature or ion saturation current bias. Other bias type will 
              trigger an error message.
          - [data_type] str
              the type of data plotted. Can be raw data or fluctuations (data minus mean value).
          - [time] float
              the time at which the image is made
          - [activated_probes] numpy array
              the list of the channels activated for the analysis
          - [v_min] float
              the minimum of the colorbar
          - [v_max] float
              the maximum of the colorbar
          - [fig] matplotlib.figure.Figure
              Create a unique figure to avoid being sorrounded by ploting window
          - [bc] float
              The boundary condition that must be applied at the walls
          - [layout] boolean
              If True (default value), the layout is plotted. Set False to have no layout
          
  Output: - [output] str
              The program plots the vessel, the activated probes in white (others in red) and the colormap 
              associated with the studied quantity and give a string character to notify that the plotting 
              procedure happened well. The legend of the image and its save is done out of the program to make it
              more polyvalent.
  
  """
  #Plot the vessel and the probes
  theta = np.linspace(0, 2*np.pi, 250)
  x_loc = 250 + 125*np.cos(theta)
  y_loc = 125*np.sin(theta)
  if layout == True:
      plt.plot(x_loc, y_loc, label='Vessel boundaries', color='black')
      
  #Plot the heating zone
  file = nptdms.TdmsFile.read(f"{path}/Data/CRIO{shot}.tdms")
  t_machine = np.array(file['Data']['Time'])*1E-3
  I_TF = file['Data']['I_TF'] # in A
  
  #Machine parameter for magnetic field calculation and resonnance layer location
  mu_0 = 4*np.pi*1E-7 # in H/m
  N_TF = 8
  N_winding = 12
  R0 = 250E-3 # in m
  f_R = 2.45E9 #MW frequency
  e = 1.6E-19 #Elementary charge
  me = 9.11E-31 #Mass of an electron

  #Convert current into the magnetic field and resonnance layer position
  ind = np.argmin(np.abs(t_machine-data[0]))
  B_0_tor = mu_0*(N_TF*N_winding)*I_TF[ind]/(2*np.pi*R0)
  P_dist = e*B_0_tor*R0/(me*2*np.pi*f_R)-R0 # in mm
  
  #Conversion to mm units and plot of the layer
  R0 = 250 # in mm
  P_dist = P_dist*1E3 # in mm
  sigma_x = 10 # in mm
  sigma_y = 62.5 # in mm
  x_abs = np.linspace(-2*sigma_x, 2*sigma_x, 250)
  y_up = 2*sigma_y*np.sqrt(1-(x_abs/(2*sigma_x))**2)
  y_down = -2*sigma_y*np.sqrt(1-(x_abs/(2*sigma_x))**2)
  if layout == True:
      plt.plot(R0+P_dist+x_abs, y_up, label='Heating zone (2$\sigma$)', color='blue', linestyle='dashed')
      plt.plot(R0+P_dist+x_abs, y_down, color='blue', linestyle='dashed')
    
  #Plot the probes and get their positions
  r, z = [], []
  for i in range(diagnostics.Probe.TOTAL_PROBES):
      probe = diagnostics.Probe(path = f"{path}/Data", shot = shot, number = i + 1, caching = True)
      x_p, y_p = probe.position['r'], probe.position['z']
      if activated_probes[i]:
          r.append(x_p)
          z.append(y_p)
      if layout == True:
          plt.plot(x_p, y_p, marker='o', markeredgecolor='k', markerfacecolor='w' if activated_probes[i] else 'r')
          plt.text(x_p, y_p-10, str(i+1), color='black', fontsize=8)
  loc_probe = np.column_stack((r, z))

  #Add all data at this time on the plot and plot only relevant data
  c = []
  for i in range(diagnostics.Probe.TOTAL_PROBES):
      if activated_probes[i]:
          c.append(data[i+1])
          
  #Enforce a Dirichlet boundary condition (like in the numerical computation) at the vessel walls
  boundary_coordinates = np.column_stack((x_loc,y_loc))
  boundary_condition = np.zeros(boundary_coordinates.shape[0])+bc #2.26e-05 on 1mm2
  all_loc = np.vstack((loc_probe, boundary_coordinates))
  all_values = np.concatenate((c, boundary_condition))

  #Creates the grid and the plot area
  grid_r, grid_z = np.meshgrid(np.linspace(250-125, 250+125, 100), np.linspace(-125, 125, 100))
  mask = (grid_r - 250)**2 + (grid_z - 0)**2 > 125**2

  #Computing the map for a contour data plot
  c_interpolator = sci.RBFInterpolator(all_loc, all_values) #, neighbors=6)
  grid_c = c_interpolator(np.column_stack((grid_r.ravel(), grid_z.ravel()))).reshape(grid_r.shape)
  grid_c[mask] = np.nan  # Set outside the circle to NaN

  #Creates the data plot with its colorbar
  if layout == True:
      contourf = plt.contourf(grid_r, grid_z, grid_c, levels = np.linspace(vmin, vmax, 50), extend='both', 
                              cmap = 'inferno_r')
      plt.colorbar(contourf, orientation='vertical', label=f"{data_type} {bias_type} SI")
    
      #Add a legend, save and close
      plt.axis('equal')
      plt.xlim((250-140, 250+140)) 
      plt.xlabel('r (mm)')
      plt.ylabel('z (mm)')
  else:
      contourf = plt.contourf(grid_r, grid_z, grid_c, levels = np.linspace(vmin, vmax, 50), extend='both', 
                              cmap = 'gray')
      plt.axis('equal')
  
  return f"image {str(time)} processed"

def video_2D(data, data_origin, shot, current_value, path, bias_type, data_type, fps, start_time):
  """
  Takes the images in the appropriate figure folder and concanate them into a video .avi file.
  This function assumes that the images have already been generated and saved by the previous function.
  Be careful that the cv2 package is sensitive to the name of the folders and do not tolerate any special
  character.
  PROBLEM: the movie created seems to shuffle images when they are concatenated.
  
  Input:  - [data] numpy array
              the full array of data collected from the .txt file with the header removed
          - [data_origin] str
              the way data were collected. 
              Can be experiment or simulation.
          - [shot] int
              the number of the studied shot
          - [current_value] float
              the value of the vertical current
          - [path] str
              the name of the folder in which all the code architecture is 
          - [bias_type] str
              the way used to bias the probes, which changes the quantity measured by it. Note 
              that it can only be a temperature or ion saturation current bias. Other bias type will 
              trigger an error message.
          - [data_type] str
              the type of data plotted. Can be raw data or fluctuations (data minus mean value).
          - [fps] int
              the frame rate of the saved movie 
          - [start] int
              the time in ms at which the movie begins
          
  Output: - [output] str
              The program plots the vessel, the activated probes in white (others in red) and the colormap 
              associated with the studied quantity and give a string character to notify that the plotting 
              procedure happened.
  """
  #Create figure 
  if data_origin=='experiment':
      image_folder = f"{path}/Figures/{shot}_{bias_type}_{data_type}/"
      video_name = f"{path}/Figures/{shot}_{bias_type}_{data_type}/{shot}_{bias_type}_{data_type}_{start_time}.avi"
  else:
      image_folder = f"{path}/Figures/{int(current_value)}_{bias_type}_{data_type}/"
      video_name = f"{path}/Figures/{int(current_value)}_{bias_type}_{data_type}/simu_{int(current_value)}_{bias_type}_{data_type}_{start_time}.avi"
  
  images = [img for img in os.listdir(image_folder) if img.endswith((".jpg", ".jpeg", ".png"))]
  images.sort()

  # Set frame from the first image
  frame = cv2.imread(os.path.join(image_folder, images[0]))
  height, width, layers = frame.shape

  # Video writer to create .avi file
  video = cv2.VideoWriter(video_name, cv2.VideoWriter_fourcc(*'DIVX'), fps, (width, height))

  # Appending images to video
  for image in images:
    video.write(cv2.imread(os.path.join(image_folder, image)))

  # Release the video file
  video.release()
  cv2.destroyAllWindows()
  return f"Video is generated in the {image_folder} folder"

def frequency_fft(data, time_step):
  """
  Takes the data from one probe at all times to calculate the modulus of the temporal FFT of the signal.
  
  Input:  - [data] numpy array
              the 1D data array
          - [time_step] float
              the physical time step between two points in the data array
          
  Output: - [freq] numpy array
              The array of positive frequencies associated with the FFT.
          - [spectrum] numpy array
              The modulus of the FFT spectrum of the signal.
  """
  #Computes the FFT
  t_spectrum = np.fft.fft(data)
  
  #Computes the associated frequencies
  n = data.size
  f = np.fft.fftfreq(n, d = time_step)
  
  freq, spectrum = f[0:int(n/2-1)], np.abs(t_spectrum[0:int(n/2-1)])
  return freq, spectrum

def spatial_fft(data, space_step):
  """
  Takes the data from all the probes at a given time to plot 2D Fourier transform
  
  Input:  - [data] numpy array
              the 1D data array
          - [space_step] float
              the physical space step between two points in the data array
          
  Output: - [k] numpy array
              The array of positive frequencies associated with the FFT.
          - [spectrum] numpy array
              The modulus of the FFT spectrum of the signal.
  """
  #Computes the FFT
  xy_spectrum = np.fft.fft2(data)
  
  #Computes the associated frequencies
  n, m = np.shape(data)
  kx = np.fft.fftfreq(n, d = space_step)
  ky = np.fft.fftfreq(m, d = space_step)
  
  k_x, k_y, spectrum = kx[0:int(n/2-1)], ky[0:int(n/2-1)], np.abs(xy_spectrum[0:int(n/2-1), 0:int(n/2-1)])
  return k_x, k_y, spectrum

def plot_spectrogram_fft(data, shot, current_value, path, studied_probe, data_origin, time_step, NFFT):
  """
  Takes the data from one probe at all times to plot a PSD colormap.
  
  Input:  - [data] numpy array
              the 1D data array
          - [shot] int
              the number of the studied shot
          - [current_value] float
              the value of the current in the vertical coils
          - [path] str
              the name of the folder in which all the code architecture is 
          - [studied_probe] numpy array
              the number of the probe on which FFT is done
          - [data_origin] str
              the way data were collected. Can be experiment or simulation
          - [time_step] float
              the physical time step between two points in the data array
          - [NFFT] int
              the window on which the signal is sampled to compute the FFT
              With great NFFT comes great frequency resolution at the cost of time resolution
          
  Output: - [output] str
              The program plots the PSD colormap associated with the studied quantity and give a string character 
              to notify that the plotting procedure happened.
  """
  #Computes the f-t spectrogram and its associated colorbar
  spectrum, freqs, t, im = plt.specgram(data, Fs=1/time_step, cmap='inferno', mode='psd', 
                                        scale='dB', NFFT=NFFT, noverlap=NFFT//2)
  plt.colorbar(im, orientation='vertical', label="ion saturation current SI")

  #Add a legend, save and clear figure
  plt.xlabel('time (s)')
  plt.ylabel('frequency (Hz)')
  #plt.yscale("log")
  plt.ylim((10,1/(2*time_step)))
  
  if data_origin == 'experiment':
      plt.title(f"Shot {shot} ion saturation current PSD spectrogram probe {studied_probe}")
      plt.savefig(f"{path}/Figures/{shot}_spectrogram/probe{studied_probe}.png")

  elif data_origin == 'simulation': 
      plt.title(f"Ion saturation current PSD spectrogram probe {studied_probe}")
      plt.savefig(f"{path}/Figures/simu_spectrogram_{int(current_value)}/probe{studied_probe}.png")
  
  plt.clf()
  return f"Probe {studied_probe} processed"

def stat_analysis(data, shot, current_value, path, data_origin, studied_probe, bias_type, fig, k):
  """
  Takes the data from all the probes to plot their distribution function. 
  Then plot a map of the fourth first statistic moments. 
  In project, not yet written.
  
  Input:  - [data] numpy array
              the 1D data array
          - [shot] int
              the number of the studied shot
          - [current_value] float
              the value of the current in the vertical coils
          - [path] str
              the name of the folder in which all the code architecture is 
          - [data_origin] str
              the way data were collected. Can be experiment or simulation
          - [studied_probe] numpy array
              the number of the probe on which FFT is done
          - [bias_type] str
              the way used to bias the probes, which changes the quantity measured by it. Note 
              that it can only be a temperature or ion saturation current bias. Other bias type will 
              trigger an error message.
          - [fig] matplotlib.figure.Figure
              Create a unique figure to avoid being sorrounded by ploting window
          - [k] int
              the extent of the colorbar in terms of standard variation of the signal
          
  Output: - [m, sigma, skew, kurt] float tupple
              The program plots the PDF of each probe for the given quantities and calculates the first 4th order 
              moments. Note that the skewness and kurtosis are normalized by the standard deviation and that the
              kurtosis is calculated with Fisher method (kurt_{fat tail}>0; kurt_{light tail}<0)
  """
  #Compute the moments
  m = np.mean(data)
  sigma = np.std(data)
  skew = scst.skew(data)
  kurt = scst.kurtosis(data)
  
  #Plot the probability density function
  plt.hist((data-m)/sigma, bins=25, range=(-k, k), density=True)
  plt.xlabel(f'{bias_type} centered reduced')
  lab = [str(i)+r'$\sigma$' for i in range(-k, k+1)]
  lab[k-1] = r'$-\sigma$'
  lab[k] = r'$0$'
  lab[k+1] = r'$\sigma$'
  plt.xticks(range(-k, k+1), labels=lab)
  plt.ylabel("Probability density function centered reduced")
  plt.title(f'Probability density function of {bias_type} of probe {studied_probe}')
  
  #Save data
  if data_origin == 'experiment':
      plt.title(f"Shot {shot} ion saturation current probability density function of {bias_type} of probe {studied_probe}")
      plt.savefig(f"{path}/Figures/{shot}_stat/probe{studied_probe}.png")

  elif data_origin == 'simulation': 
      plt.title(f"Shot {shot} ion saturation current probability density function of {bias_type} of probe {studied_probe}")
      plt.savefig(f"{path}/Figures/simu_stat_{int(current_value)}/probe{studied_probe}.png")
  
  plt.clf()
  
  print(f"Probe {studied_probe} processed")
  return m, sigma, skew, kurt

 # Programme principal
if __name__ == " __main__ ":
    pass

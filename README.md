**To Download**: Go to [releases page](https://github.com/nelsongt/mfiaDLTS/releases), current version is v0.9.1-beta. Master branch is not always in a working state.

Introduction
------------

This is a collection of tools I have developed over the course of my PhD for performing Deep Level Transient Spectroscopy (DLTS) using a Zurich Instruments MFIA. The tools are not in the most user friendly state, but this code is the most complete open-source DLTS implementation that I am aware of. I will continue to update the tools as they are actively used at the RIT lab.

I developed my own code and built a custom DLTS setup after using two commercial DLTS solutions. I settled on the MFIA, which did not have any publically available (there is labview code used at IFPAN) DLTS software specifically built for it,  as my capacitance meter. The primary drive for this work was to build an affordable, modern DLTS system with the power and flexibility that comes with having one's own source code. For example, with minimal modification one could use this code to automatically perform DLTS experiments with different pulse widths without having to repeat the temperature scan. Eg. on each temp. step, measure a transient then change pulse width and measure again, repeatedly, then move on to next temp. step. In this way, the activation energy of the trap cross-section can be extracted without much effort. Similarly, one could do depth-profiling by changing the pulse height in an automated fashion. The possibilities are almost endless because of the variety of DLTS measurements that exist and the many niche requirements of particular samples.

More about the development of these tools and some of the results can be read about in my dissertation: https://scholarworks.rit.edu/theses/10100/

There is a DLTS introduction also included on this github site, see DLTSBackground.pdf. I plan to continue to update this DLTS theory manual that is to accompany the code. I will also update with a wiring diagram, and some pdf presentations on how to use the software.

The data files are compatible with the LDLTS software database, so that Laplace DLTS can be performed if required (requires purchasing a license to use the algorithms): http://info.ifpan.edu.pl/Dodatki/WordPress/laplacedlts/

-George Nelson, 2019 (gtn1425@rit.edu)


Requirements
------------

  -Device to perform DLTS on
  
  -Cryostat with electrical leads for device connection
  
  -Lakeshore temperature controller to control cryostat temperature
  
  -Zurich Instruments MFIA
  
  -PC with MATLAB to run this code


Organization of Code
------------

In this early stage, the code is organized into 3 folders:

AcquireData - This is code to gather data for both Conventional DLTS (CDLTS) and Admittance Spectroscopy (YSpec). I will not explain the YSpec code at this time. The CDLTS code is designed to set the temperature via a lakeshore, then communicate with the MFIA to setup the DLTS experiment and then record the capacitance vs. time data. The CDLTS program outputs an .iso file for each temperature setpoint containing the digitized average transient data.

ProcessData - These tools will process the CDLTS and YSpec data files. The data folders generated by AcquireData tools should be moved to this folder, then point Transient_To_CDLTS.m file to this data folder containing all the .iso files. This tool takes the .iso transient data files and generates the conventional spectra.

SimulateData - This is a basic simulation tool to generate .iso files with transient data. Useful for debugging or predicting results or learning how DLTS works.

Software Dependencies
------------

General:

  MATLAB (with statistics and instrument control toolboxes)
  
ProcessData:

  [optional,not included] ezyfit:  http://www.fast.u-psud.fr/ezyfit/
  
AcquireData:

  [required, not included] LabOne: https://www.zhinst.com/downloads
  
  [required, not included] LabOne Matlab Driver (place in AcquireData folder): https://www.zhinst.com/downloads
  
  [required, included] lakshore driver: https://www.mathworks.com/matlabcentral/fileexchange/48366-lakeshore
  
  [required, included] cprintf: https://www.mathworks.com/matlabcentral/fileexchange/24093-cprintf-display-formatted-colored-text-in-the-command-window
  
How To Use
------------
Updated for V0.9

#### AcquireData

To take data, make sure that the software dependencies are installed correctly. You will then need to wire the MFIA properly, for this I will make a wiring diagram in the future. In the meantime:

-Connect sample leads to MFIA in 2-wire mode.

-Connect a BNC cable from Aux Out 1 of MFIA to Aux Input 1. This is for pulse generation. See: https://www.zhinst.com/blogs/tim/2018/08/22/square-pulse-c-v-measurements-for-dlts-on-the-mfia/

-Connect a BCN cable from Aux Out 2 to Trigger In 1 (on back of MFIA). This is for hardware triggering. See: https://www.zhinst.com/blogs/tim/2018/12/19/gated-data-transfer-for-increased-data-sampling-rate-on-the-mfia-and-mfli/

-Connect ethernet cable from MFIA to PC. PC and cable must be gigabit or better.

-Open CDLTS_Main.m. In this file, you will setup the experiment variables like temperature range and temperature step and sample biasing and pulsing. You will also describe the sample with a name and other key parameters.

Once setup, start the program and the output will make sure the MFIA and the lakeshore are working and configure them. After that, the temperature will be stabilized at the first temp step. Once the temperature is stable, the MFIA will be asked for data using its DAQ module. The code is setup to record data using a hardware trigger, where the voltage pulse is used as the trigger. Individual transients are recorded at the 'sampling rate' for a length of time determined by the 'individual transient sampling period'. The total number of transients recorded is determined by the 'sampling time', which is the total time you want to take data for each temperature point. Then, however many transients collected during that 'sampling time' are averaged. Finally, the averaged data is saved to disk before moving to the next temperature.

The scan process in time is: Temp set point ----- wait for temperature of cryostate to reach set point -----  temperature set wait time ----- record transients for sampling time ---- write file to disk ---- new temp set point, start over

*This code was initially inspired by software on the Zurich Instruments blogs.

#### ProcessData

After CDLTS_Main.m has collected all .iso files, take the generated data folder and move it from AcquireData to ProcessData. Then, open Transient_To_CDLTS.m and point the data folder to the folder that you just moved. It will go into that folder and process all the .iso files to generate a spectrum.

There are many complicated options in this file, such as choice of weighting function, that can be used by experts to tailor their processing. By default, the rectangular (lock-in) weighting function is used with emission rate constants from 20 Hz to 2000 Hz.

To extract the DLTS parameters, the matlab data cursor can be used on the plotted spectra to find the peaks. These peaks should be used to generate an Arrhenius plot. There is primitive support for doing this automatically, but it currently does not work well. You will have to generate most Arrhenius plots manually. This process can be understood by reading DLTSBackground.pdf.

#### SimulateData

To simulate DLTS, Transient_Sim.m can be used to generate transients for any number of traps. The trap characteristics (energy, cross-section, density) can be inputted in the simulator and these values should be recovered after processing the spectra. To create the conventional spectra from the simulated transients, use Transient_To_CDLTS.m. Take the generated transient files and place them in a folder in the ProcessData folder, then point Transient_To_CDLTS.m to that folder. Currently, the automated peak fitting is primitive, but creating the Arrhenius plots from the plotted spectra is simple to do manually.

*Initially developed by George Nelson of NanoPower Research Labs at Rochester Institute of Technology to be used in experiments funded by AF AEDC under SBIR# F151-173-0753*






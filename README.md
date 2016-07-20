# EZ_Temperature_Viewer
Allows viewing and manipulation of Luna ODiSI-B created ASCII temperature files for fiber optic temperature sensing 

#############################################################################

I. Installation Guide

Verified to work with (but not necessarily limited to): 
-	Matlab 2014a and MCR V8.3
-	Windows 7 – 64bit OS
Note: Checked with Matlab 2012b -> 
	- at least one issue because Matlab built-in function strsplit is a post (at least) 2012b addition...
	- could create alternative strsplit function or alternative fix in code
	
Decide which installation option is desired and follow guide:
a.	Use within Matlab:
	i.	Launch Matlab.
	ii.	Put EZ_Temperature_Viewer.m and EZ_Temperature_Viewer.fig files in Matlab path.
	iii.	Launch GUI by typing EZ_Temperature_Viewer in command prompt.

b.	Run from Executable (no local installation):
	i.	Requires Matlab Compiler Runtime V8.3 (MCR) to be installed on the computer.
		•	If error ‘could not find version 8.3 of the MCR’, download MCR from Mathworks website or see Option C.
	ii.	Double click the EZ_Temperature_Viewer Application file to launch
		•	Note: Opening and saving files within program will default to directory of the program.

c.	Install Application Locally with MCR (if MCR is already installed it will only install EZ_T*_V*):
	i.	Double click EXE_with_MCR_download (MCR will be automatically downloaded from the Mathworks website – file size ~1 MB) or ‘_included’ (MCR included is installer – file size ~0.6 GB). 
	ii.	Follow installation instruction. Default locations should be fine for most users.
	iii.	Locate program in installed directory and launch to use.

#############################################################################

II. User Guide:

1.	Launch EZ_Temperature_Viewer according to the chosen installation method.
2.	Press ‘Load .txt File’. Navigate to text file and follow instructions to load file. The file should output a tab delimited file with two arrays (data and textdata).
3.	Set the bulk temperature to the desired option. If loading a text file choose an Nx2 text file where the 1st column is time starting at 0 and the 2nd column is bulk temperature (°C). No headerlines in text file.
4.	Update time, temperature limits, and location limits as appropriate.
5.	Run the program. The default option is for the program to plot the data. A ‘relative time’ (time from beginning of the file) and absolute time indicators are located below the plots.
6.	Activate waterfall plots, save modified text files, or a movie as desired. Note that if the movie selection is chosen, it requires the data to be plotted and the program must be ‘in focus’. The program will first create waterfall plots, then prompt the user for the name of the modified text file and then create the file. Next the user will be prompted for the .avi file name if this option is chosen and then create the movie while updating the plot.
7.	To view the movie, click ‘Launch Movie Player’ and locate desired movie file.
8.	To view modified text file, repeat this guide at Step 2.
9.	To log event data from temperature trips, enable the event data log and run program. Tables can then be exported to tab deliminated files for viewing in Excel, Matlab, or other text editor. Push the ‘Plot Event Log’ to view the times at which the events occurred.

#############################################################################

III. Capabilities:

-	Load ODiSI-B created ASCII files of any size, recorded frequency, fiber length
-	Ability to control location in position in time and space for plotting, filtering, and output file creation. This allows the user to create small files quickly for events of interest for further analysis and viewing
-	Can pause/play and stop program during plotting and creating a movie
-	Create waterfall plots at user specified frequencies
-	Save ASCII (.txt) files of the data of interest. The output format is of the same form as the ODiSI-B output. This allows one to save their data then load in the new file for additional viewing and editing
-	Create .avi files for presentations, additional editing, faster access, etc.
-	Program automatically detects the maximum temperatures and time length of loaded files for user reference and automatically sets the location limits of files that have been modified
-	Ability to scale the raw differential data according to a constant bulk temperature, a time dependent bulk temperature (i.e. 2 column text file: 1st column with time starting at 0 and 2nd column temperature in °C), or to leave the temperature as a differential (delta) temperature
-	Bulk temperature text file does not need to be of the same length as the data file. The program will linearly interpolate the temperature based on time. If the bulk text file is shorter in time than the data file (e.g. t_bulk(end) = 10 seconds t_data(end) = 15 seconds) all data at times >= to the last bulk time will be scaled by the last bulk temperature and the user will be notified for their reference
-	Launch a video player from within the program for .avi playback
-	An indicator in the top left corner of the program lets the user know the program status
-	Ability to locate events between two temperature limits, export the events to a tab deliminated ASCII, and create a time plot of the events for quick understanding of the event occurrence.

#############################################################################

IV. Limitations:

-	Upon loading a file, the displayed file name is taken from the first headerline of the ASCII file, not the displayed filename (if the file name was changed manually) when in an explorer menu, Matlab path, etc. If the file name needs to be changed after creation, the USER may change the file name normally and then open the ASCII file in a text editor and update the file name in the first headerline as well
-	While the program can load any size file, the true limit of limit depends on the computer
-	Creating waterfall plots requires the program to store potentially significant amounts of data, therefore it is advisable to either decrease frequency of plotting or axial locations of interest if the operation is slow
-	Maximum and minimum temperatures of waterfall plot are limited to the USER specified temperature limits. All values in data range of interest falling outside of these limits are set to the temperature limit just for the waterfall plot section. This is necessary to allow the graphs to display properly.
-	It is advisable to create modified .txt files from shorter lengths of time and/or distance as it can take a significant amount of time to create the file and the only way to abort the operation is to force close the program from the Task Manager
-	Creating a movie requires the program to be ‘in focus’ as the program actually creates a plot than captures an image of the screen. This is slower than simply viewing the plots but allows very quick viewing of the entire file once created
-	Location limits will automatically be updated when the file is loaded but the user must click on one of the boxes and press enter to update the plots
-	The program must be run every time a new action (e.g. plot new data range) is to be generated (includes data extra analysis).

#############################################################################

V. Future Modification Suggestions

-	Within the .m file there are instructions of adding some sort of data analysis or filter section. The instructions/program are designed to allow the user to create their own self-contained function that can be simply added to the Matlab path. Plug and Play if you will. One major requirement is that if the user wants to work on the temperature data and create a new temperature array, the array must still be a MxN matrix where M is specified by the length of the Relative_Time matrix (Mx1) and N is specified by the Position matrix (1xN). Alternative options would be to simply create new plots of whatever outputs the user desires (e.g. frequency based analysis) within the user function such that they are completely independent of the program. Other suggestions and ideas are located in the .m file in the ‘Filter Data Section’
-	Important: The ‘Filter Data Section’ is purposefully located before any plotting, text saving, etc. such that any filters or analysis will be applied to subsequent plots, modified text files, movies, event logs, etc.

#############################################################################
Contact GitHub API Training Shop Blog About
© 2016 GitHub, Inc. Terms Privacy Security Status Help

function varargout = EZ_Temperature_Viewer(varargin)
% EZ_TEMPERATURE_VIEWER M-file for EZ_Temperature_Viewer.fig
%
% Introduction to EZ_Temperature_Viewer MATLAB GUI
%
% The following functions control the operation of the EZ_Temperature_Viewer.fig file (GUI). 
%
% User Guide:
%	1. Launch EZ_Temperature_Viewer by typing Ez_Temperature_Viewer into the command prompt.
%      EZ_Temperature_Viewer.m and EZ_Temperature_Viewer.fig files must be in the Matlab path.
%
%	2. Press ‘Load .txt File’. Navigate to text file and follow instructions to load file.
%	   The file should output a tab delimited file with two arrays (data and textdata)
%
%	3. Set the bulk temperature to the desired option. If loading a text file choose an 
%	   Nx2 text file where the 1st column is time starting at 0 and the 2nd column is 
%      bulk temperature (°C). No headerlines in text file
%
%	4.	Update time, temperature limits, and location limits as appropriate
%
%	5.	Run the program. The default option is for the program to plot the data. 
%       A ‘relative time’ (time from beginning of the file) and absolute time indicators 
%       are located below the plots
%
%	6.	Activate waterfall plots, save modified text files, or a movie as desired. Note 
%       that if the movie selection is chosen, it requires the data to be plotted and the
%       program must be ‘in focus’. The program will first create waterfall plots, then 
%       prompt the user for the name of the modified text file and then create the file. 
%       Next the user will be prompted for the .avi file name if this option is chosen and 
%       then create the movie while updating the plot.
%
%	7.	To view the movie, click ‘Launch Movie Player’ and locate desired movie file.
%
%	8.	To view modified text file, repeat this guide at Step 2.
%
%   9.	To log event data from temperature trips, enable the event data log and run program.
%       Tables can then be exported to tab deliminated files for viewing in Excel, Matlab,
%       or other text editor. Push the ‘Plot Event Log’ to view the times at which the events
%       occurred.
%
% Capabilities:
%	- Load ODiSI-B created .ascii files of any size, recorded frequency, fiber length
%	- Ability to control location in position in time and space for plotting, filtering,
%     and output file creation. This allows the user to create small files quickly for 
%     events of interest for further analysis and viewing
%	- Can pause/play and stop program during plotting and creating a movie
%	- Create waterfall plots at user specified frequencies
%	- Save .ascii (.txt) files of the data of interest. The output format is of the same
%     form as the ODiSI-B output. This allows one to save their data then load in the new
%     file for additional viewing and editing
%	- Create .avi files for presentations, additional editing, faster access, etc.
%	- Program automatically detects the maximum temperatures and time length of loaded
%     files for user reference and automatically sets the location limits of files that
%     have been modified
%	- Ability to scale the raw differential data according to a constant bulk temperature,
%     a time dependent bulk temperature (i.e. 2 column text file: 1st column with time 
%     starting at 0 and 2nd column temperature in °C), or to leave the temperature as a 
%     differential (delta) temperature
%	- Bulk temperature text file does not need to be of the same length as the data file.
%     The program will linearly interpolate the temperature based on time. If the bulk text
%     file is shorter in time than the data file (e.g. t_bulk(end) = 10 seconds 
%     t_data(end) = 15 seconds) all data at times >= to the last bulk time will be scaled by
%     the last bulk temperature and the user will be notified for their reference
%   - Launch the built in MATLAB movie viewer (if the image toolbox is available) or 
%     Windows media player for video playback
%	- An indicator in the top left corner of the program lets the user know the program status
%   - Ability to locate events between two temperature limits, export the events to a tab 
%     deliminated ASCII, and create a time plot of the events for quick understanding of the
%     event occurrence.
%
% Limitations:
%   - Upon loading a file, the displayed file name is taken from the first headerline of the 
%     ASCII file, not the displayed filename (if the file name was changed manually) when in 
%     an explorer menu, Matlab path, etc. If the file name needs to be changed after creation,
%     the USER may change the file name normally and then open the ASCII file in a text editor 
%     and update the file name in the first headerline as well.
%   - While the program can load any size file, the true limit of limit depends on the computer
%   - Creating waterfall plots requires the program to store potentially significant amounts
%     of data, therefore it is advisable to either decrease frequency of plotting or axial 
%     locations of interest if the operation is slow
%   - It is advisable to create modified .txt files from shorter lengths of time and/or distance
%     as it can take a significant amount of time to create the file and the only way to abort 
%     the operation is to force close the program from the Task Manager
%   - Creating a movie requires the program to be ‘in focus’ as the program actually creates a 
%     plot than captures an image of the screen. This is slower than simply viewing the plots
%     but allows very quick viewing of the entire file once created
%   - Location limits will automatically be updated when the file is loaded but the user must 
%     click on one of the boxes and press enter to update the plots
%   - The program must be run every time a new action (e.g. plot new data range) is to be 
%     generated (includes data extra analysis).
%
% Future Modification Suggestions
%	- Within the .m file there are instructions of adding some sort of data analysis or filter
%	  section. The instructions/program are designed to allow the user to create their own 
%	  self-contained function that can be simply added to the Matlab path. Plug and Play if 
%     you will. One major requirement is that if the user wants to work on the temperature data
%     and create a new temperature array, the array must still be a MxN matrix where M is 
%     specified by the length of the Relative_Time matrix (Mx1) and N is specified by the 
%     Position matrix (1xN). Alternative options would be to simply create new plots of whatever
%     outputs the user desires (e.g. frequency based analysis) within the user function such that
%     they are completely independent of the program. Other suggestions and ideas are located in
%     the .m file in the ‘Filter Data Section’
%	- Important: The ‘Filter Data Section’ is purposefully located before any plotting, 
%	  text saving, etc. such that any filterings or analysis will be applied to subsequent 
%	  plots, modified text files, movies, etc. 
%
%      EZ_TEMPERATURE_VIEWER, by itself, creates a new EZ_TEMPERATURE_VIEWER or raises the existing
%      singleton*.
%
%      H = EZ_TEMPERATURE_VIEWER returns the handle to a new EZ_TEMPERATURE_VIEWER or the handle to
%      the existing singleton*.
%
%      EZ_TEMPERATURE_VIEWER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EZ_TEMPERATURE_VIEWER.M with the given input arguments.
%
%      EZ_TEMPERATURE_VIEWER('Property','Value',...) creates a new EZ_TEMPERATURE_VIEWER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before EZ_Temperature_Viewer_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to EZ_Temperature_Viewer_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help EZ_Temperature_Viewer

% Author: M. Scott Greenwood
%
% Thanks to Roland Pfister for dlmcell function

% Last Modified by GUIDE v2.5 06-Aug-2014 17:19:31

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @EZ_Temperature_Viewer_OpeningFcn, ...
                   'gui_OutputFcn',  @EZ_Temperature_Viewer_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

%=========================================================================%

% --- Executes just before EZ_Temperature_Viewer is made visible.
function EZ_Temperature_Viewer_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to EZ_Temperature_Viewer (see VARARGIN)

%Centers GUI on screen
movegui(handles.Main_Figure,'center');

%Initialize plot appearance. This will change based on USER inputs
ymin = str2double(get(handles.y_min_input,'String'));
ymax = str2double(get(handles.y_max_input,'String'));
lp1_start = str2double(get(handles.Loop_1_Pos_Start,'string'));
lp1_stop  = str2double(get(handles.Loop_1_Pos_Stop,'string'));
lp2_start = str2double(get(handles.Loop_2_Pos_Start,'string'));
lp2_stop  = str2double(get(handles.Loop_2_Pos_Stop,'string'));

axes(handles.Loop_1_Temp_Plot);
title('Loop 1 Temperature Plot');
axis([lp1_start lp1_stop ymin ymax]);
xlabel('Location [m]');
ylabel('Delta Temperature [°C]');

axes(handles.Loop_2_Temp_Plot);
title('Loop 2 Temperature Plot');
axis([lp2_start lp2_stop ymin ymax]);
xlabel('Location [m]');
ylabel('Delta Temperature [°C]');

%Makes inactive Stop/Pause functions until program runs
set(handles.Stop_Button,'Enable','inactive','BackgroundColor', [0.8 0.8 0.8]);
set(handles.Pause_Button,'Enable','inactive','BackgroundColor', [0.8 0.8 0.8]);

% Choose default command line output for EZ_Temperature_Viewer
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes EZ_Temperature_Viewer wait for user response (see UIRESUME)
% uiwait(handles.Main_Figure);

%=========================================================================%

% --- Outputs from this function are returned to the command line.
function varargout = EZ_Temperature_Viewer_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

%=========================================================================%

% --- Executes on button press in Browse_Button.
function Browse_Button_Callback(hObject, eventdata, handles)

set(handles.Status_Indicator,'String','Loading File');
set(handles.Warning_Msg,'String',{});

set(handles.y_min_input,'Enable','off');
set(handles.y_max_input,'Enable','off');
set(handles.Loop_1_Pos_Start,'Enable','off');
set(handles.Loop_1_Pos_Stop,'Enable','off');
set(handles.Loop_2_Pos_Start,'Enable','off');
set(handles.Loop_2_Pos_Stop,'Enable','off');
set(handles.Time_Start,'Enable','off');
set(handles.Time_Stop,'Enable','off');
set(handles.Run_Prog_Button,'Enable','inactive','BackgroundColor',[0.8 0.8 0.8]);
set(handles.Filter_Toggle,'Enable','off');
set(handles.Plot_Toggle,'Enable','off');
set(handles.Mod_txt_Toggle,'Enable','off','Value',0);
set(handles.Save_Vid_Toggle,'Enable','off','Value',0);
set(handles.Waterfall_Toggle,'Enable','off','Value',0);
set(handles.Browse_Button,'Enable','off');
set(handles.Freq_Wplot,'Enable','off');
set(handles.Temp_bulk_browse,'Enable','off');
set(handles.Constant_Bulk_Ref,'Enable','off');
set(handles.Launch_MP_Button,'Enable','off');

S = uiimport;   %Opens GUI for user selection of text file
if isempty(S)
    %USER clicked cancel. Restore defaults and stop run.
    EZ_Temperature_Viewer('Stop_Button_Callback',hObject,[],handles);
    return;
end

%Reset Event Data Tables
row_names = {'Event #1';'Event #2'};
col_names = {'Time [sec]','Pos #1 [m]'};
table_data = [];
set(handles.Loop_1_Event_Log,'RowName',row_names,...
    'ColumnName',col_names,...
    'Data',table_data);
set(handles.Loop_2_Event_Log,'RowName',row_names,...
    'ColumnName',col_names,...
    'Data',table_data);
                         
%If an error occurs (most likely due to incorrect loaded .txt format) the
%program will jump to 'catch' and stop
try
    
    %Define variables from text file
    [rtext, ~] = size(S.textdata);
    [rdata, ~] = size(S.data);
    num_headerlines = rtext - rdata + 1;
    tstampstart = num_headerlines + 1; %+1 for line of '--' after header
    Timestamp = S.textdata(tstampstart:end);
    Temperature = S.data;
    Position = Temperature(1,:);
    Temperature(1,:) = [];
    
    %Extract file name from text file
    Split_Names = strsplit(S.textdata{1},'\');
    Split_Ext = strsplit(Split_Names{end},'.');
    File_Name = Split_Ext{1};
    set(handles.Selected_File,'String',File_Name);
    
    %Extract time portion of Timestamp and define a relative time (start t=0)
    Time_sec = zeros(length(Timestamp),1); Relative_Time = Time_sec;
    for i=1:length(Timestamp)
        Raw_Time = char(Timestamp{i});
        Day_Time = Raw_Time(12:26);
        Hour = str2double(Day_Time(1:2));
        Minute = str2double(Day_Time(4:5));
        Second = str2double(Day_Time(7:15));
        Time_sec(i) = Hour*60^2+Minute*60+Second;
        Relative_Time(i) = Time_sec(i)-Time_sec(1);
    end
    
    %The text file from ODiSI-B repeats the last line of data twice. This
    %checks to see if that has occurred and removes the line if it exists.
    if Relative_Time(end)==Relative_Time(end-1)
        Temperature(end,:)=[];
        Timestamp(end)=[];
        Time_sec(end)=[];
        Relative_Time(end)=[];
    end
    
    %Set rough initial limits based on input file locations
    dx = Position(1,2)-Position(1,1);
    Check_dx = zeros(length(Position)-1);
    for i=1:length(Position)-1
        Check_dx(i) = Position(1,i+1) - Position(1,i);
    end
    Break_Point = find(Check_dx > 1.5*dx);
    
    %If Break_Point is an empty matrix this indicates that the source file
    %has not yet been spatially modified. This effects the 'midpoint'.
    lp1_begin = num2str(Position(1,1));
    lp2_end = num2str(Position(1,end));
    set(handles.Loop_1_Pos_Start,'String',lp1_begin);
    set(handles.Loop_2_Pos_Stop,'String',lp2_end);
    if isempty(Break_Point) == 1
        %Break_Point IS an empty matrix (it is NOT a spatially modified
        %file). Define rough middle point.
        lp1_end = num2str(Position(1,end)/2);
        lp2_begin = num2str(Position(1,end)/2);
        set(handles.Loop_1_Pos_Stop,'String',lp1_end);
        set(handles.Loop_2_Pos_Start,'String',lp2_begin);
    else
        %Break_Point is NOT an empty matrix (it IS a spatially modified
        %file). Define exact midpoint limits.
        lp1_end = num2str(Position(1,Break_Point(1)));
        lp2_begin = num2str(Position(1,Break_Point(1)+1));
        set(handles.Loop_1_Pos_Stop,'String',lp1_end);
        set(handles.Loop_2_Pos_Start,'String',lp2_begin);
    end
%     set(handles.Loop_1_Pos_Start,'Enable','on','BackgroundColor', [0.5253 1.0 0.6670]);
% set(handles.Loop_1_Pos_Stop,'Enable','on','BackgroundColor', [1.0 0.4881 0.5098]);
% set(handles.Loop_2_Pos_Start,'Enable','on','BackgroundColor', [0.5253 1.0 0.6670]);
% set(handles.Loop_2_Pos_Stop,'Enable','on','BackgroundColor', [1.0 0.4881 0.5098])

    %Find and set limits on time controls
    set(handles.Time_Start,'String',num2str(Relative_Time(1)));
    set(handles.Time_Stop,'String',num2str(Relative_Time(end))); 
    
    %Set temperature limits for reference
    set(handles.Min_Temp,'String',num2str(min(Temperature(:))));
    set(handles.Max_Temp,'String',num2str(max(Temperature(:))));
        
    %Creates a modified S.textdata of the same form less the extra line at
    %the end if it existed
    mod_Stextdata = [S.textdata(1:num_headerlines);Timestamp];
    
    %Create handles to pass variables to all other functions
    handles.Temperature = Temperature;
    handles.Timestamp = Timestamp;
    handles.Relative_Time = Relative_Time;
    handles.Position = Position;
    handles.mod_Stextdata = mod_Stextdata;
    handles.num_headerlines = num_headerlines;
    handles.tstampstart = tstampstart;
    handles.dx = dx;
    guidata(hObject,handles);
    
catch
    %Error warning
    prompt = 'Error has occurred. Check that .txt file has the correct format';
    msgbox(prompt,'Error','error');
    set(handles.Status_Indicator,'String','Idle');
    EZ_Temperature_Viewer('Stop_Button_Callback',hObject,[],handles);
    return;
end

EZ_Temperature_Viewer('Stop_Button_Callback',hObject,[],handles);
%Set initial advice/warning
line1 = 'Check/change temp, pos, and time limits to display desired data.';
line2 = 'Be advised, saving large modified text files/creating long movies may take a while.';
line3 = 'The only way to abort "Saving" the modified text file is to forcefully close the program. Sorry :(';
set(handles.Warning_Msg,'String',{line1,line2,line3});

%=========================================================================%

% --- Executes on button press in Run_Prog_Button.
function Run_Prog_Button_Callback(hObject, eventdata, handles)

set(handles.Status_Indicator,'String','Getting Data');
pause(0.01);
set(handles.Warning_Msg,'String',{});

if isfield(handles,'Temperature')==0
    %USER has failed to load a .txt file. Notify USER and stop run.
    msgbox({'No file specified.','Please load a .txt file and try again.'});
    set(handles.Status_Indicator,'String','Idle');
    return;
end

%Define start/stop variables from USER input
lp1_start = str2double(get(handles.Loop_1_Pos_Start,'string'));
lp1_stop  = str2double(get(handles.Loop_1_Pos_Stop,'string'));
lp2_start = str2double(get(handles.Loop_2_Pos_Start,'string'));
lp2_stop  = str2double(get(handles.Loop_2_Pos_Stop,'string'));
t_start = str2double(get(handles.Time_Start,'string'));
t_stop = str2double(get(handles.Time_Stop,'string'));

%Locate index associated with USER input
[~, lp1_index_start] = min(abs(handles.Position-lp1_start));
[~, lp1_index_stop] = min(abs(handles.Position-lp1_stop));
[~, lp2_index_start] = min(abs(handles.Position-lp2_start));
[~, lp2_index_stop] = min(abs(handles.Position-lp2_stop));
[~, t_index_start] = min(abs(handles.Relative_Time-t_start));
[~, t_index_stop] = min(abs(handles.Relative_Time-t_stop));

%Corrects for an extra index outside of data (causes a horizontal line
%to appear on plots).
%Loop 1 Start
if (handles.Position(1,lp1_index_start+1)-handles.Position(1,lp1_index_start)) > handles.dx*1.5
    lp1_index_start = lp1_index_start + 1;
end
%Loop 1 Stop
if (handles.Position(1,lp1_index_stop)-handles.Position(1,lp1_index_stop-1)) > handles.dx*1.5
    lp1_index_stop = lp1_index_stop - 1;
end
%Loop 2 Start
if (handles.Position(1,lp2_index_start+1)-handles.Position(1,lp2_index_start)) > handles.dx*1.5
    lp2_index_start = lp2_index_start + 1;
end
%Loop 2 Stop
if (handles.Position(1,lp2_index_stop)-handles.Position(1,lp2_index_stop-1)) > handles.dx*1.5
    lp2_index_stop = lp2_index_stop - 1;
end

%Makes plot limit input cells inactive until program stops
set(handles.y_min_input,'Enable','off');
set(handles.y_max_input,'Enable','off');
set(handles.Loop_1_Pos_Start,'Enable','off');
set(handles.Loop_1_Pos_Stop,'Enable','off');
set(handles.Loop_2_Pos_Start,'Enable','off');
set(handles.Loop_2_Pos_Stop,'Enable','off');
set(handles.Time_Start,'Enable','off');
set(handles.Time_Stop,'Enable','off');
set(handles.Run_Prog_Button,'Enable','inactive','BackgroundColor',[0.8 0.8 0.8]);
set(handles.Filter_Toggle,'Enable','off');
set(handles.Plot_Toggle,'Enable','off');
set(handles.Mod_txt_Toggle,'Enable','off');
set(handles.Save_Vid_Toggle,'Enable','off');
set(handles.Temp_Popmenu,'Enable','off');
set(handles.Temp_bulk_browse,'Enable','off');
set(handles.Constant_Bulk_Ref,'Enable','off');
set(handles.Browse_Button,'Enable','off');
set(handles.Waterfall_Toggle,'Enable','off');
set(handles.Freq_Wplot,'Enable','off');

%Activates stop and pause controls
set(handles.Stop_Button,'Enable','on','BackgroundColor', 'red');
set(handles.Pause_Button,'Enable','on','BackgroundColor', 'yellow');

%><><><><><><><><><><> Base Temp Correction Section <><><><><><><><><><><>%
set(handles.Status_Indicator,'String','Scaling Data');
%Section for choice of base temperature
Temp_Choice = get(handles.Temp_Popmenu,'Value');

switch Temp_Choice
    case 1
        %Delta temperature. This is the form of the raw data.
    case 2
        %Increases raw data by the bulk reference temperature
        handles.Temperature = handles.Temperature + str2double(get(handles.Constant_Bulk_Ref,'String'));
    case 3
        try
            
            %Interpolate bulk temperature file to the correct time
            mod_Bulk_Temp = interp1(handles.Bulk_Time,handles.Bulk_Temp,handles.Relative_Time);
            
            %Check end time of files for correction of all NaN from interp1
            if handles.Relative_Time(end)>handles.Bulk_Time(end)
                line1 = 'Length of raw data file exceeds the time of the bulk temperature file.';
                line2 = sprintf('All temperatures after t = %d sec take the bulk temperature after that time as a constant.',handles.Bulk_Time(end));
                set(handles.Warning_Msg,'String',{line1,line2});
                %Replace all NaN with the last Bulk_Temp measurement
                mod_Bulk_Temp(isnan(mod_Bulk_Temp)) = handles.Bulk_Time(end);
            end
            
            %Increases raw data by time dependent interpolated bulk temperature
            [rowsTemp, ~] = size(handles.Temperature);
            for i=1:rowsTemp
                handles.Temperature(i,:) = handles.Temperature(i,:) + mod_Bulk_Temp(i);
            end
            
        catch
            %Upon error USER is prompted with following message.
            prompt1 = 'Bulk Temperature text file is of incorrect format';
            prompt2 = 'Please load a tab delimited Nx2 .txt file,';
            prompt3 = 'where N can be any length with the first column being time (starting at t=0) and the second column is temperature (°C)';
            uiwait(msgbox({prompt1,prompt2,prompt3}));
            
            %Restore defaults
            EZ_Temperature_Viewer('Stop_Button_Callback',hObject,[],handles);
            return;
        end
end
    
%^^^^^^^^^^^^^^^^^^^End Base Temp Correction Section^^^^^^^^^^^^^^^^^^^^^^%


%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%
%<><><><><><><><><><><><> Filter Data Section <><><><><><><><><><><><><><>%
if get(handles.Filter_Toggle,'Value') == 1
    set(handles.Status_Indicator,'String','Filtering Data');
    pause(0.01);
%Add in Filter Calculations here
%
%Temperature data with the USER specified adjustment to the base
%temperature can be accessed here with handles.Temperature (MxN)
%Time is accessed with: handles.Relative_Time (Mx1)
%Position is accessed with: handles.Position (1xN)
%
%Example:
%[handles.Temperature] = filterfunction(handles.Temperature,handles.Relative_Time,handles.Position)
%
%Note: the returned handles.Temperature should still be an MxN matrix or
%other variables could be created such as frequency, etc. and plotted either
%here or via the calling function and not modifying handles.Temperature
%
%Example:
%[FFT_data] = filterfunction(handles.Temperature,handles.Relative_Time,handles.Position
%
%Alternative: Instead of a USER function, a new GUI could be added such as
%a button that pops up a new GUI where many different filter options could
%be programmed and then accessed here when the program is executed
%
%!!!Important!!!
%If you want to modify the handles.Temperature and save a modified .txt
%file with notes within the text file about what filter settings were used,
%this can be accomplished by inserting additional headerlines in
%handles.mod_Stextdata. If this is done you must NOT change the first
%header line as that is explicitly used for file name and location
%purposes. It is permissible to add any number of lines of string (1
%column wide with as many characters as Matlab permits) after the first
%headerline and before the line with the dashes ('----').
%
%If this is done, additional options would be to include a new code in
%Browse_Button_Callback function that would identify filters used and
%create some sort of report so the user can load a file, plot the data, and
%know what filters were used to obtain the file, what the source file was,
%etc.
end
%^^^^^^^^^^^^^^^^^^^^^^^ End Filter Data Section ^^^^^^^^^^^^^^^^^^^^^^^^^%
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%

%<><><><><><><><><><><> Temperature Event Section <><><><><><><><><><><><>%
if get(handles.Event_Toggle,'Value') == 1
    set(handles.Status_Indicator,'String','Creating Event Log');
    pause(0.01);
    
    Event_Temp = str2double(get(handles.Event_Trip,'String'));
    Ignore_Temp = str2double(get(handles.Ignore_Temperature,'String'));
    
    %Loop 1 Event Log
    %Index location of all temperatures meeting criteria of the set trip
    lp1_Temps = handles.Temperature(t_index_start:t_index_stop,lp1_index_start:lp1_index_stop);
    [lp1_row_trip,lp1_col_trip] = find(lp1_Temps >= Event_Temp);
    
    %If triggering events for Loop 1 were detected they will be displayed
    if isempty(lp1_row_trip) == 0
        
        %Remove reference events >= the ignored temperature limit
        for i=1:length(lp1_row_trip)
            if lp1_Temps(lp1_row_trip(i),lp1_col_trip(i)) >= Ignore_Temp
                lp1_row_trip(i) = NaN;
                lp1_col_trip(i) = NaN;
            end
        end
        [I,~] = find(isnan(lp1_row_trip));
        lp1_row_trip(I,:) = [];
        lp1_col_trip(I,:) = [];
        
        %Recheck if there are still events satisfying criteria and proceed
        if isempty(lp1_row_trip) == 0
            
            %Remove reference events >= the ignored temperature limit
            for i=1:length(lp1_row_trip)
                if lp1_Temps(lp1_row_trip(i),lp1_col_trip(i)) >= Ignore_Temp
                    lp1_row_trip(i) = [];
                    lp1_col_trip(i) = [];
                end
            end
            
            %Sort rows in array in ascending order according to time of event
            lp1_rowcol = [lp1_row_trip,lp1_col_trip];
            lp1_rowcol_sorted = sortrows(lp1_rowcol,1);
            
            %Determine the number of unique times (rows) Event was tripped
            lp1_unique_rows = unique(lp1_rowcol_sorted(:,1));
            lp1_num_unique_rows = length(lp1_unique_rows);
            lp1_hist = histc(lp1_row_trip,lp1_unique_rows);
                        
            %Pre-allocate arrays with NaN
            lp1_pos_trip = nan(lp1_num_unique_rows,max(lp1_hist));
            lp1_time_trip = zeros(lp1_num_unique_rows,1);
            
            for i=1:lp1_num_unique_rows
                %Create array for row names based on unique event number
                lp1_row_names{i,1} = sprintf('Event #%d',i);
            end
            
            k = 0;
            lp1_col_names{1,1} = 'Time [sec]';
            for i=2:max(lp1_hist)+1
                %Create array for column names based on unique event number (k)
                k = k + 1;
                lp1_col_names{i,1} = sprintf('Pos #%d [m]',k);
            end
            
            for i=1:lp1_num_unique_rows
                %Create array of relavent time locations of events
                lp1_time_trip(i,1) = handles.Relative_Time(t_index_start + lp1_unique_rows(i) - 1);
            end
            
            k = 0;
            
            for i=1:lp1_num_unique_rows
                for j=1:lp1_hist(i)
                    %Loop over the unique times and place all positions
                    %associated that satisfy the temperature event at that time
                    k = k + 1;
                    lp1_pos_trip(i,j) = handles.Position(lp1_index_start + lp1_rowcol_sorted(k,2) - 1);
                end
            end
            
            %Condense result to one array and update table
            lp1_table_data = [lp1_time_trip,lp1_pos_trip];
            set(handles.Loop_1_Event_Log,'RowName',lp1_row_names,...
                'ColumnName',lp1_col_names,...
                'Data',lp1_table_data);
        else
            %Restores the default display
            lp1_row_names = {'Event #1';'Event #2'};
            lp1_col_names = {'Time [sec]','Pos #1 [m]'};
            lp1_table_data = [];
            set(handles.Loop_1_Event_Log,'RowName',lp1_row_names,...
                'ColumnName',lp1_col_names,...
                'Data',lp1_table_data);
        end
    else
        %Restores the default display
        lp1_row_names = {'Event #1';'Event #2'};
        lp1_col_names = {'Time [sec]','Pos #1 [m]'};
        lp1_table_data = [];
        set(handles.Loop_1_Event_Log,'RowName',lp1_row_names,...
            'ColumnName',lp1_col_names,...
            'Data',lp1_table_data);
    end
    
    %Loop 2 Event Log
    %Index location of all temperatures meeting criteria of the set trip
    lp2_Temps = handles.Temperature(t_index_start:t_index_stop,lp2_index_start:lp2_index_stop);
    [lp2_row_trip,lp2_col_trip] = find(lp2_Temps >= Event_Temp);
    
    %If triggering events for Loop 2 were detected they will be displayed
    if isempty(lp2_row_trip) == 0
        
        %Remove reference events >= the ignored temperature limit
        for i=1:length(lp2_row_trip)
            if lp2_Temps(lp2_row_trip(i),lp2_col_trip(i)) >= Ignore_Temp
                lp2_row_trip(i) = NaN;
                lp2_col_trip(i) = NaN;
            end
        end    
        [I,~] = find(isnan(lp2_row_trip));
        lp2_row_trip(I,:) = [];
        lp2_col_trip(I,:) = [];
            
        %Recheck if there are still events satisfying criteria and proceed
        if isempty(lp2_row_trip) == 0
                       
            %Sort rows in array in ascending order according to time of event
            lp2_rowcol = [lp2_row_trip,lp2_col_trip];
            lp2_rowcol_sorted = sortrows(lp2_rowcol,1);
            
            %Determine the number of unique times (rows) Event was tripped
            lp2_unique_rows = unique(lp2_rowcol_sorted(:,1));
            lp2_num_unique_rows = length(lp2_unique_rows);
            lp2_hist = histc(lp2_row_trip,lp2_unique_rows);
            
            %Pre-allocate arrays with NaN
            lp2_pos_trip = nan(lp2_num_unique_rows,max(lp2_hist));
            lp2_time_trip = zeros(lp2_num_unique_rows,1);
            
            for i=1:lp2_num_unique_rows
                %Create array for row names based on unique event number
                lp2_row_names{i,1} = sprintf('Event #%d',i);
            end
            
            k = 0;
            lp2_col_names{1,1} = 'Time [sec]';
            for i=2:max(lp2_hist)+1
                %Create array for column names based on unique event number (k)
                k = k + 1;
                lp2_col_names{i,1} = sprintf('Pos #%d [m]',k);
            end
            
            for i=1:lp2_num_unique_rows
                %Create array of relavent time locations of events
                lp2_time_trip(i,1) = handles.Relative_Time(t_index_start + lp2_unique_rows(i) - 1);
            end
            
            k = 0;
            
            for i=1:lp2_num_unique_rows
                for j=1:lp2_hist(i)
                    %Loop over the unique times and place all positions
                    %associated that satisfy the temperature event at that time
                    k = k + 1;
                    lp2_pos_trip(i,j) = handles.Position(lp2_index_start + lp2_rowcol_sorted(k,2) - 1);
                end
            end
            
            %Condense result to one array and update table
            lp2_table_data = [lp2_time_trip,lp2_pos_trip];
            set(handles.Loop_2_Event_Log,'RowName',lp2_row_names,...
                'ColumnName',lp2_col_names,...
                'Data',lp2_table_data);
        else
            %Restores the default display
            lp2_row_names = {'Event #1';'Event #2'};
            lp2_col_names = {'Time [sec]','Pos #1 [m]'};
            lp2_table_data = [];
            set(handles.Loop_2_Event_Log,'RowName',lp2_row_names,...
                'ColumnName',lp2_col_names,...
                'Data',lp2_table_data);
        end
    else
        %Restores the default display
        lp2_row_names = {'Event #1';'Event #2'};
        lp2_col_names = {'Time [sec]','Pos #1 [m]'};
        lp2_table_data = [];
        set(handles.Loop_2_Event_Log,'RowName',lp2_row_names,...
            'ColumnName',lp2_col_names,...
            'Data',lp2_table_data);
    end
end
%^^^^^^^^^^^^^^^^^^^^ End Temperature Event Section ^^^^^^^^^^^^^^^^^^^^^^%

%Update temperature limits for reference based on corrected values
set(handles.Min_Temp,'String',num2str(min(handles.Temperature(:))));
set(handles.Max_Temp,'String',num2str(max(handles.Temperature(:))));

%<><><><><><><><><><> Save Waterfall Plot Section <><><><><><><><><><><><>%

if get(handles.Waterfall_Toggle,'Value') == 1
    set(handles.Status_Indicator,'String','Creating Waterfalls');
    
    try
        %Establish the max and min temperatures from USER input for axis limits
        tempmin = str2double(get(handles.y_min_input,'String'));
        tempmax = str2double(get(handles.y_max_input,'String'));
        
        %Get USER defined frame rate
        freq_w = str2double(get(handles.Freq_Wplot,'String'));
        
        %Number of total frames to be plotted
        num_frames = floor(length(t_index_start:t_index_stop)/freq_w);
        
        %Create new arrays based on USER input
        lp1_wpos = handles.Position(lp1_index_start:lp1_index_stop)';
        lp2_wpos = handles.Position(lp2_index_start:lp2_index_stop)';
        
        lp1_wtemp = zeros(num_frames,length(lp1_index_start:lp1_index_stop));
        lp2_wtemp = zeros(num_frames,length(lp2_index_start:lp2_index_stop));
        lp_wtime = zeros(num_frames,1);
        for i=1:num_frames
            lp1_wtemp(i,:) = handles.Temperature(t_index_start+(i-1)*freq_w,lp1_index_start:lp1_index_stop);
            lp2_wtemp(i,:) = handles.Temperature(t_index_start+(i-1)*freq_w,lp2_index_start:lp2_index_stop);
            lp_wtime(i,1) = handles.Relative_Time(t_index_start+(i-1)*freq_w);
        end
        
        %For aesthetics, this limits the plot data to the limits of the plot
        lp1_wtemp(lp1_wtemp>tempmax)=tempmax;
        lp1_wtemp(lp1_wtemp<tempmin)=tempmin;
        lp2_wtemp(lp2_wtemp>tempmax)=tempmax;
        lp2_wtemp(lp2_wtemp<tempmin)=tempmin;
        
        %Close exisiting Eplot
        close(findobj('type','figure','Name', 'Waterfall Plot'));
        
        %Open new figure for Waterfall plots
        lp1_wplot = figure('Name','Waterfall Plot','NumberTitle','off','Visible','off');
        set(lp1_wplot,'Position',[0,0,650,250]);
        movegui(lp1_wplot,'northwest');
        set(lp1_wplot,'Visible','on');
        
        %Loop 1 Waterfall plot
        subplot(1,2,1)
        waterfall(lp1_wpos,lp_wtime,lp1_wtemp);
        axis([lp1_start lp1_stop t_start t_stop tempmin tempmax]);
        l1_c = colorbar;
        ylabel(l1_c,'Temperature [°C]');
        xlabel('Position [m]');
        ylabel('Time [sec]');
        zlabel('Temperature [°C]');
        title('Loop 1 Waterfall Plot');
        
        %Loop 2 Waterfall plot
        subplot(1,2,2)
        waterfall(lp2_wpos,lp_wtime,lp2_wtemp);
        axis([lp2_start lp2_stop t_start t_stop tempmin tempmax]);
        l2_c = colorbar;
        ylabel(l2_c,'Temperature [°C]');
        xlabel('Position [m]');
        ylabel('Time [sec]');
        zlabel('Temperature [°C]');
        title('Loop 2 Waterfall Plot');
    catch
        msgbox('Check that frame frequency is less than or equal to the number of frames avaible in the specified time span','Error','error');
        EZ_Temperature_Viewer('Stop_Button_Callback',hObject,[],handles);
        close(lp1_wplot);
        return;
    end
end
%^^^^^^^^^^^^^^^^^^^^^ End Waterfall Plot Section ^^^^^^^^^^^^^^^^^^^^^^^^%

%><><><><><><><><><><> Save Text File Section <><><><><><><><><><><><><><>%
%Routine to Check for and create modified .txt file from user inputs
if get(handles.Mod_txt_Toggle,'Value') == 1
    set(handles.Status_Indicator,'String','Saving Data');
    
    %Turn off button features as they will not work while text file is
    %created
    set(handles.Stop_Button,'Enable','inactive','BackgroundColor', [0.8 0.8 0.8]);
    set(handles.Pause_Button,'Enable','inactive','BackgroundColor', [0.8 0.8 0.8]);
        
    %Open SaveAs Dialog box. Will prompt if file already exists.
    default_txt = strcat(get(handles.Selected_File,'String'), ' modified');
    [file, path] = uiputfile('*.txt','Save Text File As',default_txt);
    if path == 0;
        %USER selected cancel. End callback.
        EZ_Temperature_Viewer('Stop_Button_Callback',hObject,[],handles);
        return;
    end
    
    %Check that source file headerline is in the correct format.
    Split_Names = strsplit(handles.mod_Stextdata{1},'\t');
    if length(Split_Names) ~= 2
        %Error due to disallowed modification to source file. Reset
        %defaults and end run.
        line1 = 'Filename from source .txt file is either missing or has been modified to an unsupported format.';
        line2 = 'Please load a new .txt file or modify the text file so the first entry is of the format:';
        line3 = sprintf('Filename \tC:\\insert_full_file_path_and_name_here');
        line4 = 'Note: The previous line is tab delimited between Filename and C:\...';
        uiwait(msgbox({line1,line2,line3,line4}));
        
        EZ_Temperature_Viewer('Stop_Button_Callback',hObject,[],handles);
        return;
    end
    
    %If file exists in searched path, the file will be deleted so a new
    %file can be written. This required USER permission.
    Txt_Label = fullfile(path,file);
    if exist(Txt_Label,'file') == 2
        delete(Txt_Label);
    end
    
    %Replaces old file path/name with the new USER defined path/name
    Split_Names{2} = Txt_Label;
    mod_File_Name = strcat(Split_Names{1},{char(9)},Split_Names{2});
    handles.mod_Stextdata{1}= mod_File_Name;
    
    %Create and write headerlines to .txt file
    headerlines = handles.mod_Stextdata(1:handles.num_headerlines-1);
    dlmcell(Txt_Label,headerlines);

    %Create modified arrays based on USER input for export
    mod_Position = [handles.Position(lp1_index_start:lp1_index_stop) handles.Position(lp2_index_start:lp2_index_stop)];
    mod_Temperature = [handles.Temperature(t_index_start:t_index_stop,lp1_index_start:lp1_index_stop) handles.Temperature(t_index_start:t_index_stop,lp2_index_start:lp2_index_stop)];
    mod_Sdata = [mod_Position; mod_Temperature];
    mod_t_index_start = handles.num_headerlines + t_index_start;
    mod_t_index_stop = handles.num_headerlines + t_index_stop;
    mod_Timestamp = [handles.mod_Stextdata(handles.num_headerlines);handles.mod_Stextdata(mod_t_index_start:mod_t_index_stop)];
    
    %Convert Timestamp and data to one cell array and append to .txt file
    cell_nums=num2cell(mod_Sdata);
    cell_array=[mod_Timestamp cell_nums];
    dlmcell(Txt_Label,cell_array,'-a');
    
    %Restore button features when the text file is completed
    set(handles.Stop_Button,'Enable','on','BackgroundColor', 'red');
    set(handles.Pause_Button,'Enable','on','BackgroundColor', 'yellow');
end
%^^^^^^^^^^^^^^^^^^^^ End Save Text File Section ^^^^^^^^^^^^^^^^^^^^^^^^^%


%><><><><><><><><><><> Save Movie .avi Section <><><><><><><><><><><><><><%
%Routine to check if USER wants to save the plot to a movie
if get(handles.Save_Vid_Toggle,'Value')==1
    set(handles.Status_Indicator,'String','Making Movie');
    %Open SaveAs Dialog box. Will prompt if file already exists.
    default_vid = strcat(get(handles.Selected_File,'String'));
    [file, path] = uiputfile('*.avi','Save Video As',default_vid);
    
    if path == 0;
        %USER selected cancel. Reset defaults and End callback.
        EZ_Temperature_Viewer('Stop_Button_Callback',hObject,[],handles);
        return;
    end
    
    %If file exists in searched path, the file will be deleted so a new
    %file can be written
    Vid_Label = fullfile(path,file);
    if exist(Vid_Label,'file') == 2
        delete(Vid_Label);
    end
    
    %Create JPEG compressed .avi file for later playback of data
    vidObj = VideoWriter(Vid_Label);
    open(vidObj);
    
    %Plot and frame capture for movie routine
    for i = t_index_start:t_index_stop
        
        %Controls response to toggle of the Stop_Button
        if get(handles.Stop_Button,'Value') == 1
            %Restores default stop value and exits run
            close(vidObj);
            EZ_Temperature_Viewer('Stop_Button_Callback',hObject,[],handles)
            set(handles.Stop_Button,'Value',0);
            return;
        end
        
        %Controls response to toggle of the Pause_Button
        while get(handles.Pause_Button,'Value') == 1
            pause(0.01);
        end
        
        %Displays relative time associated with potted data for reference
        set(handles.Timer,'String',num2str(handles.Relative_Time(i)));
        set(handles.Date,'String',handles.Timestamp(i));
        
        %Plot data for Loop 1
        set(handles.Loop_1_Temp_Plot,'NextPlot','replacechildren');
        plot(handles.Loop_1_Temp_Plot,handles.Position(lp1_index_start:lp1_index_stop),handles.Temperature(i,lp1_index_start:lp1_index_stop));
        
        %Plot data for Loop 2
        set(handles.Loop_2_Temp_Plot,'NextPlot','replacechildren');
        plot(handles.Loop_2_Temp_Plot,handles.Position(lp2_index_start:lp2_index_stop),handles.Temperature(i,lp2_index_start:lp2_index_stop));

        %Show plots
        drawnow;
        
        %Add a frame to the .avi
        writeVideo(vidObj,getframe(handles.Main_Figure, [0 0 1000 325]));
    end
    
    %Close .avi file
    close(vidObj);
end
%^^^^^^^^^^^^^^^^^^^^ End Save Movie .avi Section ^^^^^^^^^^^^^^^^^^^^^^^^%


%><><><><><><><><><><><><> Plot Data Section <><><><><><><><><><><><><><><%
%Check if USER wants to plot the data
if get(handles.Plot_Toggle,'Value') == 1 && get(handles.Save_Vid_Toggle,'Value') == 0
    set(handles.Status_Indicator,'String','Plotting Data');
    
    %Plot routine
    for i = t_index_start:t_index_stop
        
        %Controls response to toggle of the Stop_Button
        if get(handles.Stop_Button,'Value') == 1
            %Restores default stop value and exits run
            EZ_Temperature_Viewer('Stop_Button_Callback',hObject,[],handles)
            set(handles.Stop_Button,'Value',0);
            return;
        end
        
        %Controls response to toggle of the Pause_Button
        while get(handles.Pause_Button,'Value') == 1
            pause(0.01);
        end
        
        %Displays relative time associated with potted data for reference
        set(handles.Timer,'String',num2str(handles.Relative_Time(i)));
        set(handles.Date,'String',handles.Timestamp(i));
        
        %Plot option of the event trip line showing on plots
        if (get(handles.Event_Line_Toggle,'Value') == 1 && get(handles.Event_Toggle,'Value') == 1)
            Event_Line_1(1:length(handles.Position(lp1_index_start:lp1_index_stop)),1) = str2num(get(handles.Event_Trip,'String'));
            Event_Line_2(1:length(handles.Position(lp2_index_start:lp2_index_stop)),1) = str2num(get(handles.Event_Trip,'String'));
            
            %Plot data for Loop 1
            set(handles.Loop_1_Temp_Plot,'NextPlot','replacechildren');
            plot(handles.Loop_1_Temp_Plot,handles.Position(lp1_index_start:lp1_index_stop),handles.Temperature(i,lp1_index_start:lp1_index_stop),...
                                          handles.Position(lp1_index_start:lp1_index_stop),Event_Line_1,'r-');
            
            %Plot data for Loop 2
            set(handles.Loop_2_Temp_Plot,'NextPlot','replacechildren');
            plot(handles.Loop_2_Temp_Plot,handles.Position(lp2_index_start:lp2_index_stop),handles.Temperature(i,lp2_index_start:lp2_index_stop),...
                                          handles.Position(lp2_index_start:lp2_index_stop),Event_Line_2,'r-');
            
            %Show plots
            drawnow;
        else
            %Plot data for Loop 1
            set(handles.Loop_1_Temp_Plot,'NextPlot','replacechildren');
            plot(handles.Loop_1_Temp_Plot,handles.Position(lp1_index_start:lp1_index_stop),handles.Temperature(i,lp1_index_start:lp1_index_stop));
            
            %Plot data for Loop 2
            set(handles.Loop_2_Temp_Plot,'NextPlot','replacechildren');
            plot(handles.Loop_2_Temp_Plot,handles.Position(lp2_index_start:lp2_index_stop),handles.Temperature(i,lp2_index_start:lp2_index_stop));
            
            %Show plots
            drawnow;
        end
    end
end
%^^^^^^^^^^^^^^^^^^^^^^^^ End Plot Data Section ^^^^^^^^^^^^^^^^^^^^^^^^^^%

%Restore default conditions from before program was run
EZ_Temperature_Viewer('Stop_Button_Callback',hObject,[],handles);

%=========================================================================%

function Loop_1_Pos_Start_Callback(hObject, eventdata, handles)
%Check to see that min/max plot constraints are fulfilled otherwise default
%values are restored
default = '0';
[~,status] = str2num(get(hObject,'String'));
if status == 1
    ymin = str2double(get(handles.y_min_input,'String'));
    ymax = str2double(get(handles.y_max_input,'String'));
    lp1_start = str2double(get(handles.Loop_1_Pos_Start,'string'));
    lp1_stop  = str2double(get(handles.Loop_1_Pos_Stop,'string'));
    lp2_start = str2double(get(handles.Loop_2_Pos_Start,'string'));
    lp2_stop  = str2double(get(handles.Loop_2_Pos_Stop,'string'));
    
    if lp1_start >= lp1_stop
        set(hObject,'String',default);
    else
        axes(handles.Loop_1_Temp_Plot);
        axis([lp1_start lp1_stop ymin ymax]);
        
        axes(handles.Loop_2_Temp_Plot);
        axis([lp2_start lp2_stop ymin ymax]);
    end
else
    set(hObject,'String',default);
end;

%=========================================================================%

function Loop_1_Pos_Start_CreateFcn(hObject, eventdata, handles)

%=========================================================================%

function Loop_2_Start_Pos_Callback(hObject, eventdata, handles)
%Check to see that min/max plot constraints are fulfilled otherwise default
%values are restored
default = '2.5';
[~,status] = str2num(get(hObject,'String'));
if status == 1
    ymin = str2double(get(handles.y_min_input,'String'));
    ymax = str2double(get(handles.y_max_input,'String'));
    lp1_start = str2double(get(handles.Loop_1_Pos_Start,'string'));
    lp1_stop  = str2double(get(handles.Loop_1_Pos_Stop,'string'));
    lp2_start = str2double(get(handles.Loop_2_Pos_Start,'string'));
    lp2_stop  = str2double(get(handles.Loop_2_Pos_Stop,'string'));
    
    if lp2_start >= lp2_stop
        set(hObject,'String',default);
    else
        axes(handles.Loop_1_Temp_Plot);
        axis([lp1_start lp1_stop ymin ymax]);
        
        axes(handles.Loop_2_Temp_Plot);
        axis([lp2_start lp2_stop ymin ymax]);
    end
else
    set(hObject,'String',default);
end;

%=========================================================================%

function Loop_2_Start_Pos_CreateFcn(hObject, eventdata, handles)

%=========================================================================%

function Loop_1_Pos_Stop_Callback(hObject, eventdata, handles)
%Check to see that min/max plot constraints are fulfilled otherwise default
%values are restored
default = '2.5';
[~,status] = str2num(get(hObject,'String'));
if status == 1
    ymin = str2double(get(handles.y_min_input,'String'));
    ymax = str2double(get(handles.y_max_input,'String'));
    lp1_start = str2double(get(handles.Loop_1_Pos_Start,'string'));
    lp1_stop  = str2double(get(handles.Loop_1_Pos_Stop,'string'));
    lp2_start = str2double(get(handles.Loop_2_Pos_Start,'string'));
    lp2_stop  = str2double(get(handles.Loop_2_Pos_Stop,'string'));
    
    if lp1_stop <= lp1_start
        set(hObject,'String',default);
    else
        axes(handles.Loop_1_Temp_Plot);
        axis([lp1_start lp1_stop ymin ymax]);
        
        axes(handles.Loop_2_Temp_Plot);
        axis([lp2_start lp2_stop ymin ymax]);
    end
else
    set(hObject,'String',default);
end;

%=========================================================================%

function Loop_1_Pos_Stop_CreateFcn(hObject, eventdata, handles)

%=========================================================================%

function Loop_2_Pos_Stop_Callback(hObject, eventdata, handles)
%Check to see that min/max plot constraints are fulfilled otherwise default
%values are restored
default = '5';
[~,status] = str2num(get(hObject,'String'));
if status == 1
    ymin = str2double(get(handles.y_min_input,'String'));
    ymax = str2double(get(handles.y_max_input,'String'));
    lp1_start = str2double(get(handles.Loop_1_Pos_Start,'string'));
    lp1_stop  = str2double(get(handles.Loop_1_Pos_Stop,'string'));
    lp2_start = str2double(get(handles.Loop_2_Pos_Start,'string'));
    lp2_stop  = str2double(get(handles.Loop_2_Pos_Stop,'string'));
    
    if lp2_stop <= lp2_start
        set(hObject,'String',default);
    else
        axes(handles.Loop_1_Temp_Plot);
        axis([lp1_start lp1_stop ymin ymax]);
        
        axes(handles.Loop_2_Temp_Plot);
        axis([lp2_start lp2_stop ymin ymax]);
    end
else
    set(hObject,'String',default);
end;

%=========================================================================%

function Loop_2_Pos_Stop_CreateFcn(hObject, eventdata, handles)

%=========================================================================%

function y_min_input_Callback(hObject, eventdata, handles)
%Check to see that min/max plot constraints are fulfilled otherwise default
%values are restored
default = '-5';
[~,status] = str2num(get(hObject,'String'));
if status == 1
    ymin = str2double(get(handles.y_min_input,'String'));
    ymax = str2double(get(handles.y_max_input,'String'));
    lp1_start = str2double(get(handles.Loop_1_Pos_Start,'string'));
    lp1_stop  = str2double(get(handles.Loop_1_Pos_Stop,'string'));
    lp2_start = str2double(get(handles.Loop_2_Pos_Start,'string'));
    lp2_stop  = str2double(get(handles.Loop_2_Pos_Stop,'string'));
    
    if ymin >= ymax
        set(hObject,'String',default);
    else
        axes(handles.Loop_1_Temp_Plot);
        axis([lp1_start lp1_stop ymin ymax]);
        
        axes(handles.Loop_2_Temp_Plot);
        axis([lp2_start lp2_stop ymin ymax]);
    end
else
    set(hObject,'String',default);
end;

%=========================================================================%

function y_min_input_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%=========================================================================%

function y_max_input_Callback(hObject, eventdata, handles)
%Check to see that min/max plot constraints are fulfilled otherwise default
%values are restored
default = '25';
[~,status] = str2num(get(hObject,'String'));
if status == 1
    ymin = str2double(get(handles.y_min_input,'String'));
    ymax = str2double(get(handles.y_max_input,'String'));
    lp1_start = str2double(get(handles.Loop_1_Pos_Start,'string'));
    lp1_stop  = str2double(get(handles.Loop_1_Pos_Stop,'string'));
    lp2_start = str2double(get(handles.Loop_2_Pos_Start,'string'));
    lp2_stop  = str2double(get(handles.Loop_2_Pos_Stop,'string'));
    
    if ymax <= ymin
        set(hObject,'String',default);
    else
        axes(handles.Loop_1_Temp_Plot);
        axis([lp1_start lp1_stop ymin ymax]);
        
        axes(handles.Loop_2_Temp_Plot);
        axis([lp2_start lp2_stop ymin ymax]);
    end
else
    set(hObject,'String',default);
end;

%=========================================================================%

function y_max_input_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%=========================================================================%

function Timer_Callback(hObject, eventdata, handles)

%=========================================================================%

function Timer_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%=========================================================================%

function Browse_Results_Callback(hObject, eventdata, handles)

%=========================================================================%

function Browse_Results_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%=========================================================================%

function Save_Vid_Toggle_Callback(hObject, eventdata, handles)
val = get(hObject,'Value');
if val == 1
    set(handles.Plot_Toggle,'Enable','Inactive','Value',1);
else
    set(handles.Plot_Toggle,'Enable','on');
end
%=========================================================================%

function Stop_Button_Callback(hObject, eventdata, handles)
%These are the default settings to be restored when the main program is
%finished or if the stop button is activated to abort the program
set(handles.y_min_input,'Enable','on');
set(handles.y_max_input,'Enable','on');
set(handles.Loop_1_Pos_Start,'Enable','on','BackgroundColor', [0.5253 1.0 0.6670]);
set(handles.Loop_1_Pos_Stop,'Enable','on','BackgroundColor', [1.0 0.4881 0.5098]);
set(handles.Loop_2_Pos_Start,'Enable','on','BackgroundColor', [0.5253 1.0 0.6670]);
set(handles.Loop_2_Pos_Stop,'Enable','on','BackgroundColor', [1.0 0.4881 0.5098]);
set(handles.Time_Start,'Enable','on');
set(handles.Time_Stop,'Enable','on');
set(handles.Run_Prog_Button,'Enable','on','BackgroundColor',[0.0590 0.5793 0.7020]);
set(handles.Filter_Toggle,'Enable','on');
set(handles.Plot_Toggle,'Enable','on');
set(handles.Mod_txt_Toggle,'Enable','on');
set(handles.Save_Vid_Toggle,'Enable','on');
set(handles.Launch_MP_Button,'Enable','on');
set(handles.Stop_Button,'Enable','inactive','BackgroundColor', [0.8 0.8 0.8]);
set(handles.Pause_Button,'Value',0,'Enable','inactive','String','Pause','BackgroundColor', [0.8 0.8 0.8]);
set(handles.Temp_Popmenu,'Enable','on');
set(handles.Temp_bulk_browse,'Enable','on');
set(handles.Constant_Bulk_Ref,'Enable','on');
set(handles.Browse_Button,'Enable','on');
set(handles.Waterfall_Toggle,'Enable','on');
set(handles.Status_Indicator,'String','Idle');
set(handles.Freq_Wplot,'Enable','on');
fclose all;
%=========================================================================%

function Pause_Button_Callback(hObject, eventdata, handles)
if get(hObject,'Value') == 1
    set(hObject,'String','Play','BackgroundColor','green');
else
    set(hObject,'String','Pause','BackgroundColor','yellow');
end

%=========================================================================%

function Selected_File_Callback(hObject, eventdata, handles)

%=========================================================================%

function Selected_File_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%=========================================================================%

% --- Executes on button press in Launch_MP_Button.
function Launch_MP_Button_Callback(hObject, eventdata, handles)
% Check that user has the Image Processing Toolbox installed.
hasIPT = license('test', 'image_toolbox');

[file, path, ~] = uigetfile('*.avi');
if path == 0;
    %USER selected cancel. End callback.
    return;
end

if hasIPT == 1
    %Open selected file in Image Toolbox
    try
        implay(fullfile(path,file)); 
    catch
        winopen(fullfile(path,file));
    end
else
    %Open selected file in Windows Media Player
    winopen(fullfile(path,file));
end

%=========================================================================%

% --- Executes on button press in Mod_txt_Toggle.
function Mod_txt_Toggle_Callback(hObject, eventdata, handles)

%=========================================================================%

function Status_Indicator_Callback(hObject, eventdata, handles)

%=========================================================================%

function Status_Indicator_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
%=========================================================================%

function dlmcell(file,cell_array,varargin)
% <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> %
% <><><><><>     dlmcell - Write Cell Array to Text File      <><><><><> %
% <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> %
%                                                 Version:    01.06.2010 %
%                                                     (c) Roland Pfister %
%                                             roland_pfister@t-online.de %
%                        ...with many thanks to George Papazafeiropoulos %
%                        for his corrections and improvements.           %
% 1. Synopsis                                                            %
%                                                                        %
% A single cell array is written to an output file. Cells may consist of %
% any combination of (a) numbers, (b) letters, or (c) words. The inputs  %
% are as follows:                                                        %
%                                                                        %
%       - file       The output filename (string).                       %
%       - cell_array The cell array to be written.                       %
%       - delimiter  Delimiter symbol, e.g. ',' (optional;               %
%                    default: tab ('\t'}).                               %
%       - append     '-a' for appending the content to the               %
%                    output file (optional).                             %
%                                                                        %
% 2. Example                                                             %
%                                                                        %
%         mycell = {'Numbers', 'Letters', 'Words','More Words'; ...      %
%                    1, 'A', 'Apple', {'Apricot'}; ...                   %
%                    2, 'B', 'Banana', {'Blueberry'}; ...                %
%                    3, 'C', 'Cherry', {'Cranberry'}; };                 %
%         dlmcell('mytext.txt',mycell);                                  %
%                                                                        %
% <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> %

% Check input arguments
if nargin < 2
    disp('Error - Give at least two input arguments!');
    return;
elseif nargin > 4
    disp('Error - Do not give more than 4 input arguments!');
    return;
end
if ~ischar(file)
    disp(['Error - File input has to be a string (e.g. ' ...
        char(39) 'output.txt' char(39) '!']);
    return;
end;
if ~iscell(cell_array)
    disp('Error - Input cell_array not of the type "cell"!');
    return;
end;
delimiter = '\t';
append = 'w';
if nargin > 2
    for i = 1:size(varargin,2)
        if strcmp('-a',varargin{1,i}) == 1
            append = 'a';
        else
            delimiter = varargin{1,i};
        end;
    end;
end

% Open output file and prepare output array.
output_file = fopen(file,append);
output = cell(size(cell_array,1),size(cell_array,2));

% Evaluate and write input array.
for i = 1:size(cell_array,1)
    for j = 1:size(cell_array,2)
        if numel(cell_array{i,j}) == 0
            output{i,j} = '';
            % Check whether the content of cell i,j is
            % numeric and convert numbers to strings.
        elseif isnumeric(cell_array{i,j}) || islogical(cell_array{i,j})
            output{i,j} = num2str(cell_array{i,j}(1,1));
            
            % Check whether the content of cell i,j is another cell (e.g. a
            % string of length > 1 that was stored as cell. If cell sizes
            % equal [1,1], convert numbers and char-cells to strings.
            %
            % Note that any other cells-within-the-cell will produce errors
            % or wrong results.
        elseif iscell(cell_array{i,j})
            if size(cell_array{i,j},1) == 1 && size(cell_array{i,j},1) == 1
                if isnumeric(cell_array{i,j}{1,1})
                    output{i,j} = num2str(cell_array{i,j}{1,1}(1,1));
                elseif ischar(cell_array{i,j}{1,1})
                    output{i,j} = cell_array{i,j}{1,1};
                end;
            end;
            
            % If the cell already contains a string, nothing has to be done.
        elseif ischar(cell_array{i,j})
            output{i,j} = cell_array{i,j};
        end;
        
        % Cell i,j is written to the output file. A delimiter is appended for
        % all but the last element of each row. At the end of a row, a newline
        % is written to the output file.
        if j < size(cell_array,2)
            fprintf(output_file,['%s',delimiter],output{i,j});
        else
            fprintf(output_file,'%s\r\n',output{i,j});
        end
    end;
end;

% Close output file.
fclose(output_file);

%=========================================================================%

function Plot_Toggle_Callback(hObject, eventdata, handles)

%=========================================================================%

function Time_Start_Callback(hObject, eventdata, handles)
%Checks that input is convertible to number then checks for valid limits
default = '0';
[~,status] = str2num(get(hObject,'String'));
try
    if status == 1
        start_time = str2double(get(handles.Time_Start,'String'));
        stop_time = str2double(get(handles.Time_Stop,'String'));
        if start_time < handles.Relative_Time(1) || start_time >= stop_time
            set(hObject,'String',num2str(handles.Relative_Time(1)));
        end
    else
        set(hObject,'String',num2str(handles.Relative_Time(1)));
    end;
catch
    %On error (e.g. no file has yet been loaded) the default value will set
    set(hObject,'String',default);
end
%=========================================================================%
% --- Executes during object creation, after setting all properties.
function Time_Start_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%=========================================================================%

function Time_Stop_Callback(hObject, eventdata, handles)
%Checks that input is convertible to number then checks for valid limits
default = '10';
[~,status] = str2num(get(hObject,'String'));
try
    if status == 1
        start_time = str2double(get(handles.Time_Start,'String'));
        stop_time = str2double(get(handles.Time_Stop,'String'));
        if stop_time > handles.Relative_Time(end) || stop_time <= start_time
            set(hObject,'String',num2str(handles.Relative_Time(end)));
        end
    else
        set(hObject,'String',num2str(handles.Relative_Time(end)));
    end;
catch
    %On error (e.g. no file has yet been loaded) the default value will set
    set(hObject,'String',default);
end
%=========================================================================%

function Time_Stop_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%=========================================================================%

function Date_Callback(hObject, eventdata, handles)

%=========================================================================%

function Date_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%=========================================================================%

% --- Executes when user attempts to close EZ_Temperature_Viewer.
function Main_Figure_CloseRequestFcn(hObject, eventdata, handles)
fclose all;
close(findobj('type','figure','Name', 'Waterfall Plot'));
close(findobj('type','figure','Name', 'Event Plots'));
delete(hObject);

%=========================================================================%

% --- Executes on selection change in Temp_Popmenu.
function Temp_Popmenu_Callback(hObject, eventdata, handles)
set(handles.Constant_Bulk_Ref,'Visible','off');
set(handles.Temp_bulk_browse,'Visible','off')
set(handles.Bulk_Indicator,'Visible','off');
Temp_Choice = get(hObject,'Value');

switch Temp_Choice
    case 1
        %Delta temperature. This is the form of the raw data.
        axes(handles.Loop_1_Temp_Plot);
        ylabel('Delta Temperature [°C]');
        
        axes(handles.Loop_2_Temp_Plot);
        ylabel('Delta Temperature [°C]');
    case 2
        %Increases raw data by the bulk reference temperature
        set(handles.Constant_Bulk_Ref,'Visible','on');
        
        axes(handles.Loop_1_Temp_Plot);
        ylabel('Temperature [°C]');
        
        axes(handles.Loop_2_Temp_Plot);
        ylabel('Temperature [°C]');
    case 3
        %Increases raw data by a time dependent bulk temperature
        set(handles.Temp_bulk_browse,'Visible','on');
        set(handles.Bulk_Indicator,'Visible','on');
        
        axes(handles.Loop_1_Temp_Plot);
        ylabel('Temperature [°C]');
        
        axes(handles.Loop_2_Temp_Plot);
        ylabel('Temperature [°C]');
end
%=========================================================================%

function Temp_Popmenu_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%=========================================================================%

function Constant_Bulk_Ref_Callback(hObject, eventdata, handles)
%Checks that input is convertible to number then checks for valid limits
default = '25';
[~,status] = str2num(get(hObject,'String'));
try
    if status == 0
       set(hObject,'String',default);
    end
catch
    %On error (e.g. no file has yet been loaded) the default value will set
    set(hObject,'String',default);
end

%=========================================================================%

function Constant_Bulk_Ref_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%=========================================================================%

% --- Executes on button press in Temp_bulk_browse.
function Temp_bulk_browse_Callback(hObject, eventdata, handles)
B = uiimport;
if isempty(B)
    %USER clicked cancel.
    set(handles.Bulk_Indicator,'BackgroundColor','red');
    return;
end

vars = fieldnames(B);
Bulk = B.(vars{1});
if length(vars) > 1
    msgbox('Check bulk temp file for errors. It should return a Nx2 array under only one variable name');
    set(handles.Bulk_Indicator,'BackgroundColor','red');
    return;
elseif Bulk(1,1) ~= 0
    msgbox('The time dependent bulk temperature must begin at t = 0.');
    set(handles.Bulk_Indicator,'BackgroundColor','red');
    return;
end

set(handles.Bulk_Indicator,'BackgroundColor','green');
handles.Bulk_Time = Bulk(:,1);
handles.Bulk_Temp = Bulk(:,2);
guidata(hObject,handles);

%=========================================================================%

function Max_Temp_Callback(hObject, eventdata, handles)

%=========================================================================%

function Max_Temp_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%=========================================================================%

function Min_Temp_Callback(hObject, eventdata, handles)

%=========================================================================%

function Min_Temp_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%=========================================================================%

function Waterfall_Toggle_Callback(hObject, eventdata, handles)
val = get(hObject,'Value');
switch val
    case 0
    set(handles.Freq_Wplot,'Visible','off');
    set(handles.frame_freq_txt,'Visible','off');
    
    case 1
    set(handles.Freq_Wplot,'Visible','on');
    set(handles.frame_freq_txt,'Visible','on');
end
%=========================================================================%

function Freq_Wplot_Callback(hObject, eventdata, handles)
default = '100';
[val,status] = str2num(get(hObject,'String'));
try
    %Check that input can be converted to string
    if status == 0
        set(hObject,'String',default);
    else
        val = round(val);
        
        %Prevent negative frame frequency
        if val < 0
            val = 0;
        end
        
        %Update edit box
        set(hObject,'String',num2str(val));
    end
catch
    %Other error
    set(hObject,'String',default);
end

%=========================================================================%

function Freq_Wplot_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%=========================================================================%

% --- Executes on button press in Event_Toggle.
function Event_Toggle_Callback(hObject, eventdata, handles)
val = get(hObject,'Value');
if val == 1
    set(handles.Event_Trip,'Enable','on','BackgroundColor','white');
    set(handles.Ignore_Temperature,'Enable','on','BackgroundColor','white');
    set(handles.Event_Log_Plot,'Enable','on');
    set(handles.Loop_1_Log_Save,'Enable','on');
    set(handles.Loop_2_Log_Save,'Enable','on');
    set(handles.Event_Line_Toggle,'Enable','on');
elseif val == 0
    set(handles.Event_Trip,'Enable','off');
    set(handles.Ignore_Temperature,'Enable','off');
    set(handles.Event_Log_Plot,'Enable','off');
    set(handles.Loop_1_Log_Save,'Enable','off');
    set(handles.Loop_2_Log_Save,'Enable','off');
    set(handles.Event_Line_Toggle,'Enable','off');
end

%=========================================================================%

function Event_Trip_Callback(hObject, eventdata, handles)
default = '25';
[~, status] = str2num(get(hObject,'String'));
if status == 0
    %Invalid input. Value set back to default.
    set(hObject,'String',default);
end

%=========================================================================%

% --- Executes during object creation, after setting all properties.
function Event_Trip_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%=========================================================================%

function Ignore_Temperature_Callback(hObject, eventdata, handles)
default = '1000';
[val, status] = str2num(get(hObject,'String'));
if status == 0
    %Invalid input. Value set back to default.
    set(hObject,'String',default);
elseif val <= str2double(get(handles.Event_Trip,'String'));
    %Invalid input. Value set back to default.
    set(hObject,'String',default);
end

%=========================================================================%

% --- Executes during object creation, after setting all properties.
function Ignore_Temperature_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%=========================================================================%

% --- Executes on button press in Event_Log_Plot.
function Event_Log_Plot_Callback(hObject, eventdata, handles)

%Get data
Log_1_Data = get(handles.Loop_1_Event_Log,'Data');
Log_2_Data = get(handles.Loop_2_Event_Log,'Data');

%Abort plotting as user has not yet analyzed the data (will be double array
%when data has been analyzed)
if iscell(Log_1_Data);
    return;
end

%Set the time range of interest
t_start = str2double(get(handles.Time_Start,'string'));
t_stop = str2double(get(handles.Time_Stop,'string'));
[~, t_index_start] = min(abs(handles.Relative_Time-t_start));
[~, t_index_stop] = min(abs(handles.Relative_Time-t_stop));

%Initalize arrays
x_time = handles.Relative_Time(t_index_start:t_index_stop);
y1_trip = zeros(length(x_time),1);
y2_trip = zeros(length(x_time),1);

%Map time of logged events to the total time
for i=1:size(Log_1_Data,1)
    [row1,~] = find(x_time == Log_1_Data(i,1));
    y1_trip(row1,1) = 0.5;
end
for i=1:size(Log_2_Data,1)
    [row2,~] = find(x_time == Log_2_Data(i,1));
    y2_trip(row2,1) = 0.5;
end

%Close exisiting Eplot
close(findobj('type','figure','Name', 'Event Plots'));

%Initialize Event Log Plot
Eplot = figure('Name','Event Plots','NumberTitle','off','Visible','off');
set(Eplot,'Position',[0,0,600,200]);
movegui(Eplot,'northeast');
set(Eplot,'Visible','on');

%Plot time occurences of Event Log 1
subplot(1,2,1);
plot(x_time,y1_trip,'bx-')
axis([t_start t_stop 0 1]);
set(gca,'YTick',[0,0.5]);
set(gca,'YtickLabel',{'Not Tripped','Tripped'});
xlabel('Time [sec]');
title('Loop 1 Event Log');

%Plot time occurences of Event Log 2
subplot(1,2,2);
plot(x_time,y2_trip,'bx-')
axis([t_start t_stop 0 1]);
set(gca,'YTick',[0,0.5]);
set(gca,'YtickLabel',{'Not Tripped','Tripped'});
xlabel('Time [sec]');
title('Loop 2 Event Log');

%=========================================================================%

% --- Executes on button press in Loop_1_Log_Save.
function Loop_1_Log_Save_Callback(hObject, eventdata, handles)
%Prompt for file name and location
default_txt = strcat(get(handles.Selected_File,'String'), ' Loop 1 Event Log');
[file, path] = uiputfile('*.txt','Save Loop 1 Event Log File As',default_txt);
if path == 0;
    %USER selected cancel. End callback.
    return;
end

%If file exists in searched path, the file will be deleted so a new
%file can be written. This required USER permission.
Log_Label = fullfile(path,file);
if exist(Log_Label,'file') == 2
    delete(Log_Label);
end

%Get data from table
Log_Header = get(handles.Loop_1_Event_Log,'ColumnName')';
Log_Data = get(handles.Loop_1_Event_Log,'Data');

%Write data with column names to file
dlmcell(Log_Label,Log_Header)
dlmwrite(Log_Label,Log_Data,'delimiter','\t','precision','%.8f','newline','pc','-append');

%=========================================================================%

% --- Executes on button press in Loop_2_Log_Save.
function Loop_2_Log_Save_Callback(hObject, eventdata, handles)
%Prompt for file name and location
default_txt = strcat(get(handles.Selected_File,'String'), ' Loop 2 Event Log');
[file, path] = uiputfile('*.txt','Save Loop 2 Event Log File As',default_txt);
if path == 0;
    %USER selected cancel. End callback.
    return;
end

%If file exists in searched path, the file will be deleted so a new
%file can be written. This required USER permission.
Log_Label = fullfile(path,file);
if exist(Log_Label,'file') == 2
    delete(Log_Label);
end

%Get data from table
Log_Header = get(handles.Loop_2_Event_Log,'ColumnName')';
Log_Data = get(handles.Loop_2_Event_Log,'Data');

%Write data with column names to tab deliminated file
dlmcell(Log_Label,Log_Header)
dlmwrite(Log_Label,Log_Data,'delimiter','\t','precision','%.8f','newline','pc','-append');

%=========================================================================%

% --- Executes on button press in Event_Line_Toggle.
function Event_Line_Toggle_Callback(hObject, eventdata, handles)

%=========================================================================%

% --- Executes on button press in Filter_Toggle.
function Filter_Toggle_Callback(hObject, eventdata, handles)

%=========================================================================%

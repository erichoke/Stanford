% Eric Hoke 11/11/11
% This program reads all of the PL files in a folder, interpolates data to
% the desired wavelength values and then divides the spectra by a
% calibration factor that corrects for the spectrometer response. It takes
% an input of wavelength values (xvals) which the data is interpolated to.
% Subtract min (true or false or number) if true subtracts the minimum of each data series
% before doing any calibration corrections if number not equal to one will
% subtract that value
% This script requires the file SpectrometerResponse.mat in the same folder
% 10/29/13 Added normalization by integration time Eric Hoke
% 4/13/14 Modified to calculate for integrating sphere measurements
% accounting for response

function ExtractPL_Interp_Sphere(xvals,subtractmin)
xvals=xvals';
step=(max(xvals)-min(xvals)+1)/length(xvals);

% Prompts user for a file, then it looks up all files in that folder
[filename, pathname]=uigetfile('*.arc_data', 'Open any PL file in folder you want to extract from');
PLfiles=dir(strcat(pathname, '*.arc_data'));

% list of all files in folder which are not directories
PLfiles=PLfiles([PLfiles.isdir]==0);

% Loads correction factor for spectrometer response
load('SphereSpectrometerResponse','SphereSpectrometerResponse')

%initialize data holders
PLdataRaw=xvals; %raw data
RawCountsPerSec=xvals; %raw data corrected for integration time and offset value
PLdata=xvals; % corrected for mono & detector response,integration time and offset value
%PLdataNorm=xvals; %previous result normalized
%PLdataeV=[xvals 1240./xvals]; % Response corrected Converted to Energy (eV) scale accounting for Jacobian, integration time and offset value corrected
%PLdataeVNorm=[xvals 1240./xvals]; %previous result normalized

headernames{1}='Wavelength (nm)';
headernames2{1}='Wavelength (nm)';
headernames2{2}='Energy (eV)';

% Get data
for j=1:length(PLfiles)
    CurrentFile=strcat(pathname,PLfiles(j).name)
    NewData=dlmread(CurrentFile,'\t',90,0);
    smoothfactor=ceil(step/((max(NewData(:,1))-min(NewData(:,1))+1)/size(NewData,1)));
    yvals = interp1(NewData(:,1),smooth(NewData(:,2),smoothfactor),xvals);
    calib = interp1(SphereSpectrometerResponse(:,1),SphereSpectrometerResponse(:,2),xvals); % Correct for detector response
    
    PLdataRaw=[PLdataRaw, yvals];

    % Extract integration time
    fid = fopen(CurrentFile);
    temp=textscan(fid,'%s,\n','HeaderLines', 44);
    temp = temp{1};
    temp = temp{1};
    IntegerationTime=str2num(temp(18:end))/1000; % integration time in s
    fclose(fid);
        
    if subtractmin==1
        yvals=yvals-min(yvals);
    else
        yvals=yvals-subtractmin;
    end
    RawCountsPerSec=[RawCountsPerSec, yvals/IntegerationTime];
    PLdata=[PLdata, yvals./calib/IntegerationTime];
    %PLdataNorm=[PLdataNorm, yvals./calib/max(yvals./calib)];
    %PLdataeV=[PLdataeV, yvals./calib.*xvals.^2/1e10/IntegerationTime];
    %PLdataeVNorm=[PLdataeVNorm, yvals./calib.*xvals.^2/max(yvals./calib.*xvals.^2)];
    
    temp=PLfiles(j).name;
    headernames{j+1}=temp(1:end-9); % Names with extensions removed
    headernames2{j+2}=temp(1:end-9); % Names with extensions removed    
end

% write file
filepath=strcat(pathname,'SpherePLData.xls');
%xlswrite(filepath, headernames, 'RawData', 'A1')
%xlswrite(filepath, PLdataRaw, 'RawData', 'A2')
xlswrite(filepath, headernames, 'RawCountsPerSec', 'A1')
xlswrite(filepath, RawCountsPerSec, 'RawCountsPerSec', 'A2')
xlswrite(filepath, headernames, 'Sphere Response Corrected', 'A1')
xlswrite(filepath, PLdata, 'Sphere Response Corrected', 'A2') 
%xlswrite(filepath, headernames, 'Response Corrected Norm', 'A1')
%xlswrite(filepath, PLdataNorm, 'Response Corrected Norm', 'A2')
%xlswrite(filepath, headernames2, 'eV', 'A1')
%xlswrite(filepath, PLdataeV, 'eV', 'A2')
%xlswrite(filepath, headernames2, 'eV Normalized', 'A1')
%xlswrite(filepath, PLdataeVNorm, 'eV Normalized', 'A2')

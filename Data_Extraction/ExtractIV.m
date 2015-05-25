% Revised Oct 23, 2008 Eric Hoke

% This function takes all of the .txt files in a folder and generates
% a file IVdata.xls in the same folder in which each column contains the current 
% data for one of the files with a header which is the file name and
% a spreadsheet summarizing the Isc, Voc, FF, Pmax for each scan.
% When executed, the program will then prompt you to select a file 
% (doesn't matter which one) in the folder where you want to do the data
% extraction from.  This program assumes that all of the data files have 
% the exact same Voltage values data points, so REMOVE ANY FILES WHICH HAVE
% DIFFERENT VOLTAGE VALUES BEFORE RUNNING THE PROGRAM.  Any files with a
% different number of voltage points from the first file read will
% automatically be skipped.

% June 29, 2009 Added series and shunt resistance calculation
% Aug 31, 2009 Fixed bug with indicies when a file is of different size
% Nov 23, 2009 Added ideality factor calculation
% June 24, 2011 in summary, takes power conversion parameters from light
% curve and resistances, ideality factor from dark curve. For this to work,
% there must be a dark and light curve for each device, otherwise
% everything comes out wrong.
% July 31, 2013, skips over 'report.txt'

% Prompts user for a file, then it looks up all files in that folder with a
% *.txt extension
T=295; %Temperature (K) for ideality factor
Kb=1.381e-23; %Boltzmann constant in J/K
q=1.602e-19; % Elementary charge

[filename, pathname]=uigetfile('*.txt', 'Open any single .txt IV data file in the folder want to extract from');
IVfiles=dir(strcat(pathname, '*.txt'));
for n=1:length(IVfiles)
    IVfiles(n).name=strcat(pathname,IVfiles(n).name);
end

% Files should have same voltages for each IV scan.
% Extract voltage and current from first file
IVdata=dlmread(IVfiles(1).name,'\t', 27, 0);

%set up filename list
clear headernames;
headernames{1}='Voltage (V)';
temp=IVfiles(1).name;
headernames{2}=temp((length(pathname)+1):end);

% Get current data from other files  
% Does not add files which are different in length from first file
jj=0; % index to keep track of skipped files
for j=2:length(IVfiles)
    currentfile=IVfiles(j).name
    if strcmp('report.txt',currentfile((length(pathname)+1):end))
        jj=jj+1;
    else
        IVtemp=dlmread(IVfiles(j).name,'\t', 27, 0);
        if length(IVtemp) == length(IVdata(:,1))
            IVdata=cat(2,IVdata, IVtemp(:,2));
            headernames{j+1-jj}=currentfile((length(pathname)+1):end-4);
        else
            jj=jj+1;
        end
    end
end
NumData=size(IVdata,2)-1;

% Calculate Isc and Pmax
Isc = -1*interp1(IVdata(:,1),IVdata(:,2:end),0,'linear');
Pmax= max(-1*repmat(IVdata(:,1),1,NumData) .*IVdata(:,2:end),[],1);

Voc=0*Isc; %clear Voc vector
% find indicies of points where the current changes sign between
% consecutive points on IV curve
[Voc_ind,Sample_ind] = find(sign(IVdata(1:end-1,2:end) .* IVdata(2:end,2:end))==-1);

%slope of I-V curve crossing I=0 axis
m= diag((IVdata(Voc_ind+1,Sample_ind+1)-IVdata(Voc_ind,Sample_ind+1)))./(IVdata(Voc_ind+1,1)-IVdata(Voc_ind,1));
%interpolate voltage where I=0: V_int=V1-I1/m
Voc(Sample_ind)= IVdata(Voc_ind,1)-diag(IVdata(Voc_ind,Sample_ind+1))./m;

% Fill factor
FF=Pmax./(Isc.*Voc);

% Series Resistance
Rs=zeros(1,NumData); % series resistances in kohms-cm2 for all of the dark current curves
Vmax=max(IVdata(:,1));
Rs_limits=[Vmax-0.1 Vmax]; % Voltage range (V) to calculate series resistance
Rs_ind=find(IVdata(:,1)>=Rs_limits(1) & IVdata(:,1)<=Rs_limits(2)); % data range to compute series resistance
for j=1:NumData
    I=IVdata(:,j+1);
    fitcoef=polyfit(I(Rs_ind),IVdata(Rs_ind,1),1);     %Calculate series resistance Rs by linear fit between Rs_limits
    Rs(j)=fitcoef(1); % series resistance in kohms if current in mA
    if (Rs(j) > 1e6 || Rs(j) < 0) 
        Rs(j)=NaN; % if problem in fitting/ bad data, remove value
    end
end

% Shunt Resistance
Rsh=zeros(1,NumData); % Shunt resistances in kohms-cm2 for all of the dark current curves
Vmin=min(IVdata(:,1));
Rsh_limits=[Vmin Vmin+0.1]; % Voltage range (V) to calculate Shunt resistance
Rsh_ind=find(IVdata(:,1)>=Rsh_limits(1) & IVdata(:,1)<=Rsh_limits(2)); % data range to compute Shunt resistance
for j=1:NumData
    I=IVdata(:,j+1);
    fitcoef=polyfit(I(Rsh_ind),IVdata(Rsh_ind,1),1);     %Calculate Shunt resistance Rsh by linear fit between Rsh_limits
    Rsh(j)=fitcoef(1); % Shunt resistance in kohms if current in mA
    if (Rsh(j) < 0) 
        Rsh(j)=NaN; % if problem in fitting/ bad data, remove value
    end
end

n=NaN(1,NumData); % Ideality Factor
for j=find(~isnan(Rs))
    I=IVdata(:,j+1);
    Vc=IVdata(:,1)-Rs(j)*I; % Series resistance corrected voltage; this is to remove series resistance
    logI=log(abs(I)); 
    [logImin,ind] = min(logI); 
    FitInd= find(Vc>Vc(ind)+0.2); % Fit data starting at 0.2V higher than zero current crossing
    fitcoef=polyfit(logI(FitInd),Vc(FitInd),1);     %Calculate Ideality factor by linear fit of logI
    n(j)=fitcoef(1)*q/(Kb*T); % Ideality Factor
    n(n==0)=NaN; %Change zeros into NaN 
end

% Write out merged I-V file
filepath=strcat(pathname,'IVdata.xls');
xlswrite(filepath, headernames, 'IVdata', 'A1')
xlswrite(filepath, IVdata, 'IVdata', 'A2')

newheadernames{1}=' ';
for i=3:2:length(headernames)
   newheadernames{(i+1)/2} = headernames{i}(1:end-6);
end
xlswrite(filepath, newheadernames', 'Summary', 'A1')
xlswrite(filepath, {'Isc (mA/cm2)','Voc (V)','FF','Pmax (mW/cm2)', 'Rseries (ohm-cm2)','Rshunt (ohm-cm2)', 'Ideality Factor'}, 'Summary', 'B1')
xlswrite(filepath, vertcat(Isc(2:2:end),Voc(2:2:end),FF(2:2:end),Pmax(2:2:end),Rs(1:2:end)*1000,Rsh(1:2:end)*1000,n(1:2:end))', 'Summary', 'B2')

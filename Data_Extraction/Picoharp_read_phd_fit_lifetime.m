% Edited by Eric Hoke 10/09/08 and 12/18/08 to open all .phd files (from time resolved PL
% measurements in Brongersma lab) in a folder and
% export combined data as a .xls file.  The program can process files which have
% different lengths and different time resolutions. When you run the program select any file in the folder which
% you want to perform the data extraction.  Generated file is in same
% directory. Time=0 corresponds to the maximum of the PL time trace for
% each curve.

% 9/2/13 Eric Hoke. modified program so that files can have different lengths and time
% resolutions. Program shifts data curves so time=0 corresponds to the PL
% max. Program also fits first and last portion of data to an exponential
% to determine range of lifetimes. The results for each curve are plotted
% and the fitting parameters are written to the excel file.

function read_phd_Eric3(shouldNormalize,bgValue,final_res)

% CHANGE THESE VALUES TO CHANGE EXPONENTIAL FITTING RANGES
threshold=1e-3; % don't fit data which is below this threshold value of the maximum number of counts
fitfirstns=20; % (ns) fit the first few ns after maximum for initial decay
decayValue=1e-2; % fit until decay to this value
endstart=0.7; %fit this from fraction of total data to end for the final decay lifetime fit #3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%final_res=0.032; % final resolution in ns
fitfirst=round(fitfirstns/final_res);

[filename, pathname]=uigetfile('*.phd', 'Interactive mode data:');

phdfiles=dir(strcat(pathname, '*.phd')); %List of all files to process

%initialize arrays
alldata=[];
maxind=zeros(length(phdfiles),1);

for n=1:length(phdfiles)

% PicoHarp 300    File Access Demo in Matlab

% Demo access to binary PicoHarp 300 Data Files (*.phd)
% file format versions 1.0 and 1.1
% Read a PicoHarp data file and dump the contents in ASCII

% Tested with Matlab 5 and 6.
% Peter Kapusta, PicoQuant GmbH, February 2006
% This is demo code. Use at your own risk. No warranties.
% Make sure you have enough memory when loading large files!

fid=fopen(strcat(pathname,phdfiles(n).name));

fprintf(1,'\n=========================================================================== \n');
fprintf(1,'  Content of %s : \n', strcat(pathname, filename));
fprintf(1,'=========================================================================== \n');
fprintf(1,'\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% ASCII file header
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ident = char(fread(fid, 16, 'char'));
%fprintf(1,'               Ident: %s\n', Ident);

FormatVersion = deblank(char(fread(fid, 6, 'char')'));
%fprintf(1,'      Format version: %s\n', FormatVersion);

if not(strcmp(FormatVersion,'1.0') | strcmp(FormatVersion,'1.1'))
   fprintf(1,'\n\n      Warning: This program is for versions 1.0 and 1.1 only. Aborted.');
   STOP;
end;

CreatorName = char(fread(fid, 18, 'char'));
%fprintf(1,'        Creator name: %s\n', CreatorName);

CreatorVersion = char(fread(fid, 12, 'char'));
%fprintf(1,'     Creator version: %s\n', CreatorVersion);

FileTime = char(fread(fid, 18, 'char'));
%fprintf(1,'    Time of creation: %s\n', FileTime);

CRLF = char(fread(fid, 2, 'char'));

Comment = char(fread(fid, 256, 'char'));
%fprintf(1,'             Comment: %s\n', Comment);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Binary file header
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


NumberOfCurves = fread(fid, 1, 'int32');
%fprintf(1,'    Number of Curves: %d\n', NumberOfCurves);

BitsPerHistoBin = fread(fid, 1, 'int32');
%fprintf(1,'     Bits / HistoBin: %d\n', BitsPerHistoBin);

RoutingChannels = fread(fid, 1, 'int32');
%fprintf(1,'    Routing channels: %d\n', RoutingChannels);

NumberOfBoards = fread(fid, 1, 'int32');
%fprintf(1,'    Number of boards: %d\n', NumberOfBoards);

ActiveCurve = fread(fid, 1, 'int32');
%fprintf(1,'        Active Curve: %d\n', ActiveCurve);

MeasurementMode = fread(fid, 1, 'int32');
%fprintf(1,'    Measurement Mode: %d\n', MeasurementMode);

SubMode = fread(fid, 1, 'int32');
%fprintf(1,'            Sub-Mode: %d\n', SubMode);

RangeNo = fread(fid, 1, 'int32');
%fprintf(1,'            Range No: %d\n', RangeNo);

Offset = fread(fid, 1, 'int32');
%fprintf(1,'              Offset: %d\n', Offset);

Tacq = fread(fid, 1, 'int32');
%fprintf(1,'    Acquisition time: %d ms \n', Tacq);

StopAt = fread(fid, 1, 'int32');
%fprintf(1,'             Stop at: %d counts \n', StopAt);

StopOnOvfl = fread(fid, 1, 'int32');
%fprintf(1,'    Stop on Overflow: %d\n', StopOnOvfl);

Restart = fread(fid, 1, 'int32');
%fprintf(1,'             Restart: %d\n', Restart);

DispLinLog = fread(fid, 1, 'int32');
%fprintf(1,'     Display Lin/Log: %d\n', DispLinLog);

DispTimeAxisFrom = fread(fid, 1, 'int32');
%fprintf(1,'      Time Axis From: %d ns \n', DispTimeAxisFrom);

DispTimeAxisTo = fread(fid, 1, 'int32');
%fprintf(1,'        Time Axis To: %d ns \n', DispTimeAxisTo);

DispCountAxisFrom = fread(fid, 1, 'int32');
%fprintf(1,'     Count Axis From: %d\n', DispCountAxisFrom); 

DispCountAxisTo = fread(fid, 1, 'int32');
%fprintf(1,'       Count Axis To: %d\n', DispCountAxisTo);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:8
DispCurveMapTo(i) = fread(fid, 1, 'int32');
DispCurveShow(i) = fread(fid, 1, 'int32');
%fprintf(1,'-------------------------------------\n');
%fprintf(1,'            Curve No: %d\n', i-1);
%fprintf(1,'               MapTo: %d\n', DispCurveMapTo(i));
%fprintf(1,'                Show: %d\n', DispCurveShow(i));
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:3
ParamStart(i) = fread(fid, 1, 'float');
ParamStep(i) = fread(fid, 1, 'float');
ParamEnd(i) = fread(fid, 1, 'float');
%fprintf(1,'-------------------------------------\n');
%fprintf(1,'        Parameter No: %d\n', i-1);
%fprintf(1,'               Start: %d\n', ParamStart(i));
%fprintf(1,'                Step: %d\n', ParamStep(i));
%fprintf(1,'                 End: %d\n', ParamEnd(i));
end;
%fprintf(1,'-------------------------------------\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

RepeatMode = fread(fid, 1, 'int32');
%fprintf(1,'         Repeat Mode: %d\n', RepeatMode);

RepeatsPerCurve = fread(fid, 1, 'int32');
%fprintf(1,'      Repeat / Curve: %d\n', RepeatsPerCurve);

RepatTime = fread(fid, 1, 'int32');
%fprintf(1,'         Repeat Time: %d\n', RepatTime);

RepeatWaitTime = fread(fid, 1, 'int32');
%fprintf(1,'    Repeat Wait Time: %d\n', RepeatWaitTime);

ScriptName = char(fread(fid, 20, 'char'));
%fprintf(1,'         Script Name: %s\n', ScriptName);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%          Header for each board
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



for i = 1:NumberOfBoards
%fprintf(1,'-------------------------------------\n'); 
%fprintf(1,'           Board No: %d\n', i-1);

HardwareIdent(:,i) = char(fread(fid, 16, 'char'));
%fprintf(1,'Hardware Identifier: %s\n', HardwareIdent(:,i));

HardwareVersion(:,i) = char(fread(fid, 8, 'char'));
%fprintf(1,'   Hardware Version: %s\n', HardwareVersion(:,i));    
    
HardwareSerial(i) = fread(fid, 1, 'int32');
%fprintf(1,'   HW Serial Number: %d\n', HardwareSerial(i));

SyncDivider(i) = fread(fid, 1, 'int32');
%fprintf(1,'       Sync divider: %d \n', SyncDivider(i));

CFDZeroCross0(i) = fread(fid, 1, 'int32');
%fprintf(1,'    CFD0 zero cross: %d mV\n', CFDZeroCross0(i));

CFDLevel0(i) = fread(fid, 1, 'int32');
%fprintf(1,'         CFD0 level: %d mV\n', CFDLevel0(i));

CFDZeroCross1(i) = fread(fid, 1, 'int32');
%fprintf(1,'    CFD1 zero cross: %d mV\n', CFDZeroCross1(i));

CFDLevel1(i) = fread(fid, 1, 'int32');
%fprintf(1,'         CFD1 level: %d mV\n', CFDLevel1(i));

Resolution(i) = fread(fid, 1, 'float');
%fprintf(1,'         Resolution: %2.6g ns\n', Resolution(i));

end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                Headers for each histogram (curve)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:NumberOfCurves

fprintf(1,'-------------------------------------\n');
CurveIndex(i) = fread(fid, 1, 'int32');
fprintf(1,'        Curve Index: %d\n', CurveIndex(i));

TimeOfRecording(i) = fread(fid, 1, 'uint');

%  The PicoHarp software saves the time of recording
%  in a 32 bit serial time value as defined in all C libraries.
%  This equals the number of seconds elapsed since midnight
%  (00:00:00), January 1, 1970, coordinated universal time.
%  The conversion to normal date and time strings is tricky...

TimeOfRecording(i) = TimeOfRecording(i)/24/60/60+25569+693960;
fprintf(1,'  Time of Recording: %s \n', datestr(TimeOfRecording(i),'dd-mmm-yyyy HH:MM:SS'));

HardwareIdent(:,i) = char(fread(fid, 16, 'char'));
fprintf(1,'Hardware Identifier: %s\n', HardwareIdent(:,i));
    
HardwareVersion(:,i) = char(fread(fid, 8, 'char'));
fprintf(1,'   Hardware Version: %s\n', HardwareVersion(:,i));    
    
HardwareSerial(i) = fread(fid, 1, 'int32');
fprintf(1,'   HW Serial Number: %d\n', HardwareSerial(i));

SyncDivider(i) = fread(fid, 1, 'int32');
fprintf(1,'       Sync divider: %d \n', SyncDivider(i));

CFDZeroCross0(i) = fread(fid, 1, 'int32');
fprintf(1,'    CFD0 zero cross: %d mV\n', CFDZeroCross0(i));

CFDLevel0(i) = fread(fid, 1, 'int32');
fprintf(1,'         CFD0 level: %d mV\n', CFDLevel0(i));

CFDZeroCross1(i) = fread(fid, 1, 'int32');
fprintf(1,'    CFD1 zero cross: %d mV\n', CFDZeroCross1(i));

CFDLevel1(i) = fread(fid, 1, 'int32');
fprintf(1,'         CFD1 level: %d mV\n', CFDLevel1(i));

Offset(i) = fread(fid, 1, 'int32');
fprintf(1,'             Offset: %d\n', Offset(i));

RoutingChannel(i) = fread(fid, 1, 'int32');
fprintf(1,'    Routing channel: %d\n', RoutingChannel(i));

ExtDevices(i) = fread(fid, 1, 'int32');
fprintf(1,'   External Devices: %d \n', ExtDevices(i));

MeasMode(i) = fread(fid, 1, 'int32');
fprintf(1,'   Measurement Mode: %d\n', MeasMode(i));

SubMode(i) = fread(fid, 1, 'int32');
fprintf(1,'           Sub-Mode: %d\n', SubMode(i));

P1(i) = fread(fid, 1, 'float');
fprintf(1,'                 P1: %d\n', P1(i));
P2(i) = fread(fid, 1, 'float');
fprintf(1,'                 P2: %d\n', P2(i));
P3(i) = fread(fid, 1, 'float');
fprintf(1,'                 P3: %d\n', P3(i));

RangeNo(i) = fread(fid, 1, 'int32');
fprintf(1,'          Range No.: %d\n', RangeNo(i));

Resolution(n) = fread(fid, 1, 'float');
fprintf(1,'         Resolution: %2.6g ns \n', Resolution(i));

Channels(i) = fread(fid, 1, 'int32');
fprintf(1,'           Channels: %d \n', Channels(i));

Tacq(i) = fread(fid, 1, 'int32');
fprintf(1,'   Acquisition Time: %d ms \n', Tacq(i));

StopAfter(i) = fread(fid, 1, 'int32');
fprintf(1,'         Stop After: %d ms \n', StopAfter(i));

StopReason(i) = fread(fid, 1, 'int32');
fprintf(1,'        Stop Reason: %d\n', StopReason(i));

InpRate0(i) = fread(fid, 1, 'int32');
fprintf(1,'       Input Rate 0: %d Hz\n', InpRate0(i));

InpRate1(i) = fread(fid, 1, 'int32');
fprintf(1,'       Input Rate 1: %d Hz\n', InpRate1(i));

HistCountRate(i) = fread(fid, 1, 'int32');
fprintf(1,'   Hist. Count Rate: %d cps\n', HistCountRate(i));

IntegralCount(i) = fread(fid, 1, 'int64');
fprintf(1,'     Integral Count: %d\n', IntegralCount(i));

Reserved(i) = fread(fid, 1, 'int32');
fprintf(1,'           Reserved: %d\n', Reserved(i));

DataOffset(i) = fread(fid, 1, 'int32');
fprintf(1,'         DataOffset: %d\n', DataOffset(i));

end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%          Reads all histograms into one matrix
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:NumberOfCurves
    fseek(fid,DataOffset(i),'bof');
    Counts(:,i) = fread(fid, Channels(i), 'uint32');
end;

Peak=max(Counts);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%          Summary
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf(1,'\n');
fprintf(1,'\n');
fprintf(1,'=====================================================\n');
fprintf(1,'                     SUMMARY                         \n');
fprintf(1,'=====================================================\n');
fprintf(1,' Curve    Channel     Number of    Peak     Integral \n');
fprintf(1,' index   resolution   channels     count     count   \n');
fprintf(1,'=====================================================\n');

for i = 1:NumberOfCurves
fprintf(1,'  %3i       %2.6g  %10i  %10i  %10i\n', CurveIndex(i),Resolution(i), Channels(i), Peak(i), IntegralCount(i));   
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%          This is a simple display of the histogram(s)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% figure(1);
% semilogy(Counts);
% % axis([0 max(max(Channels)) 1 10*max(max(Counts))]);
% xlabel('Channel #');
% ylabel('Counts');
% 
% if NumberOfCurves<21
%    legend(num2str((1:NumberOfCurves)'),0);
% end;

fclose(fid);

% Average data to final_res
aver=round(final_res/Resolution(n));
newNumPts=floor(length(Counts)/aver);
tempdata=NaN(newNumPts,1);
for i= 1:newNumPts
    tempdata(i)=mean(Counts(aver*(i-1)+1:aver*i));
end

%remove any trailing zeros and last non-zero data point
tempdata=tempdata(find(tempdata, 1, 'first'):find(tempdata, 1, 'last'));
tempdata=tempdata(1:end-1);
%normalize
maxind(n)=find(tempdata==max(tempdata),1,'last');
if shouldNormalize
    tempdata=(tempdata-bgValue)/(max(tempdata)-bgValue);
else
    tempdata=tempdata-bgValue;
end
% add to data sheet
alldata=spreadsheetcat(alldata,tempdata);

end

% Shift data so that time=0 is at the maximum PL signal
maxshift=max(maxind)-1;
alldata2=NaN(size(alldata,1)+maxshift+1-min(maxind),size(alldata,2)+1);
alldata2(:,1)= ((1:size(alldata2,1))-maxshift-1)*final_res; %time values (in ns)
for n=1:length(phdfiles)
    alldata2(maxshift-maxind(n)+2:maxshift+1-maxind(n)+size(alldata,1),n+1)=alldata(:,n);
end

% construct header name list
headernames{1}='Time (ns)';
for n=1:length(phdfiles)
    temp=phdfiles(n).name;
    headernames{n+1}=temp(1:end-4); %file name with extension removed
end

% Perform exponential fit at beginning and end of data trace
logp1=zeros(length(phdfiles),2);
logp2=logp1;
logp3=logp1;
index_max=max(maxind);
for n=1:length(phdfiles)

    % Determine fitting range
    
    % Peform fit of first few ns
    fitTimes=alldata2(index_max:index_max+fitfirst,1);
    logp1(n,:) = polyfit(fitTimes,log(alldata2(index_max:index_max+fitfirst,n+1)),1);

    % Perform fit on decay until reaches certain level
    temp=alldata2(index_max:end,n+1);
    tempind=find([temp' 0] < max(temp)*decayValue,1);
    tempend=min([index_max+tempind-1 length(temp) find(alldata2(:,n+1)>threshold,1,'last')]);
    fitTimes3=alldata2(index_max:tempend,1);
    logp3(n,:) = polyfit(fitTimes3,log(alldata2(index_max:tempend,n+1)),1);
    
    % Peform fit of end of data trace above threshold
    stop=find(alldata2(:,n+1)>threshold,1,'last');
    start=round(endstart*stop);
    fitTimes2=alldata2(start:stop,1);
    logp2(n,:) = polyfit(fitTimes2,log(alldata2(start:stop,n+1)),1);    

    % Plot Data and fit
    figure(n)
    semilogy(alldata2(1:stop,1),alldata2(1:stop,n+1),'k')
    hold on 
    semilogy(fitTimes,exp(logp1(n,2)+logp1(n,1)*fitTimes),'r')
    semilogy(fitTimes2,exp(logp2(n,2)+logp2(n,1)*fitTimes2),'r')
    semilogy(fitTimes3,exp(logp3(n,2)+logp3(n,1)*fitTimes3),'r')
    xlabel('time (ns)')
    ylabel('PL (a.u.)')
    title(headernames{n+1})
    legend('data',strcat('t_0=',num2str(-1/logp1(n,1)),' ns'), strcat('t_n=',num2str(-1/logp2(n,1)),' ns'),strcat('t_0.01=',num2str(-1/logp3(n,1)),' ns'))
    hold off    
end
logp1(:,1)=-1./logp1(:,1); % lifetimes generated
logp2(:,1)=-1./logp2(:,1); % lifetimes generated
logp3(:,1)=-1./logp3(:,1); % lifetimes generated
logp1(:,2)=exp(logp1(:,2)); % prefactor
logp2(:,2)=exp(logp2(:,2)); % prefactor
logp3(:,2)=exp(logp3(:,2)); % prefactor

% Write files
filepath=strcat(pathname,'phd_file_data.xls');
xlswrite(filepath, headernames, 'TCSPC Data', 'A1')
xlswrite(filepath, alldata2, 'TCSPC Data', 'A2')
xlswrite(filepath, headernames, 'LifeTimeFit', 'A1')
xlswrite(filepath, {'a0*exp(-t/t0)';'t0';'a0';'tn';'an';'t3';'a3'}, 'LifeTimeFit', 'A1')
xlswrite(filepath, logp1', 'LifeTimeFit', 'B2')
xlswrite(filepath, logp2', 'LifeTimeFit', 'B4')
xlswrite(filepath, logp3', 'LifeTimeFit', 'B6')
 

function C=spreadsheetcat(A,B)
C=NaN(max(size(A,1),size(B,1)),size(A,2)+size(B,2));
C(1:size(A,1),1:size(A,2))=A;
C(1:size(B,1),size(A,2)+1:end)=B;

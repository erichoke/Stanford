% adapted from ExtractIV.m 3/16/11
% Interpolates data to xvals
% file starting at cell RowBegin, ColumnBegin (use 0 for no header text cells).  
function ExtractInterpolate(RowBegin,ColumnBegin,xvals,delim)
% for EQE RowBegin=1, ColumnBegin=0
% for EQEV RowBegin=4, ColumnBegin=0
% for PL ExtractInterpolate(90,0,xvals,'\t')

% Prompts user for a file, then it looks up all files in that folder
[filename, pathname]=uigetfile('*', 'Open any file in folder you want to extract from');
IVfiles=dir(strcat(pathname, '*'));

% list of all files in folder which are not directories
IVfiles=IVfiles([IVfiles.isdir]==0);

% Get data
IVdata=xvals';
headernames{1}='Independent Variable';
for j=1:length(IVfiles)
    CurrentFile=strcat(pathname,IVfiles(j).name)
    IVtemp=dlmread(CurrentFile,delim,RowBegin,ColumnBegin);
    yvals = interp1(IVtemp(:,1),IVtemp(:,2),xvals)';
    IVdata=[IVdata, yvals];
    
    temp=IVfiles(j).name;
    temp1=[find(temp=='.') length(temp)+1] -1;
    temp2=temp1(max(end-1,1));
    headernames{j+1}=temp(1:temp2); % Names with extensions removed
end

filepath=strcat(pathname,'Data.xls');
xlswrite(filepath, headernames, 'Data', 'A1')
xlswrite(filepath, IVdata, 'Data', 'A2') 

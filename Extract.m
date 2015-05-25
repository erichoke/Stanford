% adapted from ExtractIV.m 10/16/10
% Dumps first NumColumns in every file in a folder into one excel
% file starting at cell RowBegin, ColumnBegin (use 0 for no header text cells).  
function Extract(RowBegin,ColumnBegin,NumColumns,delim)
% for EQE: Extract(1,0,2,'\t') 
% for EQEV: Extract(4,0,2,'\t')

% Prompts user for a file, then it looks up all files in that folder
[filename, pathname]=uigetfile('*', 'Open any file in folder you want to extract from');
IVfiles=dir(strcat(pathname, '*'));

% list of all files in folder which are not directories
IVfiles=IVfiles([IVfiles.isdir]==0);

% Get data
IVdata=[];
headernames='';
for j=1:length(IVfiles)
    IVfiles(j).name
    IVtemp=dlmread(strcat(pathname,IVfiles(j).name),delim,RowBegin,ColumnBegin);
    IVdata=spreadsheetcat(IVdata, IVtemp(:,1:NumColumns));
    
    temp=IVfiles(j).name;
    temp=IVfiles(j).name;
    temp1=[find(temp=='.') length(temp)+1] -1;
    temp2=temp1(max(end-1,1));
    for k=1:NumColumns
        headernames=strcat(headernames, temp(1:temp2), ',');
    end
end

filepath=strcat(pathname,'Data.csv');
dlmwrite(filepath,headernames, '');
dlmwrite(filepath,IVdata,'-append');         
% append option does not work for older versions of MatLab

% puts matrix B horizontally after A, fills in empty spaces with NaN
function C=spreadsheetcat(A,B)
C=NaN(max(size(A,1),size(B,1)),size(A,2)+size(B,2));
C(1:size(A,1),1:size(A,2))=A;
C(1:size(B,1),size(A,2)+1:end)=B;    
%Retrieving the names of the excel files
dirdata = dir('schmidt2010.xls'); % file name of IV data

% Tao folder chua du lieu ra (neu folder chua ton tai)
mkdir('DataOutput');

% Tao folder chua du lieu cua phuong phap ortizconde
mkdir('DataOutput', 'OrtizcondeMethode'); 

for h = 1:length(dirdata)
    dirname = dirdata(h).name;
    dirnamef = dirname;
    dirname = dirname(1:length(dirname) - 4);
    
    mkdir('DataOutput\OrtizcondeMethode\', dirname);
    savepath = ['DataOutput\OrtizcondeMethode\' dirname '\'];
    excelname = [savepath 'Data' dirname '.xls'];
    rawdata = [savepath dirname '.xls'];
    
    lightfull = xlsread(dirnamef,'Light'); % Lay so lieu cua duong IV duoi su chieu sang
    datasize = size(lightfull); % Lay kich thuoc ma tran lightfull
    dnumber = datasize(2)/2; % Tinh so duong IV do duoi su chieu sang
    
    writematrix = [];
    for k = 0:(dnumber-1)
        light = [lightfull(:,(1+2*k)), lightfull(:,(2+2*k))]; % Chon duong IV thu k
        
        %Calling the main function to calculate solar cell data
        outdata = [k; Ortizconde(light)];
        
        %Writing to the write matrix
        writematrix = [writematrix outdata];
    end
    %Writing to an excel file
    valmatrix = []; %Write matrix used to write to an excel file
            i) valmatrix = [{'I-V curve'}; {'Rs'}; {'Rsh'}; {'FF'}; {'PCE'}];
    
    xlswrite(excelname, valmatrix, dirname);
    xlswrite(excelname, writematrix, dirname, 'B1');
    
    DeleteEmptyExcelSheets([cd '\' excelname]);
end

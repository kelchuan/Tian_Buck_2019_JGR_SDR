
%% Initialize variables.
filename = '/Users/Tian/Documents/MATLAB/SDR/Image_J/Workbook2.csv';
delimiter = ',';
formatSpec = '%f%f%f%f%f%f%f%[^\n\r]';
%% Open the text file.
fileID = fopen(filename,'r');
%% Read columns of data according to format string.
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter,  'ReturnOnError', false);
%% Close the text file.
fclose(fileID);
%% Allocate imported array to column variable names
VarName1 = dataArray{:, 1};
VarName2 = dataArray{:, 2};
VarName3 = dataArray{:, 3};
VarName4 = dataArray{:, 4};
VarName5 = dataArray{:, 5};
VarName6 = dataArray{:, 6};
VarName7 = dataArray{:, 7};
%% Clear temporary variables
clearvars filename delimiter formatSpec fileID dataArray ans;


yy = VarName7; %km
xx = VarName6; %km
xx_t = xx * 0.1; %Myr    0.1 for 1cm/yr of Vx
diff_yy_xx = diff(yy)./diff(xx_t); %km/Myr
Vs_flac = diff_yy_xx * 0.1; %cm/Myr
Vs_flac = [Vs_flac' 0];

hold on
plot(xx_t * 1000, Vs_flac, 'ko')



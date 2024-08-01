%write in trc file

% marker_data=all_data.data;
% C=split(all_data.textdata{3,1},',');
% 
% markers=[];
% for i=3:3:size(C,1)
%     markers=[markers; {strrep(C{i},'MarkerProtocol:','')}];
% end

nframes=size(marker_data,1);
nframe=1:nframes;
time=0:0.01:(nframes-1)/100;
nmarkers=size(markers_OpenSim,2);
% first initialise the header with a column for the Frame # and the Time
% also initialise the format for the columns of data to be written to file
dataheader1 = 'Frame#\tTime\t';
dataheader2 = '\t\t';
format_text = '%i\t%2.4f\t';
% initialise the matrix that contains the data as a frame number and time row
% data_out = [nframe; time]';
data_out =[];
data_struct.Start_Frame=nframe(1);
data_struct.End_Frame=nframe(end);
data_struct.Rate=100;
nrows=nframes;
data_struct.units='mm';
 
%% Correct it with the last unlabeled marker (or markers) does not have data
if ((size(marker_data,2)-2)/3)~=nmarkers
    nmissing=nmarkers-(size(marker_data,2)-2)/3;
    marker_data=[marker_data zeros(nframes,nmissing*3)];
end
    
    
%% now loop through each maker name and make marker name with 3 tabs for the
% first line and the X Y Z columns with the marker number on the second
% line all separated by tab delimeters
% each of the data columns (3 per marker) will be in floating format with a
% tab delimiter - also add to the data matrix
for ii = 1:nmarkers
    dataheader1 = [dataheader1 markers_OpenSim{ii} '\t\t\t'];    
    dataheader2 = [dataheader2 'X' num2str(ii) '\t' 'Y' num2str(ii) '\t'...
        'Z' num2str(ii) '\t'];
    format_text = [format_text '%f\t%f\t%f\t'];
%     data_out=[data_out marker_data(:,(ii*3-1):(ii*3+1))];
end
% marker_data(:,3:end)=marker_data(:,3:end)*1000;
data_out=marker_data;


dataheader1 = [dataheader1 '\n'];
dataheader2 = [dataheader2 '\n'];
format_text = [format_text '\n'];

disp('Writing trc file...') 

newfilename = filename;

%open the file
fid_1 = fopen([newfilename],'w');

% first write the header data
fprintf(fid_1,'PathFileType\t4\t(X/Y/Z)\t %s\n',newfilename);
fprintf(fid_1,'DataRate\tCameraRate\tNumFrames\tNumMarkers\tUnits\tOrigDataRate\tOrigDataStartFrame\tOrigNumFrames\n');
fprintf(fid_1,'%d\t%d\t%d\t%d\t%s\t%d\t%d\t%d\n', data_struct.Rate, data_struct.Rate, nrows, nmarkers, data_struct.units, data_struct.Rate,data_struct.Start_Frame,data_struct.End_Frame); 
fprintf(fid_1, dataheader1);
fprintf(fid_1, dataheader2);

% then write the output marker data
fprintf(fid_1, format_text,data_out');

% close the file
fclose(fid_1);

disp('Done.')
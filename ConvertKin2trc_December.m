filepath=fileparts(matlab.desktop.editor.getActiveFilename)
cd(filepath);

markers={'pelvis_top','hip','pelvis_bottom','knee','ankle','mtp','toe','FT_origin','ref1','ref2','ref3'};
markers_smoothed={'filtered_smoothed_pelvisTop','filtered_smoothed_hip','filtered_smoothed_pelvisBottom','filtered_smoothed_knee','filtered_smoothed_ankle','filtered_smoothed_MTP','filtered_smoothed_toe','filtered_smoothed_origin','filtered_smoothed_interface3','filtered_smoothed_interface4','filtered_smoothed_ref3'};
markers_OpenSim={'pelvis_top','hip','pelvis_bottom','knee','ankle','mtp','toe','FT_origin','ref1','ref2','ref3'};
main_folder=pwd;
addpath(main_folder);

cd('DataDecember/kinematics');

filenames=dir('*.csv');
for j=1:length(filenames)
    tablein=readtable(filenames(j).name);
%         tablein_trans=rows2vars(tablein);
    datam=[];
    % if contains(folders(i).name,{'h65','h7','h8','h9'})
    %     tablein_trans=rows2vars(tablein);
    %     timeI=find(contains(table2cell(tablein_trans(1,:)),'relative_time'));
    %     time=cell2mat(table2array(tablein_trans(2:end,timeI)));
    %     for m=1:length(markers_smoothed)
    %         datam=[datam readMarker_smoothed(tablein_trans,markers_smoothed{m})];
    %     end
    % else
        time=tablein.relative_time;
        for m=1:length(markers)
            datam=[datam readMarker(tablein,markers_smoothed{m})];
        end
    % end
    %remove repeated time stamps
    [time,datam]=RemoveRepeatedTimeStamps(time,datam);
    
    datam(:,2:3:end)=-datam(:,2:3:end); %change y coordinate to work as if it was a right foot

    nframes=size(datam,1);
    frames=0:nframes-1;
    marker_data=[frames' time datam];
    filename=strrep(filenames(j).name,'.csv','.trc');
    write_trcfile;
end

function out=readMarker(tablein,marker)
%     dataI=find(contains(table2cell(tablein_trans(1,:)),marker))

I=find(contains(tablein.Properties.VariableNames,marker));
if length(I)>3
    jj=[];
    for j=1:length(I)
        if strcmp(tablein.Properties.VariableNames(I(j)),[marker '_x'])|...
                strcmp(tablein.Properties.VariableNames(I(j)),[marker '_y'])|...
                strcmp(tablein.Properties.VariableNames(I(j)),[marker '_z']);
        else
            jj=[jj j];
        end
    end
    I(jj)=[];
end
    out=table2array(tablein(:,I:(I+2)));
end
function out=readMarker_smoothed(tablein,marker)

if contains(marker,'origin')
    I=find(contains(table2cell(tablein(1,:)),'origin'));
else
    I=find(contains(table2cell(tablein(1,:)),marker));
end
    out=cell2mat(table2array(tablein(2:end,I)));
end
function  [time,datam]=RemoveRepeatedTimeStamps(time,datam);

II=find(time(2:end)==time(1:end-1));
time(II)=[];
datam(II,:)=[];
end
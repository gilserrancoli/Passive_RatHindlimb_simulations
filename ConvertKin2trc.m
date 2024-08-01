filepath=fileparts(matlab.desktop.editor.getActiveFilename)
cd(filepath);

markers={'smoothed_pelvis top','smoothed_hip','smoothed_pelvis bottom','smoothed_knee','smoothed_ankle','smoothed_MTP','smoothed_toe','smoothed_FT_origin','smoothed_ref1','smoothed_ref2','smoothed_ref3'};
markers_OpenSim={'pelvis_top','hip','pelvis_bottom','knee','ankle','mtp','toe','FT_origin','ref1','ref2','ref3'};
main_folder=pwd;
addpath(main_folder);

cd('kinematics');
folders=dir();
folders=folders(~ismember({folders.name},{'.','..'}));
for i=2:length(folders)
    cd(folders(i).name);
    filenames=dir('*.mat');
    for j=1:length(filenames)
        load(filenames(j).name,'testTable');
        time=testTable.relative_time;
        datam=[];
        for m=1:length(markers)
            datam=[datam readMarker(testTable,markers{m})];
        end
        %remove repeated time stamps
        [time,datam]=RemoveRepeatedTimeStamps(time,datam);
        if contains(folders(i).name,["forward","backward"])
            datam=-datam; %change x and z to adapt kinematics to force coordinate system, and y to work as if it was a right foot
%             datam(:,2:3:end)=-datam(:,2:3:end); %change y coordinate to work as if it was a right foot
        end
        
        nframes=size(datam,1);
        frames=0:nframes-1;
        marker_data=[frames' time datam];
        filename=strrep(filenames(j).name,'.mat','.trc');
        write_trcfile;
    end
end
function out=readMarker(testTable,markers)
    out=[testTable.([markers '_x']) ...
        testTable.([markers '_y']) ...
        testTable.([markers '_z'])];
end
function  [time,datam]=RemoveRepeatedTimeStamps(time,datam);

II=find(time(2:end)==time(1:end-1));
time(II)=[];
datam(II,:)=[];
end
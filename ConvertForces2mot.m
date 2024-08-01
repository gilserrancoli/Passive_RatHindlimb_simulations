filepath=fileparts(matlab.desktop.editor.getActiveFilename)
cd(filepath);

colheaders={'time','FT_vx','FT_vy','FT_vz','FT_px','FT_py','FT_pz','FT_tx','FT_ty','FT_tz'};
cd('motion_force');
folders=dir();
folders=folders(~ismember({folders.name},{'.','..'}));
for ii=2:length(folders)
    cd(folders(ii).name);
    filenames=dir('*.mat');
    for j=1:length(filenames)
        load(filenames(j).name);
        dirkinfiles=dir(['..\..\kinematics\' strrep(folders(ii).name,'_mean_response','') '\*.mat']); 
        kinfilename=IdentifyKinFilename(dirkinfiles,filenames(j).name);
        kinfiledata=load([dirkinfiles(1).folder '\' kinfilename]);
        origin_FT=[kinfiledata.testTable.smoothed_FT_origin_x ...
            kinfiledata.testTable.smoothed_FT_origin_y ...
            kinfiledata.testTable.smoothed_FT_origin_z]/1000;
        kintime=kinfiledata.testTable.relative_time;
        %remove repeated time stamps
        [kintime,origin_FT]=RemoveRepeatedTimeStamps(kintime,origin_FT);
        
        time=testTable.relative_time; %% force time
        origin_FT_rs=interp1(kintime,origin_FT,time);
        %change sign in force and CoP y axis to work as if it was the right
        %foot (note that this implies to change x and z axes for torque,
        %not y)
        testTable.baseremoved_Fy=-testTable.baseremoved_Fy;
        origin_FT_rs=-origin_FT_rs; %the origin comes from kinematics, so we need to change both x and z to adapt it to force coordinate system, and y to work as if it was right leg
        testTable.baseremoved_Tx=-testTable.baseremoved_Tx;
        testTable.baseremoved_Tz=-testTable.baseremoved_Tz;
        
        data=[time testTable.baseremoved_Fx testTable.baseremoved_Fy testTable.baseremoved_Fz...
            origin_FT_rs testTable.baseremoved_Tx testTable.baseremoved_Ty testTable.baseremoved_Tz];        
        II=isnan(data);
        indnan=any(II,2);
        data(indnan,:)=[];
        q.data=data;
        q.labels=colheaders;
        forcefilename=strrep(filenames(j).name,'.mat','.mot');
        write_motionFile(q,forcefilename);
    end
end

function  kinfilename=IdentifyKinFilename(dirkinfiles,forcefilename);
C=strsplit(forcefilename,'_');
id=C{2}(8:end);

found=false;
i=1;
while (i<=length(dirkinfiles))&(~found)
    namei=dirkinfiles(i).name;
    C=strsplit(namei,'_');
    found=strcmp(C{2},['coordinates' id]);
    if found
        kinfilename=dirkinfiles(i).name;
    else
        i=i+1;
    end
end

end
function  [time,datam]=RemoveRepeatedTimeStamps(time,datam);

II=find(time(2:end)==time(1:end-1));
time(II)=[];
datam(II,:)=[];
end
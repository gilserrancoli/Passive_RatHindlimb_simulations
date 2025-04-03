filepath=fileparts(matlab.desktop.editor.getActiveFilename);
cd(filepath);

colheaders={'time','FT_vx','FT_vy','FT_vz','FT_px','FT_py','FT_pz','FT_tx','FT_ty','FT_tz'};
cd('DataMarch2025/achillescut_ForwardOnly/perturbation');

filenames=dir('*.csv');
for j=1:length(filenames)
    %load force data
    tablein=readtable(filenames(j).name);
    % tablein_trans=rows2vars(tablein);
   
    %load kinematics data
    kinfilename=['..\kinematics\' strrep(strrep(filenames(j).name,'.csv','.trc'),'motor','kinematics')];
    kinfiledata=readtable(kinfilename,'FileType','Text','NumHeaderLines',3);
    if isnan(kinfiledata.Frame_(1))
        kinfiledata(1,:)=[];
    end
    FT_originI=find(contains(kinfiledata.Properties.VariableNames,'FT_origin'));
    origin_FT=table2array(kinfiledata(:,FT_originI:FT_originI+2))/1000;
    kintime=kinfiledata.Time;
    %remove repeated time stamps
    [kintime,origin_FT]=RemoveRepeatedTimeStamps(kintime,origin_FT);
    
    %continue gathering force data
    timeI=find(contains(tablein.Properties.VariableNames,'relative_time'));
    time=table2array(tablein(:,timeI)); %% force time
    t0=max(time(1),kintime(1));
    tf=min(time(end),kintime(end));
%         time(time<t0)=[];
%         time(time>tf)=[];
    origin_FT_rs=interp1(kintime,origin_FT,time);
    
    %change sign in force and CoP y axis to work as if it was the right
    %foot (note that this implies to change x and z axes for torque,
    %not y). CoP already comes from kinematics (sign ahs already been changes. 
    %Change the sign to x, y and z to forces and moments, since these
    %will be applied to the animal (not the TF as measured)
    FxI=find(contains(tablein.Properties.VariableNames(1,:),'baseremoved_Fx'));
    FyI=find(contains(tablein.Properties.VariableNames(1,:),'baseremoved_Fy'));
    FzI=find(contains(tablein.Properties.VariableNames(1,:),'baseremoved_Fz'));
    TxI=find(contains(tablein.Properties.VariableNames(1,:),'baseremoved_Tx'));
    TyI=find(contains(tablein.Properties.VariableNames(1,:),'baseremoved_Ty'));
    TzI=find(contains(tablein.Properties.VariableNames(1,:),'baseremoved_Tz'));
    
    Fx=-table2array(tablein(:,FxI)); %change sign due to action/reaction
    Fy=-table2array(tablein(:,FyI)); %change sign due to action/reaction
    Fz=-table2array(tablein(:,FzI)); %change sign due to action/reaction
    Tx=table2array(tablein(:,TxI)); %double change sign due to action/reaction and change the y axis direction
    Ty=-table2array(tablein(:,TyI)); %change sign due to action/reaction
    Tz=table2array(tablein(:,TzI)); %double change sign due to action/reaction and change the y axis direction

    data=[time Fx Fy Fz origin_FT_rs Tx Ty Tz];        
    II=isnan(data);
    indnan=any(II,2);
    data(indnan,:)=[];
    q.data=data;
    q.labels=colheaders;
    forcefilename=strrep(filenames(j).name,'.csv','.mot');
    write_motionFile(q,forcefilename);
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
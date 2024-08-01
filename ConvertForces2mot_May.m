filepath=fileparts(matlab.desktop.editor.getActiveFilename);
cd(filepath);

colheaders={'time','FT_vx','FT_vy','FT_vz','FT_px','FT_py','FT_pz','FT_tx','FT_ty','FT_tz'};
cd('DataMay/motion&force');
folders=dir();
folders=folders(~ismember({folders.name},{'.','..'}));
for i=length(folders):-1:1
    if isfolder(folders(i).name)
    else
        folders(i)=[];
    end
end
for ii=1:length(folders)
    cd([folders(i).folder '\' folders(ii).name]);
    filenames=dir('*.csv');
    for j=1:length(filenames)
        %load force data
        tablein=readtable(filenames(j).name);
        tablein_trans=rows2vars(tablein);
       
        %load kinematics data
        kinfilename=['..\..\kinematics\' folders(ii).name '\' strrep(filenames(j).name,'.csv','.trc')];
        kinfiledata=readtable(kinfilename,'FileType','Text');
        FT_originI=find(contains(kinfiledata.Properties.VariableNames,'FT_origin'));
        origin_FT=table2array(kinfiledata(:,FT_originI:FT_originI+2))/1000;
        kintime=kinfiledata.Time;
        %remove repeated time stamps
        [kintime,origin_FT]=RemoveRepeatedTimeStamps(kintime,origin_FT);
        
        %continue gathering force data
        timeI=find(contains(table2cell(tablein_trans(1,:)),'relative_time'));
        time=cell2mat(table2array(tablein_trans(2:end,timeI))); %% force time
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
        FxI=find(contains(table2cell(tablein_trans(1,:)),'baseremoved_Fx'));
        FyI=find(contains(table2cell(tablein_trans(1,:)),'baseremoved_Fy'));
        FzI=find(contains(table2cell(tablein_trans(1,:)),'baseremoved_Fz'));
        TxI=find(contains(table2cell(tablein_trans(1,:)),'baseremoved_Tx'));
        TyI=find(contains(table2cell(tablein_trans(1,:)),'baseremoved_Ty'));
        TzI=find(contains(table2cell(tablein_trans(1,:)),'baseremoved_Tz'));
        
        Fx=-cell2mat(table2array(tablein_trans(2:end,FxI))); %change sign due to action/reaction
        Fy=-cell2mat(table2array(tablein_trans(2:end,FyI))); %change sign due to action/reaction
        Fz=-cell2mat(table2array(tablein_trans(2:end,FzI))); %change sign due to action/reaction
        Tx=cell2mat(table2array(tablein_trans(2:end,TxI))); %double change sign due to action/reaction and change the y axis direction
        Ty=-cell2mat(table2array(tablein_trans(2:end,TyI))); %change sign due to action/reaction
        Tz=cell2mat(table2array(tablein_trans(2:end,TzI))); %double change sign due to action/reaction and change the y axis direction

        data=[time Fx Fy Fz origin_FT_rs Tx Ty Tz];        
        II=isnan(data);
        indnan=any(II,2);
        data(indnan,:)=[];
        q.data=data;
        q.labels=colheaders;
        forcefilename=strrep(filenames(j).name,'.csv','.mot');
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
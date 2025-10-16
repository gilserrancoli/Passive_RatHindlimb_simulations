baseline_folder='C:\Gil\Collaborations\MatthewTresch\Zhong_ParameterEstimation\Optimization\Datarat21 - DataMay2025\baseline_ForwardOnly\kinematics';
achillescut_folder='C:\Gil\Collaborations\MatthewTresch\Zhong_ParameterEstimation\Optimization\Datarat21 - DataMay2025\achillescut_ForwardOnly\kinematics';

cd(baseline_folder);
baseline_files=dir('*.trc');
data3d_baseline=GetData(baseline_files);
pertdata_baseline=GetPertData(baseline_files);
cd(achillescut_folder);
achillescut_files=dir('*.trc');
data3d_achillescut=GetData(achillescut_files);
pertdata_achillescut=GetPertData(baseline_files);


%Plot in 2D graphs
figure(1);
plot2d(data3d_baseline,[0 0 1]);
figure(1);
plot2d(data3d_achillescut,[1 0 0]);

%% Compute mean points
data3d_all=[data3d_baseline data3d_achillescut];
[alldata, mdata]=ConcatenateAll(data3d_all);

%% Plot in 2D graphs data concatenated
figure(2);
plot2d_concatenated(alldata);

%% Compute mean points by perturbation
[alldata_bypert, mdata_bypert]=ConcatenateAll_bypert(data3d_baseline,data3d_achillescut,baseline_folder,pertdata_baseline,pertdata_achillescut);


%% Compute mean points by perturbation


%Plot3D
figure(2)
plot3d(data3d_baseline);
plot3d(data3d_achillescut);


function out=GetData(files)
    names = {files.name};
    nums = str2double(regexprep(names, '^kinematics_(\d+)_.*$', '$1'));
    [~, order] = sort(nums, 'ascend');
    files = files(order);

    for i=1:length(files)
        data=readtable([files(i).folder '\' files(i).name],'FileType','text');
        I_pelvis_top=find(contains(data.Properties.VariableNames,'pelvis_top'));
        out(i).pelvis_top=table2array(data(2:end,I_pelvis_top:I_pelvis_top+2));
        I_hip=find(contains(data.Properties.VariableNames,'hip'));
        out(i).hip=table2array(data(2:end,I_hip:I_hip+2));
        I_pelvis_bottom=find(contains(data.Properties.VariableNames,'pelvis_bottom'));
        out(i).pelvis_bottom=table2array(data(2:end,I_pelvis_bottom:I_pelvis_bottom+2));
        I_knee=find(contains(data.Properties.VariableNames,'knee'));
        out(i).knee=table2array(data(2:end,I_knee:I_knee+2));
        I_ankle=find(contains(data.Properties.VariableNames,'ankle'));
        out(i).ankle=table2array(data(2:end,I_ankle:I_ankle+2));
        I_mtp=find(contains(data.Properties.VariableNames,'mtp'));
        out(i).mtp=table2array(data(2:end,I_mtp:I_mtp+2));
        I_toe=find(contains(data.Properties.VariableNames,'toe'));
        out(i).toe=table2array(data(2:end,I_toe:I_toe+2));
    

    end
end

function plot2d(data,Color)
    try
        for i=1:length(data)
            Colorj=Color;
            if length(data)>10
                Colorj(2)=(1/length(data))*(i-1);
            else
                Colorj(2)=0.1*(i-1);
            end
            points=fieldnames(data(i));
            for j=1:length(points)
                pointj=data(i).(points{j});
                subplot(4,6,(j-1)*3+1)
                plot([0:800],pointj(:,1),'Color',Colorj);
                title(['X ' strrep(points{j},'_',' ')]);
                hold all;
                subplot(4,6,(j-1)*3+2)
                plot([0:800],pointj(:,2),'Color',Colorj);
                title('Y');
                hold all;
                subplot(4,6,(j-1)*3+3)
                plot([0:800],pointj(:,3),'Color',Colorj);
                title('Z');
                hold all;
            end
        end
    catch
        keyboard;
    end

end




function plot3d(data)


    for i=1:length(data)
        points=fieldnames(data(i));
        for j=1:length(points)
            for k=1:length(data(i).(points{j}))
                pointj=data(i).(points{j});
                plot3(pointj(k,1),pointj(k,2),pointj(k,3),'o');
                hold all;
            end
        end
    end
end

function [alldata mdata]=ConcatenateAll(data3d_all);

    markers={'pelvis_top','hip','pelvis_bottom','ankle','mtp','toe'};
    for i=1:length(markers);
        alldata.(markers{i})=[];
        for j=1:length(data3d_all)
            alldata.(markers{i})=[alldata.(markers{i}); data3d_all(j).(markers{i})];
        end
        if strcmp(markers{i},'ankle')||strcmp(markers{i},'mtp')||strcmp(markers{i},'toe')
            nframes=length(alldata.(markers{i}));
            nframesxtrial=size(data3d_all(j).(markers{i}),1);
            mdata.(markers{i})=mean(alldata.(markers{i})(1:nframesxtrial:end,:));
        else
            mdata.(markers{i})=mean(alldata.(markers{i}));
        end

        

    end
end

function plot2d_concatenated(alldata)

    markers={'pelvis_top','hip','pelvis_bottom','ankle'};
    for i=1:length(markers)
        subplot(4,3,(i-1)*3+1);
        plot(alldata.(markers{i})(:,1));
        title(['X ' strrep(markers{i},'_',' ')]);
        ylabel('[mm]');
        subplot(4,3,(i-1)*3+2);
        plot(alldata.(markers{i})(:,2));
        title(['Y ' strrep(markers{i},'_',' ')]);
        ylabel('[mm]');
        subplot(4,3,(i-1)*3+3);
        plot(alldata.(markers{i})(:,3));
        title(['Z ' strrep(markers{i},'_',' ')]);
        ylabel('[mm]');
    end

end

function [alldata mdata]=ConcatenateAll_bypert(data3d_baseline,data3d_achillescut,main_folder,pertdata_baseline,pertdata_achillescut)
    markers={'pelvis_top','hip','pelvis_bottom','ankle','mtp','toe'};
    if contains(main_folder,'May2025')
        allperts=[30 20 10 0 0 10 20 30];
    elseif contains(main_folder,'rat22')
        allperts=[30 20 10 0 -10 -20 -20 -10 0 10 20 30];
    else
        allperts=[30 20 10 0 -10 -10 0 10 20 30];
    end
    unique_perts=unique(allperts);
    for i=1:length(markers)
        for pert=1:length(unique_perts)
            alldata(pert).(markers{i})=[];
            alldata(pert).perts=unique_perts(pert);
            I=find(allperts==unique_perts(pert));
            for j=1:length(I);
                alldata(pert).(markers{i})=[alldata(pert).(markers{i}); data3d_baseline(I(j)).(markers{i}); data3d_achillescut(I(j)).(markers{i})];
            end
            if strcmp(markers{i},'ankle')||strcmp(markers{i},'mtp')||strcmp(markers{i},'toe')
                nframesxtrial=length(data3d_baseline(I(j)).(markers{i}));
                mdata(pert).(markers{i})=mean(alldata(pert).(markers{i})(1:nframesxtrial:end,:));

                ncases=size(alldata(pert).(markers{i}),1)/nframesxtrial;
                for k=1:ncases
                    pos_0=alldata(pert).(markers{i})((k-1)*nframesxtrial+1:k*nframesxtrial,:);
                    p_centered=pos_0-mean(pos_0,1);
                    [~,~,V]=svd(p_centered,'econ');
                    vs.(markers{i})(:,k)=V(:,1);
                end
                mdata(pert).(['v_' markers{i}])=mean(vs.(markers{i}),2)./norm(mean(vs.(markers{i}),2));

            else
                mdata(pert).(markers{i})=mean(alldata(pert).(markers{i}));
            end
            mdata(pert).perts=unique_perts(pert);
        end
    end

    %motor data
    for pert=1:length(unique_perts)
        alldatamotor(pert).pertdata=[];
        alldatamotor(pert).perts=unique_perts(pert);
        I=find(allperts==unique_perts(pert));
        for j=1:length(I);
            alldatamotor(pert).pertdata=[alldatamotor(pert).pertdata pertdata_baseline(:,I(j)) pertdata_achillescut(:,I(j))];
        end
        mdata(pert).pert=mean(alldatamotor(pert).pertdata,2);
    end

end


function pertdata=GetPertData(files)
    names = {files.name};
    nums = str2double(regexprep(names, '^kinematics_(\d+)_.*$', '$1'));
    [~, order] = sort(nums, 'ascend');
    files = files(order);
for i=1:length(files)
    filename=strrep(files(i).name,'kinematics','motor');
    filename=strrep(filename,'.trc','.csv');
    folder=files(i).folder;
    folder=strrep(folder,'kinematics','perturbation');
    data_pert=importdata([folder '\' filename]);
    I_motorpert=find(contains(data_pert.colheaders,'filtered_position'));
    pertdata(:,i)=data_pert.data(:,I_motorpert);
end

end

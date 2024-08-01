%LoadDataMay

main_folder=pwd;
cd('DataMay/kinematics')
pkl_files=dir('*.pkl');

for i=1:length(pkl_files)
    
    fid=py.open(pkl_files(i).name,'rb');
    data=py.pickle.load(fid);
    datacell=cell(data);
    sufix=pkl_files(i).name(12:end);
    trial=strrep(sufix,'.pkl','');
    
    if ~isdir(trial)
        mkdir(trial);
    end
    for j=1:length(datacell)
        datacell{j}.T.to_csv([trial '\perturb' num2str(j) '.csv']);
    end
end

cd(main_folder);
cd('DataMay/motion&force')
pkl_files=dir('*.pkl');

for i=1:length(pkl_files)
    
    fid=py.open(pkl_files(i).name,'rb');
    data=py.pickle.load(fid);
    datacell=cell(data);
    trial=strrep(pkl_files(i).name,'_trns.pkl','');
        
    if ~isdir(trial)
        mkdir(trial);
    end
    for j=1:length(datacell)
        datacell{j}.T.to_csv([trial '\perturb' num2str(j) '.csv']);
    end
end

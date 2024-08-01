files=dir('*/**');
for i=1:length(files)
    filepath=[files(i).folder '\' files(i).name];
    if exist(filepath)==2&&strcmp(files(i).name(end-3:end),'.pkl')
        fid=py.open(filepath,'rb');
        data=py.pickle.load(fid);
        testTable = df2t(data);   
        save([files(i).folder '\' strrep(files(i).name,'.pkl','.mat')]);
    end
end


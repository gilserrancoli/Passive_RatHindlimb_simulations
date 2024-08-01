[range_data colheaders]=LoadData();

ndata=5000;
data_all(:,1)=0:0.01:(ndata-1)/100;
data_all(:,8)=3.7;
% data_all(:,9:15)=rand(ndata,1)*diff((range_data(:,8:end)))*1.5-(diff((range_data(:,8:end)))*1.5)/2+mean(range_data(:,8:end));
% just taking ranges from exp data too small variability

ranges_custom=ones(1,7)*30;
ranges_custom([1 4])=60;
for j=8:14
    data_all(:,j+1)=rand(ndata,1).*ranges_custom(:,j-7)-ranges_custom(:,j-7)/2+mean(range_data(:,j));
end

q.data=data_all;
q.labels=colheaders;

write_motionFile(q,'dummy_motion.mot');


function  [range_data colheaders]=LoadData()
cd('..');
range_data=[];
    current_folder=pwd;
    kinfiles=dir('kinematics/forward/*.mot');
    forcefiles=dir('motion_force/forward_mean_response/*.mot');
    for i=1:length(kinfiles);
       kinfilename=kinfiles(i).name; 
       kindata=importdata([kinfiles(i).folder '/' kinfiles(i).name]);
       kindata=ProcessKinematics(kindata);
       if isempty(range_data)
           range_data(1,:)=min(kindata.data(:,2:end));
           range_data(2,:)=max(kindata.data(:,2:end));
       else
            range_data(1,:)=min([range_data(1,:);kindata.data(:,2:end)]); 
            range_data(2,:)=max([range_data(2,:);kindata.data(:,2:end)]); 
       end
       
    end
    colheaders=kindata.colheaders;
    cd('Polynomials');
end
function kindata_out=ProcessKinematics(kindata)
kindata_out=kindata;
v=diff(kindata.data(:,2:end))./diff(kindata.data(:,1));
for i=1:size(v,2)
    
    [out,I]=rmoutliers(v(:,i));
    if any(I)
        aux=kindata.data(~I,i+1);
        time=kindata.data(~I,1);
        kindata_out.data(:,i+1)=interp1(time,aux,kindata.data(:,1));
    end
end



end
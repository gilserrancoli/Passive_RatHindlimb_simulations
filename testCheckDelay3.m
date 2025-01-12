kin=importdata('C:\Gil\Collaborations\MatthewTresch\Zhong_ParameterEstimation\Optimization\DataSeptember\kinematics\data_kin.csv');
pert=readtable('C:\Gil\Collaborations\MatthewTresch\Zhong_ParameterEstimation\Optimization\DataSeptember\perturbation\data.csv','FileType','Text');

perttime=pert.relative_time;
pertpos=pert.position;
pertpos=pert.filtered_position;

kintime=kin.data(:,2);
kintoe=kin.data(:,56);
rangekintoe=max(kintoe)-min(kintoe);
midkintoe=max(kintoe)-rangekintoe/2;
plot(kintime,(kintoe-midkintoe)/rangekintoe);


hold all;
rangepertpos=max(pertpos)-min(pertpos);
midpertpos=max(pertpos)-rangepertpos/2;
plot(perttime,(pertpos-midpertpos)/rangepertpos);
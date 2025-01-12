pert=importdata('DataJuly\perturbation\config0.csv');
kin=readtable('DataJuly\kinematics\config0_kin.trc','FileType','Text');

perttime=pert.data(:,1);
pertpos=pert.data(:,2);

kintime=kin.Time(2:end);
kintoe=kin.toe(2:end);
rangekintoe=max(kintoe)-min(kintoe);
midkintoe=max(kintoe)-rangekintoe/2;
plot(kintime,(kintoe-midkintoe)/rangekintoe);


hold all;
rangepertpos=max(pertpos)-min(pertpos);
midpertpos=max(pertpos)-rangepertpos/2;
plot(perttime,(pertpos-midpertpos)/rangepertpos);
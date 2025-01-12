force=importdata('DataJuly\perturbation\config0.csv');

kin=readtable('DataJuly\kinematics\config0_kin.trc','FileType','text');

meanf=mean(force.data(:,2));
range=max(force.data(:,2))-min(force.data(:,2));
plot(force.data(:,1),(force.data(:,2)-meanf)/(range/2));
hold all;
meank=mean(table2array(kin(:,15)),'omitnan');
rangek=max(table2array(kin(:,15)))-min(table2array(kin(:,15)));
plot(table2array(kin(:,2)),(table2array(kin(:,15))-meank)/(rangek/2));
legend('force x normalized','x data normalized')



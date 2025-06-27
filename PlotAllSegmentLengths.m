kinX_ac{1}=importdata("DataMarch2025\achillescut_ForwardOnly\kinematics\kinematics_0_pos30.0.csv");
kinX_ac{2}=importdata("DataMarch2025\achillescut_ForwardOnly\kinematics\kinematics_1_pos15.0.csv");
kinX_ac{3}=importdata("DataMarch2025\achillescut_ForwardOnly\kinematics\kinematics_2_pos0.0.csv");
kinX_ac{4}=importdata("DataMarch2025\achillescut_ForwardOnly\kinematics\kinematics_3_pos-15.0.csv");
kinX_ac{5}=importdata("DataMarch2025\achillescut_ForwardOnly\kinematics\kinematics_4_pos0.0.csv");
kinX_ac{6}=importdata("DataMarch2025\achillescut_ForwardOnly\kinematics\kinematics_5_pos15.0.csv");
kinX_ac{7}=importdata("DataMarch2025\achillescut_ForwardOnly\kinematics\kinematics_6_pos30.0.csv");

kinX{1}=importdata("DataMarch2025\baseline_ForwardOnly\kinematics\kinematics_0_pos30.0.csv");
kinX{2}=importdata("DataMarch2025\baseline_ForwardOnly\kinematics\kinematics_1_pos15.0.csv");
kinX{3}=importdata("DataMarch2025\baseline_ForwardOnly\kinematics\kinematics_2_pos0.0.csv");
kinX{4}=importdata("DataMarch2025\baseline_ForwardOnly\kinematics\kinematics_3_pos-15.0.csv");
kinX{5}=importdata("DataMarch2025\baseline_ForwardOnly\kinematics\kinematics_4_pos0.0.csv");
kinX{6}=importdata("DataMarch2025\baseline_ForwardOnly\kinematics\kinematics_5_pos15.0.csv");
kinX{7}=importdata("DataMarch2025\baseline_ForwardOnly\kinematics\kinematics_6_pos30.0.csv");

figure
for i=1:7
    c=[1 0+0.1*i 0+0.1*i];
    PlotLengths(kinX_ac{i},c);
    hold all;
end

figure
for i=1:7
    c=[1 0+0.1*i 0+0.1*i];
    PlotLengths(kinX{i},c);
    hold all;
end


function PlotLengths(kin,c);
time=kin.data(:,1);
pelvis_topI=find(contains(kin.colheaders,'filtered_smoothed_pelvis top'));
pelvis_bottomI=find(contains(kin.colheaders,'filtered_smoothed_pelvis bottom'));
hipI=find(contains(kin.colheaders,'filtered_smoothed_hip'));
kneeI=find(contains(kin.colheaders,'filtered_smoothed_knee'));
ankleI=find(contains(kin.colheaders,'filtered_smoothed_ankle'));
mtpI=find(contains(kin.colheaders,'filtered_smoothed_MTP'));
toeI=find(contains(kin.colheaders,'filtered_smoothed_toe'));

femur_l=sqrt(sum((kin.data(:,hipI)-kin.data(:,kneeI)).^2,2));
tibia_l=sqrt(sum((kin.data(:,kneeI)-kin.data(:,ankleI)).^2,2));
foot_l=sqrt(sum((kin.data(:,ankleI)-kin.data(:,toeI)).^2,2));

subplot(1,3,1);
plot(time,femur_l,'Color',c);
hold all;
xlabel('time [s]');
ylabel('length [mm]');
title('femur length');

subplot(1,3,2);
plot(time,tibia_l,'Color',c);
hold all;
xlabel('time [s]');
ylabel('length [mm]');
title('tibia length');

subplot(1,3,3);
plot(time,foot_l,'Color',c);
hold all;
xlabel('time [s]');
ylabel('length [mm]');
title('foot length');

end
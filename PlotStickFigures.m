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
    c=[0 0+0.1*i 0+0.1*i];
    PlotStickFigure(kinX_ac{i},c);
    hold all;
end
xlim([-160 -60]);
ylim([-30 70]);

figure
for i=1:7
    c=[1 0+0.1*i 0+0.1*i];
    PlotStickFigure(kinX{i},c);
    hold all;
end
xlim([-160 -60]);
ylim([-30 70]);


function PlotStickFigure(kin,c)

pelvis_topI=find(contains(kin.colheaders,'filtered_smoothed_pelvis top'));
pelvis_bottomI=find(contains(kin.colheaders,'filtered_smoothed_pelvis bottom'));
hipI=find(contains(kin.colheaders,'filtered_smoothed_hip'));
kneeI=find(contains(kin.colheaders,'filtered_smoothed_knee'));
ankleI=find(contains(kin.colheaders,'filtered_smoothed_ankle'));
mtpI=find(contains(kin.colheaders,'filtered_smoothed_MTP'));
toeI=find(contains(kin.colheaders,'filtered_smoothed_toe'));

plot([kin.data(1,pelvis_topI(1)) kin.data(1,pelvis_bottomI(1))],...
    [kin.data(1,pelvis_topI(3)) kin.data(1,pelvis_bottomI(3))],'Color',c);
hold all;
plot([kin.data(1,hipI(1)) kin.data(1,kneeI(1)) kin.data(1,ankleI(1)) kin.data(1,mtpI(1)) kin.data(1,toeI(1))],...
    [kin.data(1,hipI(3)) kin.data(1,kneeI(3)) kin.data(1,ankleI(3)) kin.data(1,mtpI(3)) kin.data(1,toeI(3))],'Color',c);


end
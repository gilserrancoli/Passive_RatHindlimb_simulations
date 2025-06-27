kinX_ac_mot{1}=importdata("DataMarch2025\achillescut_ForwardOnly\perturbation\motor_0_pos30.0.mot");
kinX_ac_mot{2}=importdata("DataMarch2025\achillescut_ForwardOnly\perturbation\motor_1_pos15.0.mot");
kinX_ac_mot{3}=importdata("DataMarch2025\achillescut_ForwardOnly\perturbation\motor_2_pos0.0.mot");
kinX_ac_mot{4}=importdata("DataMarch2025\achillescut_ForwardOnly\perturbation\motor_3_pos-15.0.mot");
kinX_ac_mot{5}=importdata("DataMarch2025\achillescut_ForwardOnly\perturbation\motor_4_pos0.0.mot");
kinX_ac_mot{6}=importdata("DataMarch2025\achillescut_ForwardOnly\perturbation\motor_5_pos15.0.mot");
kinX_ac_mot{7}=importdata("DataMarch2025\achillescut_ForwardOnly\perturbation\motor_6_pos30.0.mot");

kinX_mot{1}=importdata("DataMarch2025\baseline_ForwardOnly\perturbation\motor_0_pos30.0.mot");
kinX_mot{2}=importdata("DataMarch2025\baseline_ForwardOnly\perturbation\motor_1_pos15.0.mot");
kinX_mot{3}=importdata("DataMarch2025\baseline_ForwardOnly\perturbation\motor_2_pos0.0.mot");
kinX_mot{4}=importdata("DataMarch2025\baseline_ForwardOnly\perturbation\motor_3_pos-15.0.mot");
kinX_mot{5}=importdata("DataMarch2025\baseline_ForwardOnly\perturbation\motor_4_pos0.0.mot");
kinX_mot{6}=importdata("DataMarch2025\baseline_ForwardOnly\perturbation\motor_5_pos15.0.mot");
kinX_mot{7}=importdata("DataMarch2025\baseline_ForwardOnly\perturbation\motor_6_pos30.0.mot");


legend_text={'30','15','0','-15','0','15','30'};
%% with achiles tendon cut

for i=1:7
    c_mot=[0 0+0.1*i 0+0.1*i];
    Plot_comparison(kinX_ac_mot{i},c_mot);

end
legend_text={'30','15','0','-15','0','15','30'};
legend(legend_text);

%% baseline

for i=1:7
    c_mot=[0 0+0.1*i 0+0.1*i];
    Plot_comparison(kinX_ac_mot{i},c_mot);

end
legend(legend_text);


function Plot_comparison(kinX_ac_mot,c)
time=kinX_ac_mot.data(:,1);
FxI=find(contains(kinX_ac_mot.colheaders,'FT_vx'));
FyI=find(contains(kinX_ac_mot.colheaders,'FT_vy'));
FzI=find(contains(kinX_ac_mot.colheaders,'FT_vz'));
PxI=find(contains(kinX_ac_mot.colheaders,'FT_px'));
PyI=find(contains(kinX_ac_mot.colheaders,'FT_py'));
PzI=find(contains(kinX_ac_mot.colheaders,'FT_pz'));
TxI=find(contains(kinX_ac_mot.colheaders,'FT_tx'));
TyI=find(contains(kinX_ac_mot.colheaders,'FT_ty'));
TzI=find(contains(kinX_ac_mot.colheaders,'FT_tz'));

subplot(3,3,1);
plot(time,kinX_ac_mot.data(:,FxI),'Color',c);
hold all;
title('Fx');
ylabel('Force [N]');
subplot(3,3,2);
plot(time,kinX_ac_mot.data(:,FyI),'Color',c);
hold all;
title('Fy');
ylabel('Force [N]');
subplot(3,3,3);
plot(time,kinX_ac_mot.data(:,FzI),'Color',c);
hold all;
title('Fz');
ylabel('Force [N]');

subplot(3,3,4);
plot(time,kinX_ac_mot.data(:,PxI),'Color',c);
hold all;
title('Px');
ylabel('pos [m]');
subplot(3,3,5);
plot(time,kinX_ac_mot.data(:,PyI),'Color',c);
hold all;
title('Py');
ylabel('pos [m]');
subplot(3,3,6);
plot(time,kinX_ac_mot.data(:,PzI),'Color',c);
hold all;
title('Pz');
ylabel('pos [m]');

subplot(3,3,7);
plot(time,kinX_ac_mot.data(:,TxI),'Color',c);
hold all;
title('Tx');
ylabel('Mon [Nm]');
xlabel('time [s]');
subplot(3,3,8);
plot(time,kinX_ac_mot.data(:,TyI),'Color',c);
hold all;
title('Ty');
ylabel('Mon [Nm]');
xlabel('time [s]');
subplot(3,3,9);
plot(time,kinX_ac_mot.data(:,TzI),'Color',c);
hold all;
title('Tz');
ylabel('Mon [Nm]');
xlabel('time [s]');



end
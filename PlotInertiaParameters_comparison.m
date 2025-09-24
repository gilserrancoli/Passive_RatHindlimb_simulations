S1=load('sol_val_optJointPassive_DataMay2025baselineF_tolem5_Jminres777_JpenKDT1e-3_JminKDT1e-5_JminInertiaP1e-2_nokneemarker2D_pert2_wKDinteractionsnoanklehip_ParamID.mat');
S2=load('sol_val_optJointPassive_DataMay2025achillescutF_tolem5_Jminres777_JpenKDT1e-3_JminKDT1e-5_JminInertiaP1e-2_nokneemarker2D_pert2_wKDinteractionsnoanklehip_ParamID.mat');

bounds=S1.sol_val.bounds;
scaling=S1.sol_val.scaling;

%% Plot inertia parameters
c=[0 0 0; 0 1 0; 1 0 0; 0 0 0];
subplot(3,3,1);
h=bar([bounds.inertiaParam.lower(1)' S1.sol_val.inertiaParam(1) S2.sol_val.inertiaParam(1) bounds.inertiaParam.upper(1)'].*scaling.inertiaparam(1));
h.FaceColor='flat';
h.CData(1,:)=c(1,:);
h.CData(2,:)=c(2,:);
h.CData(3,:)=c(3,:);
h.CData(4,:)=c(4,:);
set(gca,'XTickLabels',{'LB','m','UB'});
title('mass femur');
ylabel('mass [kg]')
subplot(3,3,2);
h=bar([bounds.inertiaParam.lower(4)' S1.sol_val.inertiaParam(4) S2.sol_val.inertiaParam(4) bounds.inertiaParam.upper(4)'].*scaling.inertiaparam(4));
title('mass tibia');
ylabel('mass [kg]')
h.FaceColor='flat';
h.CData(1,:)=c(1,:);
h.CData(2,:)=c(2,:);
h.CData(3,:)=c(3,:);
h.CData(4,:)=c(4,:);
set(gca,'XTickLabels',{'LB','m','UB'});
subplot(3,3,3);
h=bar([bounds.inertiaParam.lower(5)' S1.sol_val.inertiaParam(7) S2.sol_val.inertiaParam(7) bounds.inertiaParam.upper(7)'].*scaling.inertiaparam(7));
title('mass foot');
ylabel('mass [kg]')
h.FaceColor='flat';
h.CData(1,:)=c(1,:);
h.CData(2,:)=c(2,:);
h.CData(3,:)=c(3,:);
h.CData(4,:)=c(4,:);
set(gca,'XTickLabels',{'LB','m','UB'});

subplot(3,3,4);
h=bar([bounds.inertiaParam.lower(2)' S1.sol_val.inertiaParam(2) S2.sol_val.inertiaParam(2) bounds.inertiaParam.upper(2)'].*scaling.inertiaparam(2));
h.FaceColor='flat';
h.CData(1,:)=c(1,:);
h.CData(2,:)=c(2,:);
h.CData(3,:)=c(3,:);
h.CData(4,:)=c(4,:);
set(gca,'XTickLabels',{'LB','I_z','UB'});
title('I_z femur');
ylabel('I_z [kg m^2]')
subplot(3,3,5);
h=bar([bounds.inertiaParam.lower(5)' S1.sol_val.inertiaParam(5) S2.sol_val.inertiaParam(5) bounds.inertiaParam.upper(5)'].*scaling.inertiaparam(5));
h.FaceColor='flat';
h.CData(1,:)=c(1,:);
h.CData(2,:)=c(2,:);
h.CData(3,:)=c(3,:);
h.CData(4,:)=c(4,:);
set(gca,'XTickLabels',{'LB','I_z','UB'});
title('I_z tibia');
ylabel('I_z [kg m^2]')
subplot(3,3,6);
h=bar([bounds.inertiaParam.lower(8)' S1.sol_val.inertiaParam(8) S2.sol_val.inertiaParam(8) bounds.inertiaParam.upper(8)'].*scaling.inertiaparam(8));
h.FaceColor='flat';
h.CData(1,:)=c(1,:);
h.CData(2,:)=c(2,:);
h.CData(3,:)=c(3,:);
h.CData(4,:)=c(4,:);
set(gca,'XTickLabels',{'LB','I_z','UB'});
title('I_z foot');
ylabel('I_z [kg m^2]')

subplot(3,3,7);
h=bar([bounds.inertiaParam.lower(3)' S1.sol_val.inertiaParam(3) S2.sol_val.inertiaParam(3) bounds.inertiaParam.upper(3)'].*scaling.inertiaparam(3));
h.FaceColor='flat';
h.CData(1,:)=c(1,:);
h.CData(2,:)=c(2,:);
h.CData(3,:)=c(3,:);
h.CData(4,:)=c(4,:);
set(gca,'XTickLabels',{'LB','y','UB'});
title('y_{COM} femur');
ylabel('y_{COM} [m]')
subplot(3,3,8);
h=bar([bounds.inertiaParam.lower(6)' S1.sol_val.inertiaParam(6) S2.sol_val.inertiaParam(6) bounds.inertiaParam.upper(6)'].*scaling.inertiaparam(6));
h.FaceColor='flat';
h.CData(1,:)=c(1,:);
h.CData(2,:)=c(2,:);
h.CData(3,:)=c(3,:);
h.CData(4,:)=c(4,:);
set(gca,'XTickLabels',{'LB','y','UB'});
title('y_{COM} tibia');
ylabel('y_{COM} [m]')
subplot(3,3,9);
h=bar([bounds.inertiaParam.lower(9)' S1.sol_val.inertiaParam(9) S2.sol_val.inertiaParam(9) bounds.inertiaParam.upper(9)'].*scaling.inertiaparam(9));
h.FaceColor='flat';
h.CData(1,:)=c(1,:);
h.CData(2,:)=c(2,:);
h.CData(3,:)=c(3,:);
h.CData(4,:)=c(4,:);
set(gca,'XTickLabels',{'LB','y','UB'});
title('y_{COM} foot');
ylabel('y_{COM} [m]')

hold on;

dummyBars = gobjects(4, 1); % Preallocate
for i = 1:4
    dummyBars(i) = bar(nan, nan, 'FaceColor', c(i,:));
end
hold off;
legend(dummyBars,{'LB','intact','achilles cut','UB'})





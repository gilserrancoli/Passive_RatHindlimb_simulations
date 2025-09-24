S1=load('sol_val_optJointPassive_DataAugust2025baselineF_tolem5_Jminres777_JpenKDT1e-3_JminKDT1e-5_JminInertiaP1e-2_nokneemarker2D_pert2_ParamID_openbounds4KD.mat');
S2=load('sol_val_optJointPassive_DataAugust2025achillescutF_tolem5_Jminres777_JpenKDT1e-3_JminKDT1e-5_JminInertiaP1e-2_nokneemarker2D_pert2_ParamID_openbounds4KD.mat');


for i=1:length(S1.sol_val.QsQdots_col_unsc)
    figure(i);
    set(gcf,'Position',[360,97.7,745,237.3]);
    subplot(1,3,1)
    plot(S1.sol_val.tgrid{i}(S1.sol_val.t2plot{i}),S1.sol_val.QsQdots_col_unsc{i}(:,1)*180/pi);
    hold all;
    plot(S2.sol_val.tgrid{i}(S2.sol_val.t2plot{i}),S2.sol_val.QsQdots_col_unsc{i}(:,1)*180/pi);
    title('hip flexion');
    ylabel('[º]');
    xlabel('time [s]');

    subplot(1,3,2)
    plot(S1.sol_val.tgrid{i}(S1.sol_val.t2plot{i}),S1.sol_val.QsQdots_col_unsc{i}(:,7)*180/pi);
    hold all;
    plot(S2.sol_val.tgrid{i}(S2.sol_val.t2plot{i}),S2.sol_val.QsQdots_col_unsc{i}(:,7)*180/pi);
    title('knee flexion');
    ylabel('[º]');
    xlabel('time [s]');

    subplot(1,3,3)
    plot(S1.sol_val.tgrid{i}(S1.sol_val.t2plot{i}),S1.sol_val.QsQdots_col_unsc{i}(:,9)*180/pi);
    hold all;
    plot(S2.sol_val.tgrid{i}(S2.sol_val.t2plot{i}),S2.sol_val.QsQdots_col_unsc{i}(:,9)*180/pi);
    title('ankle flexion');
    ylabel('[º]');
    xlabel('time [s]');
end

for i=1:length(S1.sol_val.forces_prescribed)
    figure(i);
    set(gcf,'Position',[360,97.7,745,237.3]);
    subplot(1,3,1)
    plot(S1.sol_val.tgrid{i}(S1.sol_val.t2plot{i}),S1.sol_val.forces_prescribed{i}(:,1));
    hold all;
    plot(S2.sol_val.tgrid{i}(S2.sol_val.t2plot{i}),S2.sol_val.forces_prescribed{i}(:,1));
    title('F_x');
    ylabel('[N]');
    xlabel('time [s]');

    subplot(1,3,2)
    plot(S1.sol_val.tgrid{i}(S1.sol_val.t2plot{i}),S1.sol_val.forces_prescribed{i}(:,2));
    hold all;
    plot(S2.sol_val.tgrid{i}(S2.sol_val.t2plot{i}),S2.sol_val.forces_prescribed{i}(:,2));
    title('F_y');
    ylabel('[N]');
    xlabel('time [s]');

    subplot(1,3,3)
    plot(S1.sol_val.tgrid{i}(S1.sol_val.t2plot{i}),S1.sol_val.forces_prescribed{i}(:,3));
    hold all;
    plot(S2.sol_val.tgrid{i}(S2.sol_val.t2plot{i}),S2.sol_val.forces_prescribed{i}(:,3));
    title('F_z');
    ylabel('[N]');
    xlabel('time [s]');
end
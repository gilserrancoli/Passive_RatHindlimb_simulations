S1=load('sol_val_optKinematics_optJointPassive_optInertia_DataMarch2025IntactForward&Backward_tolem5_Jqtrack100_Jqtrackdot10_Jtrackqd2dot001_JpenKDT1e-3_JminKDT1e-5_JminInertiaP1e-2_N40_WRONG_noCONTCONST.mat');
S2=load('sol_val_optKinematics_optJointPassive_optInertia_DataMarch2025achillescutForward&Backward_tolem5_Jqtrack100_Jqtrackdot10_Jtrackqd2dot001_JpenKDT1e-3_JminKDT1e-5_JminInertiaP1e-2_N40_WRONG_noCONTCONST.mat');


for i=1:length(S1.sol_val.out_opt)
    figure(i)
    set(gcf,'Position',[360,97.7,749.7,216]);
    subplot(1,3,1)
    plot(S1.sol_val.tgrid_col{i},S1.sol_val.out_opt{i}(:,8));
    hold all;
    plot(S2.sol_val.tgrid_col{i},S2.sol_val.out_opt{i}(:,8));
    title('hip flexion');
    ylabel('M [Nm]')
    xlabel('time [s]')
    subplot(1,3,2)
    plot(S1.sol_val.tgrid_col{i},S1.sol_val.out_opt{i}(:,11));
    hold all;
    plot(S2.sol_val.tgrid_col{i},S2.sol_val.out_opt{i}(:,11));
    title('knee flexion')
    ylabel('M [Nm]')
    xlabel('time [s]')

    subplot(1,3,3)
    plot(S1.sol_val.tgrid_col{i},S1.sol_val.out_opt{i}(:,12));
    hold all;
    plot(S2.sol_val.tgrid_col{i},S2.sol_val.out_opt{i}(:,12));
    title('ankle flexion')
    ylabel('M [Nm]')
    xlabel('time [s]')
    legend({'intact','achilles cut'},'box','off','Position',[0.843,0.1591,0.14,0.15])
end

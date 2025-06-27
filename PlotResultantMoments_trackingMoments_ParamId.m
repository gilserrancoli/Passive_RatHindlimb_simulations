S1=load('sol_val_optJointPassive_optInertia_DataMay2025IntactFB_JpenKDT1e-3_JminKDT1e-5_JminInertiaP1e-2_optInertiaPassprop_ParamId.mat');
S2=load('sol_val_optJointPassive_optInertia_DataMay2025achillescutFB_JpenKDT1e-3_JminKDT1e-5_JminInertiaP1e-2_optInertiaPassprop_ParamId.mat');


for i=1:length(S1.sol_val.out)
    figure(i)
    set(gcf,'Position',[360,97.7,749.7,216]);
    subplot(1,3,1)
    plot(S1.sol_val.tgrid{i}(S1.sol_val.t2plot{i}),S1.sol_val.out{i}(:,8),'b');
    hold all;
    plot(S1.sol_val.tgrid{i}(S1.sol_val.t2plot{i}),S1.sol_val.out_model{i}(:,8),'--b');
    plot(S2.sol_val.tgrid{i}(S1.sol_val.t2plot{i}),S2.sol_val.out{i}(:,8),'r');
    plot(S1.sol_val.tgrid{i}(S1.sol_val.t2plot{i}),S2.sol_val.out_model{i}(:,8),'--r');
    title('hip flexion');
    ylabel('M [Nm]')
    xlabel('time [s]')

    subplot(1,3,2)
    plot(S1.sol_val.tgrid{i}(S1.sol_val.t2plot{i}),S1.sol_val.out{i}(:,11),'b');
    hold all;
    plot(S1.sol_val.tgrid{i}(S1.sol_val.t2plot{i}),S1.sol_val.out_model{i}(:,11),'--b');
    plot(S2.sol_val.tgrid{i}(S1.sol_val.t2plot{i}),S2.sol_val.out{i}(:,11),'r');
    plot(S1.sol_val.tgrid{i}(S1.sol_val.t2plot{i}),S2.sol_val.out_model{i}(:,11),'--r');    
    title('knee flexion')
    ylabel('M [Nm]')
    xlabel('time [s]')

    subplot(1,3,3)
    plot(S1.sol_val.tgrid{i}(S1.sol_val.t2plot{i}),S1.sol_val.out{i}(:,12),'b');
    hold all;
    plot(S1.sol_val.tgrid{i}(S1.sol_val.t2plot{i}),S1.sol_val.out_model{i}(:,12),'--b');
    plot(S2.sol_val.tgrid{i}(S1.sol_val.t2plot{i}),S2.sol_val.out{i}(:,12),'r');
    plot(S1.sol_val.tgrid{i}(S1.sol_val.t2plot{i}),S2.sol_val.out_model{i}(:,12),'--r');    
    title('ankle flexion')
    ylabel('M [Nm]')
    xlabel('time [s]')
    legend({'intact exp','intact model','achilles cut exp','achilles cut model'},'box','off','Position',[0.843,0.1591,0.14,0.15])
end

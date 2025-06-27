S1=load('sol_val_optKinematics_optJointPassive_optInertia_DataMarch2025IntactForward&Backward_tolem5_Jqtrack100_Jqtrackdot01_Jtrackqd2dot1em4_JpenKDT1e-5_JminInertiaP1e-2_N40.mat');
S2=load('sol_val_optKinematics_optJointPassive_optInertia_DataMarch2025achillescutForward&Backward_tolem5_Jqtrack100_Jqtrackdot01_Jtrackqd2dot1em4_JpenKDT1e-5_JminInertiaP1e-2_N40.mat');

x=[30 15 0 -15 0 15 30];

for i=1:length(S1.sol_val.QsQdots)
    t=S1.sol_val.tgrid{i};
    if i>7
        k=i-7;
    else
        k=i;
    end
    if k<=4
        kk=k;
    elseif k==5
        kk=3;
    elseif k==6
        kk=2;
    elseif k==7
        kk=1;
    end

    pert=x(k);

    if i>7
        direct='backward';
    else
        direct='forward';
    end
    titleplot=[direct '-' num2str(pert)];

    figure(i);
    set(gcf,'Position',[71,241.7,829.3,298])
    subplot(1,3,1);
    green=[0.47 0.67 0.19];
    plot(t(S1.sol_val.t2plot{i}),S1.sol_val.QsQdots_col_unsc{i}(1,:)*180/pi,'Color',green);
    hold all;
    plot(t(S1.sol_val.t2plot{i}),S1.sol_val.guess.QsQdots{i}(S1.sol_val.t2plot{i},1)*S1.sol_val.scaling.QsQdots(1)*180/pi,'--','Color',[0.47 0.67 0.19]);
    plot(t(S1.sol_val.t2plot{i}),S2.sol_val.QsQdots_col_unsc{i}(1,:)*180/pi,'Color','r');
    plot(t(S1.sol_val.t2plot{i}),S2.sol_val.guess.QsQdots{i}(S2.sol_val.t2plot{i},1)*S2.sol_val.scaling.QsQdots(1)*180/pi,'--','Color','r');
    title('hip angle');
    xlabel('time [s]');
    ylabel('angle [º]');

    subplot(1,3,2);
    green=[0.47 0.67 0.19];
    plot(t(S1.sol_val.t2plot{i}),S1.sol_val.QsQdots_col_unsc{i}(7,:)*180/pi,'Color',green);
    hold all;
    plot(t(S1.sol_val.t2plot{i}),S1.sol_val.guess.QsQdots{i}(S1.sol_val.t2plot{i},7)*S1.sol_val.scaling.QsQdots(1)*180/pi,'--','Color',[0.47 0.67 0.19]);
    plot(t(S1.sol_val.t2plot{i}),S2.sol_val.QsQdots_col_unsc{i}(7,:)*180/pi,'Color','r');
    plot(t(S1.sol_val.t2plot{i}),S2.sol_val.guess.QsQdots{i}(S2.sol_val.t2plot{i},7)*S2.sol_val.scaling.QsQdots(1)*180/pi,'--','Color','r');
    title('knee angle');
    xlabel('time [s]');
    ylabel('angle [º]');

    subplot(1,3,3);
    green=[0.47 0.67 0.19];
    plot(t(S1.sol_val.t2plot{i}),S1.sol_val.QsQdots_col_unsc{i}(9,:)*180/pi,'Color',green);
    hold all;
    plot(t(S1.sol_val.t2plot{i}),S1.sol_val.guess.QsQdots{i}(S1.sol_val.t2plot{i},9)*S1.sol_val.scaling.QsQdots(1)*180/pi,'--','Color',[0.47 0.67 0.19]);
    plot(t(S1.sol_val.t2plot{i}),S2.sol_val.QsQdots_col_unsc{i}(9,:)*180/pi,'Color','r');
    plot(t(S1.sol_val.t2plot{i}),S2.sol_val.guess.QsQdots{i}(S2.sol_val.t2plot{i},9)*S2.sol_val.scaling.QsQdots(1)*180/pi,'--','Color','r');
    title('ankle angle');
    xlabel('time [s]');
    ylabel('angle [º]');

    legend({'intact model','intact exp','achilles cut model','achilles cut exp'},'Position',[0.82,0.7176,0.165,0.206]);
    annotation('textbox',[0 1 0 0],'String',titleplot, 'FitBoxToText', 'on','EdgeColor', 'none');
end



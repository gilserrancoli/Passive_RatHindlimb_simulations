S1=load('sol_val_optKinematics_optJointPassive_optInertia_DataMay2025IntactForward&Backward_tolem5_Jqtrack100_Jqtrackdot10_Jtrackqd2dot001_JpenKDT1e-3_JminKDT1e-5_JminInertiaP1e-2_N40.mat');
S2=load('sol_val_optKinematics_optJointPassive_optInertia_DataMay2025achillescutForward&Backward_tolem5_Jqtrack100_Jqtrackdot10_Jtrackqd2dot001_JpenKDT1e-3_JminKDT1e-5_JminInertiaP1e-2_N40.mat');

% x=[30 15 0 -15 0 15 30];
k=4; %trial k
figure
set(gcf,'Position',[47,55.7,1174,580]);

for i=1:3
    subplot(2,6,i);
    bar(1,S1.sol_val.Kstiff_unsc{k}(i));
    hold all;
    bar(2,S2.sol_val.Kstiff_unsc{k}(i),'r');

    set(gca,'XTick',[1 2],'XTickLabel',{'Intact','Achilles Cut'});
    if i==1
        title('K hip');
    elseif i==2
        title('K knee');
    elseif i==3
        title('K ankle');
    end
    ylabel('K [Nm/rad]');
end

for i=1:3
    subplot(2,6,i+3);
    bar(1,S1.sol_val.Ddamp_unsc{k}(i));
    hold all;
    bar(2,S2.sol_val.Ddamp_unsc{k}(i),'r');

    set(gca,'XTick',[1 2],'XTickLabel',{'Intact','Achilles Cut'});
    if i==1
        title('D hip');
    elseif i==2
        title('D knee');
    elseif i==3
        title('D ankle');
    end
    ylabel('D [Nm/(rad/s)]');
end



%hip moment
blue=[0,0.45,0.74];
subplot(2,3,4);
plot(S1.sol_val.tgrid_col{k},S1.sol_val.out_opt{k}(:,8),'LineWidth',2,'Color',blue);
hold all;
plot(S1.sol_val.tgrid_col{k},-S1.sol_val.Kstiff_unsc{k}(1)*S1.sol_val.QsQdots_col_unsc{k}(1,:),'--','Color',blue);
plot(S1.sol_val.tgrid_col{k},repmat(-S1.sol_val.Kstiff_unsc{k}(1)*(-S1.sol_val.theta0_unsc{k}(1)),1,length(S1.sol_val.tgrid_col{k})),'.','Color',blue);
plot(S1.sol_val.tgrid_col{k},-S1.sol_val.Ddamp_unsc{k}(1)*S1.sol_val.QsQdots_col_unsc{k}(2,:),'Color',blue);

plot(S1.sol_val.tgrid_col{k},S2.sol_val.out_opt{k}(:,8),'LineWidth',2,'Color','r');
hold all;
plot(S1.sol_val.tgrid_col{k},-S2.sol_val.Kstiff_unsc{k}(1)*S2.sol_val.QsQdots_col_unsc{k}(1,:),'--','Color','r');
plot(S1.sol_val.tgrid_col{k},repmat(-S2.sol_val.Kstiff_unsc{k}(1)*(-S2.sol_val.theta0_unsc{k}(1)),1,length(S2.sol_val.tgrid_col{k})),'.','Color','r');
plot(S1.sol_val.tgrid_col{k},-S2.sol_val.Ddamp_unsc{k}(1)*S2.sol_val.QsQdots_col_unsc{k}(2,:),'-','Color','r');
title('hip moment');
legend({'total intact', '$-K\theta$', '$K\theta_0$', '$-D\dot{\theta}$','total achilles cut','$-K\theta$', '$K\theta_0$', '$-D\dot{\theta}$'}, 'Interpreter', 'latex','Orientation','horizontal','box','off','Position',[0.22,0.00757,0.597,0.0336]);
ylabel('[Nm]');

subplot(2,3,5);
plot(S1.sol_val.tgrid_col{k},S1.sol_val.out_opt{k}(:,11),'LineWidth',2,'Color',blue);
hold all;
plot(S1.sol_val.tgrid_col{k},-S1.sol_val.Kstiff_unsc{k}(2)*S1.sol_val.QsQdots_col_unsc{k}(7,:),'--','Color',blue);
plot(S1.sol_val.tgrid_col{k},repmat(-S1.sol_val.Kstiff_unsc{k}(2)*(-S1.sol_val.theta0_unsc{k}(2)),1,length(S1.sol_val.tgrid_col{k})),'.','Color',blue);
plot(S1.sol_val.tgrid_col{k},-S1.sol_val.Ddamp_unsc{k}(2)*S1.sol_val.QsQdots_col_unsc{k}(8,:),'-','Color',blue);

plot(S1.sol_val.tgrid_col{k},S2.sol_val.out_opt{k}(:,11),'LineWidth',2,'Color','r');
hold all;
plot(S1.sol_val.tgrid_col{k},-S2.sol_val.Kstiff_unsc{k}(2)*S2.sol_val.QsQdots_col_unsc{k}(7,:),'--','Color','r');
plot(S1.sol_val.tgrid_col{k},repmat(-S2.sol_val.Kstiff_unsc{k}(2)*(-S2.sol_val.theta0_unsc{k}(2)),1,length(S1.sol_val.tgrid_col{k})),'.','Color','r');
plot(S1.sol_val.tgrid_col{k},-S2.sol_val.Ddamp_unsc{k}(2)*S2.sol_val.QsQdots_col_unsc{k}(8,:),'-','Color','r');
title('knee moment');

subplot(2,3,6);
plot(S1.sol_val.tgrid_col{k},S1.sol_val.out_opt{k}(:,12),'LineWidth',2,'Color',blue);
hold all;
plot(S1.sol_val.tgrid_col{k},-S1.sol_val.Kstiff_unsc{k}(3)*S1.sol_val.QsQdots_col_unsc{k}(9,:),'--','Color',blue);
plot(S1.sol_val.tgrid_col{k},-S1.sol_val.Kstiff_unsc{k}(3)*(-S1.sol_val.theta0_unsc{k}(3)),'.','Color',blue);
plot(S1.sol_val.tgrid_col{k},-S1.sol_val.Ddamp_unsc{k}(3)*S1.sol_val.QsQdots_col_unsc{k}(10,:),'-','Color',blue);

plot(S1.sol_val.tgrid_col{k},S2.sol_val.out_opt{k}(:,12),'LineWidth',2,'Color','r');
hold all;
plot(S1.sol_val.tgrid_col{k},-S2.sol_val.Kstiff_unsc{k}(3)*S2.sol_val.QsQdots_col_unsc{k}(9,:),'--','Color','r');
plot(S1.sol_val.tgrid_col{k},-S2.sol_val.Kstiff_unsc{k}(3)*(-S2.sol_val.theta0_unsc{k}(3)),'.','Color','r');
plot(S1.sol_val.tgrid_col{k},-S2.sol_val.Ddamp_unsc{k}(3)*S2.sol_val.QsQdots_col_unsc{k}(10,:),'-','Color','r');
title('ankle moment');

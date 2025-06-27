S1=load('sol_val_optKinematics_optJointPassive_optInertia_DataMay2025IntactForward&Backward_tolem5_Jqtrack100_Jqtrackdot10_Jtrackqd2dot001_JpenKDT1e-3_JminKDT1e-5_JminInertiaP1e-2_N40.mat');
S2=load('sol_val_optKinematics_optJointPassive_optInertia_DataMay2025achillescutForward&Backward_tolem5_Jqtrack100_Jqtrackdot10_Jtrackqd2dot001_JpenKDT1e-3_JminKDT1e-5_JminInertiaP1e-2_N40.mat');

scaling=S1.sol_val.scaling;
% x=[30 15 0 -15 0 15 30];
k=1; %trial k
figure
set(gcf,'Position',[47,55.7,1174,580]);

green=[0 0.66 0.33];

hip_mom=-S1.sol_val.Kstiff_unsc{1}(1)*(S1.sol_val.QsQdots_col_unsc{1}(1,:)-S1.sol_val.theta0_unsc{1}(1))-S1.sol_val.Ddamp_unsc{1}(1)*S1.sol_val.QsQdots_col_unsc{1}(2,:);
knee_mom=-S1.sol_val.Kstiff_unsc{1}(2)*(S1.sol_val.QsQdots_col_unsc{1}(7,:)-S1.sol_val.theta0_unsc{1}(2))-S1.sol_val.Ddamp_unsc{1}(2)*S1.sol_val.QsQdots_col_unsc{1}(8,:);
ankle_mom=-S1.sol_val.Kstiff_unsc{1}(3)*(S1.sol_val.QsQdots_col_unsc{1}(9,:)-S1.sol_val.theta0_unsc{1}(3))-S1.sol_val.Ddamp_unsc{1}(3)*S1.sol_val.QsQdots_col_unsc{1}(10,:);

hip_mom_IG=-S1.sol_val.Kstiff_unsc{1}(1)*(S1.sol_val.guess.QsQdots{1}(:,1)*scaling.q-S1.sol_val.theta0_unsc{1}(1))-S1.sol_val.Ddamp_unsc{1}(1)*S1.sol_val.guess.QsQdots{1}(:,2)*scaling.qdot;
knee_mom_IG=-S1.sol_val.Kstiff_unsc{1}(2)*(S1.sol_val.guess.QsQdots{1}(:,7)*scaling.q-S1.sol_val.theta0_unsc{1}(2))-S1.sol_val.Ddamp_unsc{1}(2)*S1.sol_val.guess.QsQdots{1}(:,8)*scaling.qdot;
ankle_mom_IG=-S1.sol_val.Kstiff_unsc{1}(3)*(S1.sol_val.guess.QsQdots{1}(:,9)*scaling.q-S1.sol_val.theta0_unsc{1}(3))-S1.sol_val.Ddamp_unsc{1}(3)*S1.sol_val.guess.QsQdots{1}(:,10)*scaling.qdot;
t2plot=1:161;
t2plot(1:4:end)=[];

%hip moment
subplot(3,2,1);
plot(S1.sol_val.tgrid_col{k},hip_mom,'LineWidth',2,'Color',blue);
hold all;
plot(S1.sol_val.tgrid_col{k},-S1.sol_val.Kstiff_unsc{k}(1)*S1.sol_val.QsQdots_col_unsc{k}(1,:),'--','Color',blue);
plot(S1.sol_val.tgrid_col{k},repmat(-S1.sol_val.Kstiff_unsc{k}(1)*(-S1.sol_val.theta0_unsc{k}(1)),1,length(S1.sol_val.tgrid_col{k})),'.','Color',blue);
plot(S1.sol_val.tgrid_col{k},-S1.sol_val.Ddamp_unsc{k}(1)*S1.sol_val.QsQdots_col_unsc{k}(2,:),'Color',blue);
title('hip moment optimal');
legend({'total intact', '$-K\theta$', '$K\theta_0$', '$-D\dot{\theta}$'}, 'Interpreter', 'latex','Orientation','horizontal','box','off','Position',[0.22,0.00757,0.597,0.0336]);
ylabel('[Nm]');
ylim([-0.03 0.02]);

subplot(3,2,2);
plot(S1.sol_val.tgrid_col{k},hip_mom_IG(t2plot),'LineWidth',2,'Color',green);
hold all;
plot(S1.sol_val.tgrid_col{k},-S1.sol_val.Kstiff_unsc{1}(1)*S1.sol_val.guess.QsQdots{1}((t2plot),1)*scaling.q,'--','Color',green);
plot(S1.sol_val.tgrid_col{k},repmat(-S1.sol_val.Kstiff_unsc{k}(1)*(-S1.sol_val.theta0_unsc{k}(1)),1,length(S1.sol_val.tgrid_col{k})),'.','Color',green);
plot(S1.sol_val.tgrid_col{k},-S1.sol_val.Ddamp_unsc{1}(1)*S1.sol_val.guess.QsQdots{1}((t2plot),2)*scaling.qdot,'Color',green);
ylim([-0.03 0.02]);
title('hip moment with kinematics IG');

%knee moment
subplot(3,2,3);
plot(S1.sol_val.tgrid_col{k},knee_mom,'LineWidth',2,'Color',blue);
hold all;
plot(S1.sol_val.tgrid_col{k},-S1.sol_val.Kstiff_unsc{k}(2)*S1.sol_val.QsQdots_col_unsc{k}(7,:),'--','Color',blue);
plot(S1.sol_val.tgrid_col{k},repmat(-S1.sol_val.Kstiff_unsc{k}(2)*(-S1.sol_val.theta0_unsc{k}(2)),1,length(S1.sol_val.tgrid_col{k})),'.','Color',blue);
plot(S1.sol_val.tgrid_col{k},-S1.sol_val.Ddamp_unsc{k}(2)*S1.sol_val.QsQdots_col_unsc{k}(8,:),'Color',blue);
title('knee moment optimal');
ylabel('[Nm]');
ylim([-0.06 0.04]);

subplot(3,2,4);
plot(S1.sol_val.tgrid_col{k},knee_mom_IG(t2plot),'LineWidth',2,'Color',green);
hold all;
plot(S1.sol_val.tgrid_col{k},-S1.sol_val.Kstiff_unsc{1}(2)*S1.sol_val.guess.QsQdots{1}((t2plot),7)*scaling.q,'--','Color',green);
plot(S1.sol_val.tgrid_col{k},repmat(-S1.sol_val.Kstiff_unsc{k}(2)*(-S1.sol_val.theta0_unsc{k}(2)),1,length(S1.sol_val.tgrid_col{k})),'.','Color',green);
plot(S1.sol_val.tgrid_col{k},-S1.sol_val.Ddamp_unsc{1}(2)*S1.sol_val.guess.QsQdots{1}((t2plot),8)*scaling.qdot,'Color',green);
ylim([-0.06 0.04]);
title('knee moment with kinematics IG');

%ankle moment
subplot(3,2,5);
plot(S1.sol_val.tgrid_col{k},ankle_mom,'LineWidth',2,'Color',blue);
hold all;
plot(S1.sol_val.tgrid_col{k},-S1.sol_val.Kstiff_unsc{k}(3)*S1.sol_val.QsQdots_col_unsc{k}(9,:),'--','Color',blue);
plot(S1.sol_val.tgrid_col{k},repmat(-S1.sol_val.Kstiff_unsc{k}(3)*(-S1.sol_val.theta0_unsc{k}(3)),1,length(S1.sol_val.tgrid_col{k})),'.','Color',blue);
plot(S1.sol_val.tgrid_col{k},-S1.sol_val.Ddamp_unsc{k}(3)*S1.sol_val.QsQdots_col_unsc{k}(10,:),'Color',blue);
title('ankle moment optimal');
ylabel('[Nm]');
ylim([-0.02 0.01]);

subplot(3,2,6);
plot(S1.sol_val.tgrid_col{k},ankle_mom_IG(t2plot),'LineWidth',2,'Color',green);
hold all;
plot(S1.sol_val.tgrid_col{k},-S1.sol_val.Kstiff_unsc{1}(3)*S1.sol_val.guess.QsQdots{1}((t2plot),9)*scaling.q,'--','Color',green);
plot(S1.sol_val.tgrid_col{k},repmat(-S1.sol_val.Kstiff_unsc{k}(3)*(-S1.sol_val.theta0_unsc{k}(3)),1,length(S1.sol_val.tgrid_col{k})),'.','Color',green);
plot(S1.sol_val.tgrid_col{k},-S1.sol_val.Ddamp_unsc{1}(3)*S1.sol_val.guess.QsQdots{1}((t2plot),10)*scaling.qdot,'Color',green);
ylim([-0.02 0.01]);
title('ankle moment with kinematics IG');

% I=1;
% 
% %hip moment
% blue=[0,0.45,0.74];
% subplot(2,3,4);
% plot(S1.sol_val.tgrid_col{k},S1.sol_val.out_opt{k}(:,8),'LineWidth',2,'Color',blue);
% hold all;
% plot(S1.sol_val.tgrid_col{k},-S1.sol_val.Kstiff_unsc{k}(1)*S1.sol_val.QsQdots_col_unsc{k}(1,:),'--','Color',blue);
% plot(S1.sol_val.tgrid_col{k},repmat(-S1.sol_val.Kstiff_unsc{k}(1)*(-S1.sol_val.theta0_unsc{k}(1)),1,length(S1.sol_val.tgrid_col{k})),'.','Color',blue);
% plot(S1.sol_val.tgrid_col{k},-S1.sol_val.Ddamp_unsc{k}(1)*S1.sol_val.QsQdots_col_unsc{k}(2,:),'Color',blue);
% title('hip moment');
% legend({'total intact', '$-K\theta$', '$K\theta_0$', '$-D\dot{\theta}$'}, 'Interpreter', 'latex','Orientation','horizontal','box','off','Position',[0.22,0.00757,0.597,0.0336]);
% ylabel('[Nm]');
% 
% subplot(2,3,5);
% plot(S1.sol_val.tgrid_col{k},S1.sol_val.out_opt{k}(:,11),'LineWidth',2,'Color',blue);
% hold all;
% plot(S1.sol_val.tgrid_col{k},-S1.sol_val.Kstiff_unsc{k}(2)*S1.sol_val.QsQdots_col_unsc{k}(7,:),'--','Color',blue);
% plot(S1.sol_val.tgrid_col{k},repmat(-S1.sol_val.Kstiff_unsc{k}(2)*(-S1.sol_val.theta0_unsc{k}(2)),1,length(S1.sol_val.tgrid_col{k})),'.','Color',blue);
% plot(S1.sol_val.tgrid_col{k},-S1.sol_val.Ddamp_unsc{k}(2)*S1.sol_val.QsQdots_col_unsc{k}(8,:),'-','Color',blue);
% 
% plot(S1.sol_val.tgrid_col{k},S2.sol_val.out_opt{k}(:,11),'LineWidth',2,'Color','r');
% hold all;
% plot(S1.sol_val.tgrid_col{k},-S2.sol_val.Kstiff_unsc{k}(2)*S2.sol_val.QsQdots_col_unsc{k}(7,:),'--','Color','r');
% plot(S1.sol_val.tgrid_col{k},repmat(-S2.sol_val.Kstiff_unsc{k}(2)*(-S2.sol_val.theta0_unsc{k}(2)),1,length(S1.sol_val.tgrid_col{k})),'.','Color','r');
% plot(S1.sol_val.tgrid_col{k},-S2.sol_val.Ddamp_unsc{k}(2)*S2.sol_val.QsQdots_col_unsc{k}(8,:),'-','Color','r');
% title('knee moment');
% 
% subplot(2,3,6);
% plot(S1.sol_val.tgrid_col{k},S1.sol_val.out_opt{k}(:,12),'LineWidth',2,'Color',blue);
% hold all;
% plot(S1.sol_val.tgrid_col{k},-S1.sol_val.Kstiff_unsc{k}(3)*S1.sol_val.QsQdots_col_unsc{k}(9,:),'--','Color',blue);
% plot(S1.sol_val.tgrid_col{k},-S1.sol_val.Kstiff_unsc{k}(3)*(-S1.sol_val.theta0_unsc{k}(3)),'.','Color',blue);
% plot(S1.sol_val.tgrid_col{k},-S1.sol_val.Ddamp_unsc{k}(3)*S1.sol_val.QsQdots_col_unsc{k}(10,:),'-','Color',blue);
% 
% plot(S1.sol_val.tgrid_col{k},S2.sol_val.out_opt{k}(:,12),'LineWidth',2,'Color','r');
% hold all;
% plot(S1.sol_val.tgrid_col{k},-S2.sol_val.Kstiff_unsc{k}(3)*S2.sol_val.QsQdots_col_unsc{k}(9,:),'--','Color','r');
% plot(S1.sol_val.tgrid_col{k},-S2.sol_val.Kstiff_unsc{k}(3)*(-S2.sol_val.theta0_unsc{k}(3)),'.','Color','r');
% plot(S1.sol_val.tgrid_col{k},-S2.sol_val.Ddamp_unsc{k}(3)*S2.sol_val.QsQdots_col_unsc{k}(10,:),'-','Color','r');
% title('ankle moment');

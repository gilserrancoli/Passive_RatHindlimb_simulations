%% Results Tracking moment (sagittal plane)

subplot(3,1,1)
%hip flexion
plot(sol_val.tgrid_col{1},out_opt{i}(:,8),'LineWidth',2);
hold all;
if Options.optimizeMuscleProp
    plot(sol_val.tgrid_col{1},sum(T_hip_flx_opt,2),'Color',[0.85 0.33 0.1])
end
if Options.optimizePassiveJointEl
    plot(sol_val.tgrid_col{1},PassiveM_hip_flx_opt,'Color',[0.93,0.69,0.13]);
end
if Options.optimizeMuscleProp&&Options.optimizePassiveJointEl
    plot(sol_val.tgrid_col{1},sum(T_hip_flx_opt,2)+PassiveM_hip_flx_opt,'LineWidth',2,'Color',[0.49, 0.18,0.56])
end
xlabel('time [s]');
ylabel('Moment [Nm]')
title('Hip')

subplot(3,1,2)
%knee flexion
plot(sol_val.tgrid_col{1},out_opt{i}(:,11),'LineWidth',2)
hold all
if Options.optimizeMuscleProp
    plot(sol_val.tgrid_col{1},sum(T_knee_flx_opt,2),'Color',[0.85 0.33 0.1])
end
if Options.optimizePassiveJointEl
    plot(sol_val.tgrid_col{1},PassiveM_knee_flx_opt,'Color',[0.93,0.69,0.13])
end
% plot(sol_val.res_col_unsc{i}(2,:)')
if Options.optimizeMuscleProp&&Options.optimizePassiveJointEl
    plot(sol_val.tgrid_col{1},sum(T_knee_flx_opt,2)+PassiveM_knee_flx_opt,'LineWidth',2,'Color',[0.49, 0.18,0.56])
end
ylabel('Moment [Nm]')
title('Knee')

subplot(3,1,3)
%ankle flexion
plot(sol_val.tgrid_col{1},out_opt{i}(:,12),'LineWidth',2)
hold all
if Options.optimizeMuscleProp
    plot(sol_val.tgrid_col{1},sum(T_ankle_flx_opt,2),'Color',[0.85 0.33 0.1])
end
if Options.optimizePassiveJointEl
    plot(sol_val.tgrid_col{1},PassiveM_ankle_flx_opt,'Color',[0.93,0.69,0.13])
end
% plot(sol_val.res_col_unsc{i}(2,:)')
if Options.optimizeMuscleProp&&Options.optimizePassiveJointEl
    plot(sol_val.tgrid_col{1},sum(T_ankle_flx_opt,2)+PassiveM_ankle_flx_opt,'LineWidth',2,'Color',[0.49, 0.18,0.56])
end
ylabel('Moment [Nm]')
title('Ankle')
if Options.optimizeMuscleProp&&Options.optimizePassiveJointEl
    legend({'Joint moment ID','Muscle moment','Passive moment','Muscle+Passive moment'},'Orientation','horizontal','Box','off');
elseif Options.optimizeMuscleProp&&~Options.optimizePassiveJointEl
    legend({'Joint moment ID','Muscle moment'},'Orientation','horizontal','Box','off');
elseif ~Options.optimizeMuscleProp&&Options.optimizePassiveJointEl
    legend({'Joint moment ID','Passive moment'},'Orientation','horizontal','Box','off');
end

%% Kinematics sagittal plane
% Angles
subplot(3,1,1)
plot(sol_val.tgrid_col{1},sol_val.QsQdots_col_unsc{i}(:,1*2-1)*180/pi,'LineWidth',2);
xlabel('time [s]');
ylabel('hip angle [º]');

subplot(3,1,2)
plot(sol_val.tgrid_col{1},sol_val.QsQdots_col_unsc{i}(:,4*2-1)*180/pi,'LineWidth',2);
xlabel('time [s]');
ylabel('knee angle [º]');

subplot(3,1,3)
plot(sol_val.tgrid_col{1},sol_val.QsQdots_col_unsc{i}(:,5*2-1)*180/pi,'LineWidth',2);
xlabel('time [s]');
ylabel('ankle angle [º]');

% Velocities
subplot(3,1,1)
plot(sol_val.tgrid_col{1},sol_val.QsQdots_col_unsc{i}(:,1*2),'LineWidth',2);
xlabel('time [s]');
ylabel('hip vel [rad/s]');

subplot(3,1,2)
plot(sol_val.tgrid_col{1},sol_val.QsQdots_col_unsc{i}(:,4*2),'LineWidth',2);
xlabel('time [s]');
ylabel('knee vel [rad/s]');

subplot(3,1,3)
plot(sol_val.tgrid_col{1},sol_val.QsQdots_col_unsc{i}(:,5*2),'LineWidth',2);
xlabel('time [s]');
ylabel('ankle vel [rad/s]');

%% Accelerations
subplot(3,1,1)
plot(sol_val.tgrid_col{1},sol_val.Qd2dot_col_unsc{i}(:,1),'LineWidth',2)
xlabel('time [s]');
ylabel('hip acc [rad/s^2]');

subplot(3,1,2)
plot(sol_val.tgrid_col{1},sol_val.Qd2dot_col_unsc{i}(:,4),'LineWidth',2)
xlabel('time [s]');
ylabel('knee acc [rad/s^2]');

subplot(3,1,3)
plot(sol_val.tgrid_col{1},sol_val.Qd2dot_col_unsc{i}(:,5),'LineWidth',2)
xlabel('time [s]');
ylabel('ankle acc [rad/s^2]');

%% Forces
subplot(3,1,1)
h=plot(sol_val.tgrid_col{1},forces_prescribed{i}(:,1:3),'LineWidth',2)
set(h(1),'Color','r')
set(h(2),'Color','g')
set(h(3),'Color','b')
xlabel('time [s]')
ylabel('force [N]')

subplot(3,1,2)
h=plot(sol_val.tgrid_col{1},forces_prescribed{i}(:,7:9),'LineWidth',2)
set(h(1),'Color','r')
set(h(2),'Color','g')
set(h(3),'Color','b')
xlabel('time [s]')
ylabel('Moment [Nm]')

subplot(3,1,3)
h=plot(sol_val.tgrid_col{1},forces_prescribed{i}(:,4:6),'LineWidth',2)
set(h(1),'Color','r')
set(h(2),'Color','g')
set(h(3),'Color','b')
xlabel('time [s]')
ylabel({'Point of force application'; 'in G frame [m]'})

legend({'X','Y','Z'},'Orientation','horizontal','Box','off');

%% Forces and moments splitted
%% Forces
subplot(3,3,1)
h=plot(sol_val.tgrid_col{1},forces_prescribed{i}(:,1),'LineWidth',2)
set(h(1),'Color','r')
set(gca,'XGrid','on')
ylabel('force [N]')
title('X')
subplot(3,3,2)
h=plot(sol_val.tgrid_col{1},forces_prescribed{i}(:,2),'LineWidth',2)
set(h(1),'Color','g')
set(gca,'XGrid','on')
title('Y')
subplot(3,3,3)
h=plot(sol_val.tgrid_col{1},forces_prescribed{i}(:,3),'LineWidth',2)
set(h(1),'Color','b')
xlabel('time [s]')
title('Z')
set(gca,'XGrid','on')

subplot(3,3,4)
h=plot(sol_val.tgrid_col{1},forces_prescribed{i}(:,7),'LineWidth',2)
set(h(1),'Color','r')
ylabel('Moment [Nm]')
set(gca,'XGrid','on')
subplot(3,3,5)
h=plot(sol_val.tgrid_col{1},forces_prescribed{i}(:,8),'LineWidth',2)
set(h(1),'Color','g')
set(gca,'XGrid','on')
subplot(3,3,6)
h=plot(sol_val.tgrid_col{1},forces_prescribed{i}(:,9),'LineWidth',2)
set(h(1),'Color','b')
xlabel('time [s]')
set(gca,'XGrid','on')


subplot(3,3,7)
h=plot(sol_val.tgrid_col{1},forces_prescribed{i}(:,4),'LineWidth',2)
set(h(1),'Color','r')
ylabel({'Point of force application'; 'in G frame [m]'})
set(gca,'XGrid','on')
subplot(3,3,8)
h=plot(sol_val.tgrid_col{1},forces_prescribed{i}(:,5),'LineWidth',2)
set(h(1),'Color','g')
set(gca,'XGrid','on')
subplot(3,3,9)
h=plot(sol_val.tgrid_col{1},forces_prescribed{i}(:,6),'LineWidth',2)
set(h(1),'Color','b')
xlabel('time [s]')
set(gca,'XGrid','on')
% legend({'X','Y','Z'},'Orientation','horizontal','Box','off');

%% Moments and test indentifying bumps
i=1;
for k=1:size(all_Qd2dot_opt{i},1)
        if Options.optInertiaParam
            out_opt{i}(k,:)=full(F([all_QsQdot_opt{i}(k,:)';all_Qd2dot_opt{i}(k,:)';forces_prescribed{i}(k,:)';sol_val.inertiaParam]));
        else
            out_opt{i}(k,:)=full(F([all_QsQdot_opt{i}(k,:)';all_Qd2dot_opt{i}(k,:)';forces_prescribed{i}(k,:)']));
        end
end

i=1;
all_QsQdot_opt_no_vel=all_QsQdot_opt{i};
all_QsQdot_opt_no_vel(:,2:2:end)=0;
for k=1:size(all_Qd2dot_opt{i},1)
    if Options.optInertiaParam
        out_opt_novel{i}(k,:)=full(F([all_QsQdot_opt_no_vel(k,:)';all_Qd2dot_opt{i}(k,:)';forces_prescribed{i}(k,:)';sol_val.inertiaParam]));
    else
        out_opt_novel{i}(k,:)=full(F([all_QsQdot_opt_no_vel(k,:)';all_Qd2dot_opt{i}(k,:)';forces_prescribed{i}(k,:)']));
    end
end

for k=1:size(all_Qd2dot_opt{i},1)
    if Options.optInertiaParam
        out_opt_noacc{i}(k,:)=full(F([all_QsQdot_opt{i}(k,:)';zeros(14,1);forces_prescribed{i}(k,:)';sol_val.inertiaParam]));
    else
        out_opt_noacc{i}(k,:)=full(F([all_QsQdot_opt{i}(k,:)';zeros(14,1);forces_prescribed{i}(k,:)']));
    end
end

for k=1:size(all_Qd2dot_opt{i},1)
    if Options.optInertiaParam
        out_opt_noforce{i}(k,:)=full(F([all_QsQdot_opt{i}(k,:)';all_Qd2dot_opt{i}(k,:)';zeros(9,1);sol_val.inertiaParam]));
    else
        out_opt_noforce{i}(k,:)=full(F([all_QsQdot_opt{i}(k,:)';all_Qd2dot_opt{i}(k,:)';zeros(9,1)]));
    end
end

for k=1:size(all_Qd2dot_opt{i},1)
    if Options.optInertiaParam
        out_opt_noacc{i}(k,:)=full(F([all_QsQdot_opt{i}(k,:)';zeros(14,1);forces_prescribed{i}(k,:)';sol_val.inertiaParam]));
    else
        out_opt_noacc{i}(k,:)=full(F([all_QsQdot_opt{i}(k,:)';zeros(14,1);forces_prescribed{i}(k,:)']));
    end
end

for k=1:size(all_Qd2dot_opt{i},1)
    if Options.optInertiaParam
        out_opt_noforce{i}(k,:)=full(F([all_QsQdot_opt{i}(k,:)';all_Qd2dot_opt{i}(k,:)';zeros(9,1);sol_val.inertiaParam]));
    else
        out_opt_noforce{i}(k,:)=full(F([all_QsQdot_opt{i}(k,:)';all_Qd2dot_opt{i}(k,:)';zeros(9,1)]));
    end
end


subplot(3,1,1)
%hip flexion
plot(sol_val.tgrid_col{1},out_opt{i}(:,8),'LineWidth',4);
hold all;
plot(sol_val.tgrid_col{1},out_opt_novel{i}(:,8),'LineWidth',2);
plot(sol_val.tgrid_col{1},out_opt_noacc{i}(:,8),'LineWidth',2);
plot(sol_val.tgrid_col{1},out_opt_noforce{i}(:,8),'LineWidth',2);
xlabel('time [s]')
ylabel('hip moment [Nm]')

subplot(3,1,2)
%knee flexion
plot(sol_val.tgrid_col{1},out_opt{i}(:,11),'LineWidth',4);
hold all;
plot(sol_val.tgrid_col{1},out_opt_novel{i}(:,11),'LineWidth',2);
plot(sol_val.tgrid_col{1},out_opt_noacc{i}(:,11),'LineWidth',2);
plot(sol_val.tgrid_col{1},out_opt_noforce{i}(:,11),'LineWidth',2);
xlabel('time [s]')
ylabel('knee moment [Nm]')

subplot(3,1,3)
%ankle flexion
plot(sol_val.tgrid_col{1},out_opt{i}(:,12),'LineWidth',4);
hold all;
plot(sol_val.tgrid_col{1},out_opt_novel{i}(:,12),'LineWidth',2);
plot(sol_val.tgrid_col{1},out_opt_noacc{i}(:,12),'LineWidth',2);
plot(sol_val.tgrid_col{1},out_opt_noforce{i}(:,12),'LineWidth',2);
xlabel('time [s]')
ylabel('ankle moment [Nm]')
legend('moment','without vel','without acc','without ext. force','Orientation','horizontal','box','off')

%% Moments and test indentifying bumps for different values of accelerations
i=1;
for k=1:size(all_Qd2dot_opt{i},1)
        out_opt{i}(k,:)=full(F([all_QsQdot_opt{i}(k,:)';all_Qd2dot_opt{i}(k,:)';forces_prescribed{i}(k,:)']));
        out_opt125{i}(k,:)=full(F([all_QsQdot_opt{i}(k,:)';all_Qd2dot_opt{i}(k,:)'*1.25;forces_prescribed{i}(k,:)']));
        out_opt150{i}(k,:)=full(F([all_QsQdot_opt{i}(k,:)';all_Qd2dot_opt{i}(k,:)'*1.5;forces_prescribed{i}(k,:)']));
        out_opt175{i}(k,:)=full(F([all_QsQdot_opt{i}(k,:)';all_Qd2dot_opt{i}(k,:)'*1.75;forces_prescribed{i}(k,:)']));
        out_opt200{i}(k,:)=full(F([all_QsQdot_opt{i}(k,:)';all_Qd2dot_opt{i}(k,:)'*2;forces_prescribed{i}(k,:)']));
end

subplot(3,1,1)
%hip flexion
plot(sol_val.tgrid_col{1},out_opt{i}(:,8),'LineWidth',4);
hold all;
plot(sol_val.tgrid_col{1},out_opt125{i}(:,8),'LineWidth',1);
plot(sol_val.tgrid_col{1},out_opt150{i}(:,8),'LineWidth',1);
plot(sol_val.tgrid_col{1},out_opt175{i}(:,8),'LineWidth',1);
plot(sol_val.tgrid_col{1},out_opt200{i}(:,8),'LineWidth',1);
xlabel('time [s]')
ylabel('hip moment [Nm]')

subplot(3,1,2)
%knee flexion
plot(sol_val.tgrid_col{1},out_opt{i}(:,11),'LineWidth',4);
hold all;
plot(sol_val.tgrid_col{1},out_opt125{i}(:,11),'LineWidth',1);
plot(sol_val.tgrid_col{1},out_opt150{i}(:,11),'LineWidth',1);
plot(sol_val.tgrid_col{1},out_opt175{i}(:,11),'LineWidth',1);
plot(sol_val.tgrid_col{1},out_opt200{i}(:,11),'LineWidth',1);
xlabel('time [s]')
ylabel('knee moment [Nm]')

subplot(3,1,3)
%ankle flexion
plot(sol_val.tgrid_col{1},out_opt{i}(:,12),'LineWidth',4);
hold all;
plot(sol_val.tgrid_col{1},out_opt125{i}(:,12),'LineWidth',1);
plot(sol_val.tgrid_col{1},out_opt150{i}(:,12),'LineWidth',1);
plot(sol_val.tgrid_col{1},out_opt175{i}(:,12),'LineWidth',1);
plot(sol_val.tgrid_col{1},out_opt200{i}(:,12),'LineWidth',1);

xlabel('time [s]')
ylabel('ankle moment [Nm]')
legend({'moment','moment w. acc*1.25','moment w. acc*1.50','moment w. acc*1.75','moment w. acc*2.00'},'Orientation','vertical','box','off')


%% Plot passive properties

subplot(3,1,1)
bar(sol_val.Kstiff);
set(gca,'XTick',1:3,'XTickLabels',{'hip','knee','ankle'});
ylabel('K [Nm/rad]');
subplot(3,1,2)
bar(sol_val.theta0);
set(gca,'XTick',1:3,'XTickLabels',{'hip','knee','ankle'});
ylabel('\theta_0 [rad]');
subplot(3,1,3)
bar(sol_val.Ddamp);
set(gca,'XTick',1:3,'XTickLabels',{'hip','knee','ankle'});
ylabel('D [Nm/(rad/s)]');

%% Plot inertia parameters
c=[0 0 0; 0 0 1; 0 0 0];
subplot(3,3,1);
h=bar([bounds.inertiaParam.lower(1)' sol_val.inertiaParam(1) bounds.inertiaParam.upper(1)']);
h.FaceColor='flat';
h.CData(1,:)=c(1,:);
h.CData(2,:)=c(2,:);
h.CData(3,:)=c(3,:);
set(gca,'XTickLabels',{'LB','m','UB'});
title('mass femur');
ylabel('mass [kg]')
subplot(3,3,2);
h=bar([bounds.inertiaParam.lower(4)' sol_val.inertiaParam(4) bounds.inertiaParam.upper(4)']);
title('mass tibia');
ylabel('mass [kg]')
h.FaceColor='flat';
h.CData(1,:)=c(1,:);
h.CData(2,:)=c(2,:);
h.CData(3,:)=c(3,:);
set(gca,'XTickLabels',{'LB','m','UB'});
subplot(3,3,3);
h=bar([bounds.inertiaParam.lower(5)' sol_val.inertiaParam(7) bounds.inertiaParam.upper(7)']);
title('mass foot');
ylabel('mass [kg]')
h.FaceColor='flat';
h.CData(1,:)=c(1,:);
h.CData(2,:)=c(2,:);
h.CData(3,:)=c(3,:);
set(gca,'XTickLabels',{'LB','m','UB'});

subplot(3,3,4);
h=bar([bounds.inertiaParam.lower(2)' sol_val.inertiaParam(2) bounds.inertiaParam.upper(2)']);
h.FaceColor='flat';
h.CData(1,:)=c(1,:);
h.CData(2,:)=c(2,:);
h.CData(3,:)=c(3,:);
set(gca,'XTickLabels',{'LB','I_z','UB'});
title('I_z femur');
ylabel('I_z [kg m^2]')
subplot(3,3,5);
h=bar([bounds.inertiaParam.lower(5)' sol_val.inertiaParam(5) bounds.inertiaParam.upper(5)']);
h.FaceColor='flat';
h.CData(1,:)=c(1,:);
h.CData(2,:)=c(2,:);
h.CData(3,:)=c(3,:);
set(gca,'XTickLabels',{'LB','I_z','UB'});
title('I_z tibia');
ylabel('I_z [kg m^2]')
subplot(3,3,6);
h=bar([bounds.inertiaParam.lower(8)' sol_val.inertiaParam(8) bounds.inertiaParam.upper(8)']);
h.FaceColor='flat';
h.CData(1,:)=c(1,:);
h.CData(2,:)=c(2,:);
h.CData(3,:)=c(3,:);
set(gca,'XTickLabels',{'LB','I_z','UB'});
title('I_z foot');
ylabel('I_z [kg m^2]')

subplot(3,3,7);
h=bar([bounds.inertiaParam.lower(3)' sol_val.inertiaParam(3) bounds.inertiaParam.upper(3)']);
h.FaceColor='flat';
h.CData(1,:)=c(1,:);
h.CData(2,:)=c(2,:);
h.CData(3,:)=c(3,:);
set(gca,'XTickLabels',{'LB','y','UB'});
title('y_{COM} femur');
ylabel('y_{COM} [m]')
subplot(3,3,8);
h=bar([bounds.inertiaParam.lower(6)' sol_val.inertiaParam(6) bounds.inertiaParam.upper(6)']);
h.FaceColor='flat';
h.CData(1,:)=c(1,:);
h.CData(2,:)=c(2,:);
h.CData(3,:)=c(3,:);
set(gca,'XTickLabels',{'LB','y','UB'});
title('y_{COM} tibia');
ylabel('y_{COM} [m]')
subplot(3,3,9);
h=bar([bounds.inertiaParam.lower(9)' sol_val.inertiaParam(9) bounds.inertiaParam.upper(9)']);
h.FaceColor='flat';
h.CData(1,:)=c(1,:);
h.CData(2,:)=c(2,:);
h.CData(3,:)=c(3,:);
set(gca,'XTickLabels',{'LB','y','UB'});
title('y_{COM} foot');
ylabel('y_{COM} [m]')






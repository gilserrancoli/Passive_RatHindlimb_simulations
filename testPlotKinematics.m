subplot(3,3,1)
plot(t,expdata.(trial_name).kinematics(:,9),'LineWidth',2);
title('hip');
ylabel('angle [deg]')
hold all;

subplot(3,3,4)
plot(t,expdata.(trial_name).kinematics(:,12),'LineWidth',2);
title('knee');
ylabel('angle [deg]')
hold all;

subplot(3,3,7)
plot(t,expdata.(trial_name).kinematics(:,13),'LineWidth',2);
title('ankle');
ylabel('angle [deg]')
hold all;



subplot(3,3,2)
plot(t,expdata.(trial_name).kinematics_v(:,9)*pi/180,'LineWidth',2);
title('hip');
ylabel('angle vel [rad/s]')
hold all;

subplot(3,3,5)
plot(t,expdata.(trial_name).kinematics_v(:,12)*pi/180,'LineWidth',2);
title('knee');
ylabel('angle vel [rad/s]')
hold all;

subplot(3,3,8)
plot(t,expdata.(trial_name).kinematics_v(:,13)*pi/180,'LineWidth',2);
title('ankle');
ylabel('angle vel [rad/s]')
hold all;



subplot(3,3,3)
plot(t,expdata.(trial_name).kinematics_a(:,9)*pi/180,'LineWidth',2);
title('hip');
ylabel('angle acc [rad/s^2]')
xlabel('time[s]')
hold all;

subplot(3,3,6)
plot(t,expdata.(trial_name).kinematics_a(:,12)*pi/180,'LineWidth',2);
title('knee');
ylabel('angle acc [rad/s^2]')
xlabel('time[s]')
hold all;

subplot(3,3,9)
plot(t,expdata.(trial_name).kinematics_a(:,13)*pi/180,'LineWidth',2);
title('ankle');
ylabel('angle acc [rad/s^2]')
xlabel('time[s]')
hold all;



csvdata=csvread('..\..\rat files\kinematics\170503_files\17050301.csv',7,0);
spine_top=csvdata(:,3:5);
spine_bottom=csvdata(:,6:8);
hip_top=csvdata(:,9:11);
hip_middle=csvdata(:,12:14);
hip_bottom=csvdata(:,15:17);
femur_mid=csvdata(:,18:20);
knee=csvdata(:,21:23);
tibia_mid=csvdata(:,24:26);
heel=csvdata(:,27:29);
foot_mid=csvdata(:,30:32);
toe=csvdata(:,33:35);
ref_a=csvdata(:,36:38);

for i=1:size(spine_top,1)
    plot3(spine_top(i,1),spine_top(i,2),spine_top(i,3),'o');
    hold all;
    plot3(spine_bottom(i,1),spine_bottom(i,2),spine_bottom(i,3),'o');
    plot3(hip_top(i,1),hip_top(i,2),hip_top(i,3),'o');
    plot3(hip_middle(i,1),hip_middle(i,2),hip_middle(i,3),'o');
    plot3(hip_bottom(i,1),hip_bottom(i,2),hip_bottom(i,3),'o');
    plot3(femur_mid(i,1),femur_mid(i,2),femur_mid(i,3),'o');
    plot3(knee(i,1),knee(i,2),knee(i,3),'o');
    plot3(tibia_mid(i,1),tibia_mid(i,2),tibia_mid(i,3),'o');
    plot3(heel(i,1),heel(i,2),heel(i,3),'o');
    plot3(foot_mid(i,1),foot_mid(i,2),foot_mid(i,3),'o');
    plot3(toe(i,1),toe(i,2),toe(i,3),'o');
%     plot3(ref_a(i,1),ref_a(i,2),ref_a(i,3),'o');
    plot3([spine_top(i,1) spine_bottom(i,1) hip_top(i,1) hip_bottom(i,1) hip_middle(i,1) ...
        femur_mid(i,1) knee(i,1) tibia_mid(i,1) heel(i,1) foot_mid(i,1) toe(i,1)],...
        [spine_top(i,2) spine_bottom(i,2) hip_top(i,2) hip_bottom(i,2) hip_middle(i,2)...
        femur_mid(i,2) knee(i,2) tibia_mid(i,2) heel(i,2) foot_mid(i,2) toe(i,2)],...
        [spine_top(i,3) spine_bottom(i,3) hip_top(i,3) hip_bottom(i,3) hip_middle(i,3)...
        femur_mid(i,3) knee(i,3) tibia_mid(i,3) heel(i,3) foot_mid(i,3) toe(i,3)],'k');
    pause(0.01);
    clf;
end
xlim([600 1200]);
legend({'spine top','spine bottom','hip top','hip bottom','hip middle','femur mid','knee','tibia mid','heel','foot mid','toe'});

    
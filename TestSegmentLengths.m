aaa=readtable('kinematics_0_pos10.0.csv');
hip=aaa(:,56:58);
knee=aaa(:,62:64);
ankle=aaa(:,65:67);
bbb=readtable('kinematics_1_pos0.0.csv');
hip1=bbb(:,56:58);
knee1=bbb(:,62:64);
ankle1=bbb(:,65:67);
ccc=readtable('kinematics_2_pos-10.0.csv');
hip2=ccc(:,56:58);
knee2=ccc(:,62:64);
ankle2=ccc(:,65:67);

femur0=sqrt(sum((table2array(hip)-table2array(knee)).^2,2));
femur1=sqrt(sum((table2array(hip1)-table2array(knee1)).^2,2));
femur2=sqrt(sum((table2array(hip2)-table2array(knee2)).^2,2));
plot(femur0);
hold all;
plot(femur1);
plot(femur2);

tibia0=sqrt(sum((table2array(knee)-table2array(ankle)).^2,2));
tibia1=sqrt(sum((table2array(knee1)-table2array(ankle1)).^2,2));
tibia2=sqrt(sum((table2array(knee2)-table2array(ankle2)).^2,2));
plot(tibia0);
hold all;
plot(tibia1);
plot(tibia2);


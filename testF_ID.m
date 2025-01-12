IK0=importdata('C:\Gil\Collaborations\MatthewTresch\Zhong_ParameterEstimation\Optimization\DataMay\kinematics\h8\perturb1.mot');

IK.data=repmat(IK0.data(1,:),10,1);
IK.data(:,1)=0:0.0001:0.0009;
IK.labels=IK0.colheaders;
write_motionFile(IK,'test.mot');

import casadi.*
q_out(:,1:2:28)=IK.data(:,2:end)*pi/180;
q_out(:,28)=0;
F = external('F','RightRatHindlimb_Zhong.dll');   
for i=1:size(IK.data,1)
    aux=F([q_out(i,:)'; zeros(14,1);zeros(9,1)]);
    outtt(i,:)=full(aux);
end

%%test 2 with forces
forces=importdata('C:\Gil\Collaborations\MatthewTresch\Zhong_ParameterEstimation\Optimization\DataMay\motion&force\h8\perturb1.mot');
F = external('F','RightRatHindlimb_Zhong.dll');   
for i=1:size(IK.data,1)
    aux=F([q_out(i,:)'; zeros(14,1);forces.data(i,2:end)']);
    outtt(i,:)=full(aux);
end

%%test 3 with only moment
forces_onlyMoment=forces;
forces_onlyMoment.data(:,2:7)=0;
q.data=forces_onlyMoment.data;
q.labels=forces_onlyMoment.colheaders;
q.data(:,1)=0:0.0001:(size(q.data,1)-1)*0.0001;
write_motionFile(q,'test_forces_onlyMoment.mot');

F = external('F','RightRatHindlimb_Zhong.dll');   
for i=1:size(q_out,1)
    aux=F([q_out(i,:)'; zeros(14,1);q.data(i,2:end)']);
    outtt(i,:)=full(aux);
end

%%test 3 with only force
forces_onlyForce=forces;
forces_onlyForce.data(:,8:10)=0;
q.data=forces_onlyForce.data;
q.labels=forces_onlyForce.colheaders;
q.data(:,1)=0:0.0001:(size(q.data,1)-1)*0.0001;
write_motionFile(q,'test_forces_onlyForce.mot');

F = external('F','RightRatHindlimb_Zhong.dll');   
for i=1:size(q_out,1)
    aux=F([q_out(i,:)'; zeros(14,1);q.data(i,2:end)']);
    outtt(i,:)=full(aux);
end

%% Test zero position
q.data=zeros(20,size(IK.data,2));
q.data(1:20,1)=0:0.0001:0.0019;
q.labels=IK0.colheaders;
write_motionFile(q,'test_zeroAngles.mot');

F = external('F','RightRatHindlimb_Zhong.dll');  
clear q_out;
q_out=zeros(20,28);
q_out(:,1:2:28)=q.data(:,2:end);
for i=1:size(q_out,1)
    aux=F([q_out(i,:)'; zeros(14,1);zeros(9,1)]);
    outtt(i,:)=full(aux);
end

%% Test zero position with Y, X or Z force
q.data=zeros(20,10);
q.data(1:20,1)=0:0.0001:0.0019;
q.labels=forces_onlyForce.colheaders;
q.data(1:5,2:4)=repmat([0 1 0],5,1);
q.data(6:10,[2:4])=repmat([1 0 0],5,1);
q.data(11:20,[2:4])=repmat([0 0 1],10,1);
write_motionFile(q,'test_forces_YXZdirection.mot');

F = external('F','RightRatHindlimb_Zhong.dll');  
clear q_out;

for i=1:size(q_out,1)
    aux=F([q_out(i,:)'; zeros(14,1);q.data(i,2:end)']);
    outtt(i,:)=full(aux);
end




%% Test with original data perturb1
IK0=importdata('C:\Gil\Collaborations\MatthewTresch\Zhong_ParameterEstimation\Optimization\DataMay\kinematics\h8\perturb1.mot');

IK.data=IK0.data;
IK.labels=IK0.colheaders;

for i=1:size(IK.data,2)-1
    IK.datafilt=filtfilt(B,A,IK.data(:,i+1));
    qsplines(i)=spline(IK.data(:,1),IK.datafilt);
    qdotsplines(i)=fnder(qsplines(i),1);
    qd2dotsplines(i)=fnder(qsplines(i),2);
end


clear q_out;
q_out(:,1:2:28)=IK.data(:,2:end);
q_out(:,2:2:28)=
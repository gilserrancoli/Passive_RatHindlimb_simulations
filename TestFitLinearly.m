% kindata_file=importdata('C:\Gil\Collaborations\MatthewTresch\Zhong_ParameterEstimation\Optimization\DataSeptember\kinematics\data6_kin.mot');
% pertdata_file=importdata('C:\Gil\Collaborations\MatthewTresch\Zhong_ParameterEstimation\Optimization\DataSeptember\perturbation\data6.mot');
kindata_file=importdata('C:\Gil\Collaborations\MatthewTresch\Zhong_ParameterEstimation\Optimization\DataNovember\kinematics\kinematics_1_pos0.0.mot');
pertdata_file=importdata('C:\Gil\Collaborations\MatthewTresch\Zhong_ParameterEstimation\Optimization\DataNovember\perturbation\motor_1_pos0.0.mot');

N=20;
d=3;
method='legendre';
[tau_root,C,D,B] = CollocationScheme(d,method);

kindata=ProcessKinematics(kindata_file);
% C=strrep(kinfilename,'perturb','');
% sufix=strrep(C,'.mot','');
% 
% trial_name=[strrep(kinfilename,'.mot','')];
%parameterize with splines
% trial_name='trialSept';
trial_name='kinematics_1_pos0_0';
t=kindata.data(1,1):0.0002:kindata.data(end,1);
expdata.(trial_name).kinematics(:,1)=t;
expdata.(trial_name).kinematics_v(:,1)=t;
expdata.(trial_name).kinematics_a(:,1)=t;
[B,A]=butter(3,500/(5000/2));
for j=2:size(kindata.data,2)
    intdata=interp1(kindata.data(:,1),kindata.data(:,j),t,'spline');
    smoothed_kin=smooth(t,intdata,0.5,'rloess');
    smoothed_filt_kin=filtfilt(B,A,smoothed_kin);

    kindata_spline(j-1)=spline(t,smoothed_filt_kin);

    expdata.(trial_name).kinematics(:,j)=ppval(kindata_spline(j-1),t);
    expdata.(trial_name).kinematics_v(:,j)=ppval(fnder(kindata_spline(j-1),1),t);
    expdata.(trial_name).kinematics_a(:,j)=ppval(fnder(kindata_spline(j-1),2),t);

end
expdata.(trial_name).kinematics_labels=kindata.colheaders;

expdata.(trial_name).forces=pertdata_file.data;
expdata.(trial_name).forces_labels=pertdata_file.colheaders;

t0=max(pertdata_file.data(1,1),expdata.(trial_name).kinematics(1,1));
tf=min(pertdata_file.data(end,1),expdata.(trial_name).kinematics(end,1));
tgrid(1:(d+1):(N+1)*(d+1))=t0:((tf-t0)/N):tf;
deltat=tgrid(d+1+1)-tgrid(1);
tgrid(2:(d+1):N*(d+1))=tgrid(1:(d+1):N*(d+1))+tau_root(2)*deltat;
tgrid(3:(d+1):N*(d+1))=tgrid(1:(d+1):N*(d+1))+tau_root(3)*deltat;
tgrid(4:(d+1):N*(d+1))=tgrid(1:(d+1):N*(d+1))+tau_root(4)*deltat;
expdata.(trial_name).tgrid=tgrid;
h=(tgrid(end)-tgrid(1))/N;
expdata.(trial_name).h=h;
expdata.(trial_name).q=interp1(expdata.(trial_name).kinematics(:,1),expdata.(trial_name).kinematics,tgrid);
expdata.(trial_name).q(:,[2:4 8:end])=expdata.(trial_name).q(:,[2:4 8:end])*pi/180;
expdata.(trial_name).qdot=interp1(expdata.(trial_name).kinematics(:,1),expdata.(trial_name).kinematics_v,tgrid);
expdata.(trial_name).qdot(:,[2:4 8:end])=expdata.(trial_name).qdot(:,[2:4 8:end])*pi/180;
expdata.(trial_name).qd2dot=interp1(expdata.(trial_name).kinematics(:,1),expdata.(trial_name).kinematics_a,tgrid);
expdata.(trial_name).qd2dot(:,[2:4 8:end])=expdata.(trial_name).qd2dot(:,[2:4 8:end])*pi/180;
expdata.(trial_name).f=interp1(pertdata_file.data(:,1),pertdata_file.data,tgrid);

% 
% A=[expdata.trialSept.q(:,[9 12 13]) expdata.trialSept.qdot(:,[9 12 13])...
%     expdata.trialSept.qd2dot(:,[9 12 13])];
% b=expdata.(trial_name).f(:,2);
% x=A\b;
% bpred=A*x;
% 
% plot(expdata.trialSept.tgrid,b);
% hold all;
% plot(expdata.trialSept.tgrid,bpred);
% 
% load('sol_val.mat');

Moments=[sol_val.out_opt{1}(:,8) sol_val.out_opt{1}(:,11) sol_val.out_opt{1}(:,12)];
A=[ones(size(expdata.(trial_name).q,1),1) expdata.(trial_name).q(:,[9 12 13]) expdata.(trial_name).qdot(:,[9 12 13])...
    expdata.(trial_name).qd2dot(:,[9 12 13])];
A(1:4:end,:)=[];
b2=Moments;
x2=A\b2;
b2pred=A*x2;

%% Dependency on angle and velocity only
% Cross-dependence

% Solved using least squared minimum
Moments=[sol_val.out_opt{1}(:,8) sol_val.out_opt{1}(:,11) sol_val.out_opt{1}(:,12)];
A=[ones(size(expdata.(trial_name).q,1),1) expdata.(trial_name).q(:,[9 12 13]) expdata.(trial_name).qdot(:,[9 12 13])];
A(1:4:end,:)=[];
b2=Moments;
x2=A\b2;
b2pred=A*x2;

% Solved using lsqlin (with bounds)
C=A;
lb=[-10 -10*ones(1,3) -10*ones(1,3)];
ub=[ pi  zeros(1,3)    10*ones(1,3)];
x0=[  0   0   0  0.1*ones(1,3) 0.1*ones(1,3)];
options=optimset('Display','iter');
for i=1:3
    x = lsqlin(C,b2(:,i),[],[],[],[],lb,ub,[],options);
    sol(i).vals=x;
end
b3pred=A*[sol(1).vals sol(2).vals sol(3).vals];


for i=1:3
    subplot(3,1,i)
    plot(sol_val.tgrid_col{1},b2(:,i),'LineWidth',2);
    hold all;
    plot(sol_val.tgrid_col{1},b2pred(:,i),'LineWidth',2);
    plot(sol_val.tgrid_col{1},b3pred(:,i),'LineWidth',2);
    r=corrcoef(b2(:,i),b2pred(:,i));
    r3=corrcoef(b2(:,i),b3pred(:,i));
    if i==1
        tt=-0.01;
        title('hip')
    elseif i==2
        tt=-0.007;
        title('knee')
    elseif i==3
        tt=-0.006;
        title('ankle')
    end
    text(0,tt,['r = ' , num2str(r(2))])
    text(0,tt+0.001,['r = ' , num2str(r(2))])
end
legend('exp','model unbounded','model bounded');

% Solved with stiffness dependent on one joint
% Moments=[sol_val.out_opt{1}(:,8) sol_val.out_opt{1}(:,11) sol_val.out_opt{1}(:,12)];
% A=[ones(size(expdata.(trial_name).q,1),1) expdata.(trial_name).q(:,[9 12 13]) expdata.(trial_name).qdot(:,[9 12 13])];
Moments=[sol_val.out_opt{1}(:,8)];
A=[ones(size(expdata.(trial_name).q,1),1) expdata.(trial_name).q(:,[9]) expdata.(trial_name).qdot(:,[9])];
A(1:4:end,:)=[];
b2=Moments;
x2(:,1)=A\b2;
b2pred(:,1)=A*x2(:,1);

Moments=[sol_val.out_opt{1}(:,11)];
A=[ones(size(expdata.(trial_name).q,1),1) expdata.(trial_name).q(:,[12]) expdata.(trial_name).qdot(:,[12])];
A(1:4:end,:)=[];
b2=Moments;
x2(:,2)=A\b2;
b2pred(:,2)=A*x2(:,2);

Moments=[sol_val.out_opt{1}(:,12)];
A=[ones(size(expdata.(trial_name).q,1),1) expdata.(trial_name).q(:,[13]) expdata.(trial_name).qdot(:,[13])];
A(1:4:end,:)=[];
b2=Moments;
x2(:,3)=A\b2;
b2pred(:,3)=A*x2(:,3);

% Solved using lsqlin (with bounds)
Moments=[sol_val.out_opt{1}(:,8)];
C=[ones(size(expdata.(trial_name).q,1),1) expdata.(trial_name).q(:,[9]) expdata.(trial_name).qdot(:,[9])];
C(1:4:end,:)=[];
lb=[-10 -10 -10];
ub=[ 10   0   0];
x0=[  0   0   0];
options=optimset('Display','iter');
x3(:,1) = lsqlin(C,Moments,[],[],[],[],lb,ub,[],options);
b3pred(:,1)=A*x3(:,1);

Moments=[sol_val.out_opt{1}(:,11)];
C=[ones(size(expdata.(trial_name).q,1),1) expdata.(trial_name).q(:,[12]) expdata.(trial_name).qdot(:,[12])];
C(1:4:end,:)=[];
lb=[-10 -10 -10];
ub=[ 10   0   0];
x0=[  0   0   0];
options=optimset('Display','iter');
x3(:,2) = lsqlin(C,Moments,[],[],[],[],lb,ub,[],options);
b3pred=A*x3(:,2);

Moments=[sol_val.out_opt{1}(:,12)];
C=[ones(size(expdata.(trial_name).q,1),1) expdata.(trial_name).q(:,[13]) expdata.(trial_name).qdot(:,[13])];
C(1:4:end,:)=[];
lb=[-10 -10 -10];
ub=[ 10   0   0];
x0=[  0   0   0];
options=optimset('Display','iter');
x3(:,3) = lsqlin(C,Moments,[],[],[],[],lb,ub,[],options);
b3pred=A*x3(:,3);


function kindata_out=ProcessKinematics(kindata)
kindata_out=kindata;
v=diff(kindata.data(:,2:end))./diff(kindata.data(:,1));
for i=1:size(v,2)
    
    [out,I]=rmoutliers(v(:,i));
    if any(I)
        aux=kindata.data(~I,i+1);
        aux(end+1)=kindata.data(end,i+1);
        time=kindata.data(~I,1);
        time(end+1)=kindata.data(end,1);
        kindata_out.data(:,i+1)=interp1(time,aux,kindata.data(:,1));
    end
end



end
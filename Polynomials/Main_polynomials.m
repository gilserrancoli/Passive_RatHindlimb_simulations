% This function generates polynomials to approximate muscle-tendon lengths
% and moment arms. The code is from Wouter Aerts and is adapted to be 
% used with CasADi.
%
% Author: Antoine Falisse
%
% Datum: 03/04/2018
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc

saveQdot=1;
savePolynomials=1;

%% Extract time and angles from dummy motion

pathmain = pwd;
name_dummymotion = 'dummy_motion.mot';
path_polynomials = [pathmain];
path_resultsMA = [path_polynomials,'\MomentArms\'];

dummy_motion = importdata([path_polynomials,'\',name_dummymotion]);
% 7 dofs 
% Order of dofs: hip flx, hip add, hip int, knee flx, ankle flx, 
% ankle add, ankle int
q = dummy_motion.data(:,9:15).*(pi/180);

% Generate random numbers between -1000 and 1000 (°/s) 
if saveQdot
    a = -1000;
    b = 1000;
    r1 = (b-a).*rand(size(q,1),1) + a;
    r2 = (b-a).*rand(size(q,1),1) + a;
    r3 = (b-a).*rand(size(q,1),1) + a;
    r4 = (b-a).*rand(size(q,1),1) + a;
    r5 = (b-a).*rand(size(q,1),1) + a;
    r6 = (b-a).*rand(size(q,1),1) + a;
    r7 = (b-a).*rand(size(q,1),1) + a;
    
    r = [r1,r2,r3,r4,r5,r6,r7];
    qdot = zeros(size(q));
    qdot = r.*(pi/180);
    dummy_qdot = qdot;
    save([path_polynomials,'dummy_qdot.mat'],'dummy_qdot');
end
load([path_polynomials,'dummy_qdot.mat']);
qdot = dummy_qdot(:,:);


%% Import data
% lMT
lMT = importdata([path_resultsMA,'RatHindlimbRight-scaled_MuscleAnalysis_Length.sto']);
% hip flexion
MA.hip.flex = importdata([path_resultsMA,'RatHindlimbRight-scaled_MuscleAnalysis_MomentArm_hip_flx.sto']);
% hip adduction
MA.hip.add = importdata([path_resultsMA,'RatHindlimbRight-scaled_MuscleAnalysis_MomentArm_hip_add.sto']);
% hip int
MA.hip.int = importdata([path_resultsMA,'RatHindlimbRight-scaled_MuscleAnalysis_MomentArm_hip_int.sto']);
% knee flx 
MA.knee.flex = importdata([path_resultsMA,'RatHindlimbRight-scaled_MuscleAnalysis_MomentArm_knee_flx.sto']);
% ankle flx
MA.ankle.flex = importdata([path_resultsMA,'RatHindlimbRight-scaled_MuscleAnalysis_MomentArm_ankle_flx.sto']);
% ankle add
MA.ankle.add = importdata([path_resultsMA,'RatHindlimbRight-scaled_MuscleAnalysis_MomentArm_ankle_add.sto']);
% ankle int
MA.ankle.int = importdata([path_resultsMA,'RatHindlimbRight-scaled_MuscleAnalysis_MomentArm_ankle_int.sto']);

%% Organize MuscleData
    MuscleData.dof_names = dummy_motion.colheaders([9:15]); 
    muscleNames = {'BFa','BFp','CF','STp',...
        'STa','SM','QF','TFL','GMa','GMe','GMi','Pir',...
        'GP','GA','AL','AM','AB','Pec','IP',...
        'OE','OI','GS','GI','RF','VL','VI',...
        'VM','MG','LG','Pla','Sol','TA','EDL',...
        'TP','FDL','FHL','Per','Pop'};
   
    MuscleData.muscle_names = muscleNames;
    for m = 1:length(muscleNames)
        MuscleData.lMT(:,m)     = lMT.data(:,strcmp(lMT.colheaders,muscleNames{m}));       % lMT    
        MuscleData.dM(:,m,1)    = MA.hip.flex.data(:,strcmp(lMT.colheaders,muscleNames{m}));    % hip_flex
        MuscleData.dM(:,m,2)    = MA.hip.add.data(:,strcmp(lMT.colheaders,muscleNames{m}));     % hip_add
        MuscleData.dM(:,m,3)    = MA.hip.int.data(:,strcmp(lMT.colheaders,muscleNames{m}));     % hip_rot
        MuscleData.dM(:,m,4)    = MA.knee.flex.data(:,strcmp(lMT.colheaders,muscleNames{m}));   % knee flx
        MuscleData.dM(:,m,5)    = MA.ankle.flex.data(:,strcmp(lMT.colheaders,muscleNames{m}));  % ankle flx
        MuscleData.dM(:,m,6)    = MA.ankle.add.data(:,strcmp(lMT.colheaders,muscleNames{m}));   % ankle add
        MuscleData.dM(:,m,7)    = MA.ankle.int.data(:,strcmp(lMT.colheaders,muscleNames{m}));   % ankle int
    end
    MuscleData.q = q;
    MuscleData.qdot = qdot;


%% Call PolynomialFit
[muscle_spanning_joint_INFO,MuscleInfo] = ...
    PolynomialFit(MuscleData);

if savePolynomials
    save MuscleData_subject MuscleData
    save muscle_spanning_joint_INFO_subject muscle_spanning_joint_INFO
    save MuscleInfo_subject MuscleInfo
end



%% Create CasADi functions
import casadi.*
% Order mobilities: hip_flex, hip_add, hip_rot, knee_angle, ankle-angle, 
load muscle_spanning_joint_INFO_subject.mat
load MuscleInfo_subject.mat
NMuscle = length(MuscleInfo.muscle);
q_leg = 7;
qin     = SX.sym('qin',1,q_leg);
qdotin  = SX.sym('qdotin',1,q_leg);
lMT     = SX(NMuscle,1);
vMT     = SX(NMuscle,1);
dM      = SX(NMuscle,q_leg);
for i=1:NMuscle     
    index_dof_crossing           = find(muscle_spanning_joint_INFO(i,1:7)==1);
    order                        = MuscleInfo.muscle(i).order;
    [mat,diff_mat_q]             = n_art_mat_3_cas_SX(qin(1,index_dof_crossing),order);
    lMT(i,1)                     = mat*MuscleInfo.muscle(i).coeff;
    vMT(i,1)                     = 0;
    dM(i,1:q_leg)          = 0;
    nr_dof_crossing              = length(index_dof_crossing); 
    for dof_nr = 1:nr_dof_crossing
        dM(i,index_dof_crossing(dof_nr)) = (-(diff_mat_q(:,dof_nr)))'*MuscleInfo.muscle(i).coeff;
        vMT(i,1) = vMT(i,1) + (-dM(i,index_dof_crossing(dof_nr))*qdotin(1,index_dof_crossing(dof_nr)));
    end 
end
f_lMT_vMT_dM = Function('f_lMT_vMT_dM',{qin,qdotin},{lMT,vMT,dM});

%% Check results
load MuscleData_subject.mat
lMT_out = zeros(size(q,1),NMuscle);
vMT_out = zeros(size(q,1),NMuscle);
dM_out = zeros(size(q,1),NMuscle,q_leg);
for i = 1:size(q,1)
    [out1_r,out2_r,out3_r] = f_lMT_vMT_dM(MuscleData.q(i,:),MuscleData.qdot(i,:));
    lMT_out(i,:) = full(out1_r);
    vMT_out(i,:) = full(out2_r);
    dM_out(i,:,1) = full(out3_r(:,1));
    dM_out(i,:,2) = full(out3_r(:,2));
    dM_out(i,:,3) = full(out3_r(:,3));
    dM_out(i,:,4) = full(out3_r(:,4));
    dM_out(i,:,5) = full(out3_r(:,5));   
    dM_out(i,:,6) = full(out3_r(:,6));
    dM_out(i,:,7) = full(out3_r(:,7));
end

%% lMT

figure()
subplot(4,4,1)
scatter(MuscleData.q(:,4),lMT_out(:,10)); hold on;
scatter(MuscleData.q(:,4),MuscleData.lMT(:,10));
xlabel('q knee');
title('GMe');
subplot(4,4,2)
scatter(MuscleData.q(:,4),lMT_out(:,29)); hold on;
scatter(MuscleData.q(:,4),MuscleData.lMT(:,29));
xlabel('q knee');
title('LG');
subplot(4,4,3)
scatter(MuscleData.q(:,4),lMT_out(:,30)); hold on;
scatter(MuscleData.q(:,4),MuscleData.lMT(:,30));
xlabel('q knee');
title('PLA');


subplot(4,4,4)
scatter(MuscleData.q(:,5),lMT_out(:,31)); hold on;
scatter(MuscleData.q(:,5),MuscleData.lMT(:,31));
xlabel('q ankle');
title('Sol');
subplot(4,4,5)
scatter(MuscleData.q(:,5),lMT_out(:,34)); hold on;
scatter(MuscleData.q(:,5),MuscleData.lMT(:,34));
ylabel('q ankle');
title('TP');
subplot(4,4,6)
scatter3(MuscleData.q(:,5),MuscleData.q(:,6),lMT_out(:,35)); hold on;
scatter3(MuscleData.q(:,5),MuscleData.q(:,6),MuscleData.lMT(:,35));
xlabel('q knee');
ylabel('q ankle');
title('FDL');
subplot(4,4,7)
scatter3(MuscleData.q(:,5),MuscleData.q(:,6),lMT_out(:,36)); hold on;
scatter3(MuscleData.q(:,5),MuscleData.q(:,6),MuscleData.lMT(:,36));
xlabel('q knee');
ylabel('q ankle');
title('FHL');
subplot(4,4,8)
scatter3(MuscleData.q(:,5),MuscleData.q(:,6),lMT_out(:,37)); hold on;
scatter3(MuscleData.q(:,5),MuscleData.q(:,6),MuscleData.lMT(:,37));
xlabel('q knee');
ylabel('q ankle');
title('Per');
subplot(4,4,9)
scatter3(MuscleData.q(:,5),MuscleData.q(:,6),lMT_out(:,38)); hold on;
scatter3(MuscleData.q(:,5),MuscleData.q(:,6),MuscleData.lMT(:,38));
xlabel('q knee');
ylabel('q ankle');
title('Pop');

subplot(4,4,14)
scatter3(MuscleData.q(:,5),MuscleData.q(:,6),lMT_out(:,4)); hold on;
scatter3(MuscleData.q(:,5),MuscleData.q(:,6),MuscleData.lMT(:,43));
xlabel('q knee');
ylabel('q ankle');
title('STp');
legend('Polynomial','Model');
title('lMT');

%% Assert results
for i = 1:NMuscle  
    assertLMT(:,i) = abs(lMT_out(:,i) - MuscleData.lMT(:,i));
    assertdM.hip.flex(:,i) = abs(dM_out(:,i,1) - MuscleData.dM(:,i,1));
    assertdM.hip.add(:,i) = abs(dM_out(:,i,2) - MuscleData.dM(:,i,2));
    assertdM.hip.int(:,i) = abs(dM_out(:,i,3) - MuscleData.dM(:,i,3));
    assertdM.knee_flex(:,i) = abs(dM_out(:,i,4) - MuscleData.dM(:,i,4));
    assertdM.ankle.flex(:,i) = abs(dM_out(:,i,5) - MuscleData.dM(:,i,5));
    assertdM.ankle.add(:,i) = abs(dM_out(:,i,6) - MuscleData.dM(:,i,6));
    assertdM.ankle.int(:,i) = abs(dM_out(:,i,7) - MuscleData.dM(:,i,7));
end

assertLMTmax_r = max(max(assertLMT));
assertdM.hip.flexmax = max(max(assertdM.hip.flex));
assertdM.hip.addmax = max(max(assertdM.hip.add));
assertdM.hip.rotmax = max(max(assertdM.hip.int));
assertdM.kneemax = max(max(assertdM.knee_flex));
assertdM.ankle.flexmax = max(max(assertdM.ankle.flex));
assertdM.ankle.addmax = max(max(assertdM.ankle.add));
assertdM.ankle.intmax = max(max(assertdM.ankle.int));

 
import casadi.*
method='legendre';

%% Options
Options.useRigidTendon=0;
Options.optimizelM0=1;
Options.optimizelTs=1;
Options.optimizefiberdamping=1;
Options.optimizetendondamping=0;
Options.optimizePassiveJointEl=0;
Options.dofs_to_track=[1 1 1 1 1 1 1]; %1. hip flex 2. hip add 3. hip int 4. knee flex 5. ankle flex 6. ankle add 7. ankle int
Options.useOptimizedIG=0;
Options.trialstotrack='all'; %if 'all', all trials are taken
Options.penalizeHighFTtilde=1;
Options.penalizeFTtot=0;
Options.penalizeoutoflMtilde1=0;

%% 
N=6;
d=3;
[tau_root,C,D,B] = CollocationScheme(d,method);

trials={'h8'};
expdata=LoadData(N,d,tau_root,trials);
nmuscles=38;
muscle_names={'BFa','BFp','CF','STp','STa','SM','QF','TFL','GMa','GMe',...
    'GMi','Pir','GP','GA','AL','AM','AB','Pec','IP','OE','OI','GS','GI',...
    'RF','VL','VI','VM','MG','LG','Pla','Sol','TA','EDL','TP','FDL','FHL',...
    'Per','Pop'};

%% Create CasADi function for Hill-muscle relations
% Function for Hill-equilibrium
if Options.useRigidTendon
    FM          = SX(nmuscles,1); %muscle fiber forces
    lMtilde     = SX.sym('lMtilde',nmuscles); % Normalized fiber lengths
    FMactFL     = SX.sym('FMactFL',nmuscles); % Normalized force-length term
    FMactFV     = SX.sym('FMactFV',nmuscles); % Normalized force-velocity term
else
    FTtilde     = SX.sym('FTtilde',nmuscles); % Normalized tendon forces
    dFTtilde    = SX.sym('dFTtilde',nmuscles); % Time derivative tendon forces
    % tension_SX  = SX.sym('tension',NMuscle); % Tensions not used here
    vMmax       = SX(nmuscles,1); % Maximum contraction velocities
    Fpetilde    = SX(nmuscles,1); % Normalized passive forces
    lMtilde     = SX(nmuscles,1); % Normalized fiber lengths

end
a           = SX.sym('a',nmuscles); % Muscle activations
lMT         = SX.sym('lMT',nmuscles); % Muscle-tendon lengths
vMT         = SX.sym('vMT',nmuscles); % Muscle-tendon velocities
Hilldiff    = SX(nmuscles,1); % Hill-equilibrium   
Fce         = SX(nmuscles,1); % Contractile element forces
Fiso        = SX(nmuscles,1); % Normalized forces from force-length curve 
FT          = SX(nmuscles,1); % Tendon forces    
lTtilde     = SX(nmuscles,1); % Normalized tendon lengths

% Parameters of force-length-velocity curves
load Fvparam
load Fpparam
load Faparam
load MTparam %Load parameter values from osim model (FM0, lM0, alphao, vMax)

lTs     = SX.sym('lTs',1,nmuscles); %tendon slack length
if Options.optimizelM0
    lM0     = SX.sym('lM0',1,nmuscles); %optimal fiber length
else
    lM0     = MTparam(2,:);
end
fiber_damping=SX.sym('fiber_damping',1,1);
tendon_damping=SX.sym('tendon_damping',1,1);
kT=35; %tendon stiffness, not optimized for now
shift=0; %not shifted, if so, apply getShift

current_folder=pwd;
cd('MuscleModel');
if Options.useRigidTendon
    [FT, FM, lMtilde, FMactFL, FMactFV, FMpas, cos_alpha] = ...
            HillModel_RigidTendon(a',lMT',vMT',MTparam,lM0,lTs,fiber_damping);
else
    for m = 1:nmuscles
        [Hilldiff(m),FT(m),Fce(m),Fiso(m),vMmax(m),Fpetilde(m),lMtilde(m),lTtilde(m)] = ...
            ForceEquilibrium_FtildeState(a(m),FTtilde(m),dFTtilde(m),...
            lMT(m),vMT(m),MTparam(1,m),lM0(m),lTs(m),MTparam(4,m),MTparam(5,m),Fvparam,Fpparam,Faparam,fiber_damping,tendon_damping);
            %tension_SX(m));
    end
end
cd(current_folder);
if Options.useRigidTendon
    if (Options.optimizelM0)&&(Options.optimizelTs)
        if Options.optimizefiberdamping
            f_forceRigidTendon=Function('f_forceRigidTendon',...
                {a,lMT,vMT,lM0,lTs,fiber_damping},{FT,FM,lMtilde});
        else
            f_forceRigidTendon=Function('f_forceRigidTendon',...
                {a,lMT,vMT,lM0,lTs},{FT,FM,lMtilde});
        end
    else
        %to check if it is enough the previous function
        keyboard;
    end
else
    if (Options.optimizelM0)&&(Options.optimizelTs)
        f_forceEquilibrium_FtildeState = ...
            Function('f_forceEquilibrium_FtildeState',{a,FTtilde,dFTtilde,...
            lMT,vMT,lTs,lM0,fiber_damping,tendon_damping},{Hilldiff,FT,Fce,Fiso,vMmax,lMtilde,lTtilde});
    elseif (~Options.optimizelM0)&&(Options.optimizelTs)
        f_forceEquilibrium_FtildeState = ...
            Function('f_forceEquilibrium_FtildeState',{a,FTtilde,dFTtilde,...
            lMT,vMT,lTs,fiber_damping,tendon_damping},{Hilldiff,FT,Fce,Fiso,vMmax,lMtilde,lTtilde});
    else
    end
end

%% Load external function
F = external('F','RightRatHindlimb_Zhong.dll');   

%% Prepare polynomials
% Indices of the muscles actuating the different joints for later use
pathpolynomial = [pwd,'/Polynomials'];
addpath(genpath(pathpolynomial));
load([pathpolynomial,'/muscle_spanning_joint_INFO_subject.mat']);
load([pathpolynomial,'/MuscleInfo_subject.mat']);
pathmusclefunctions=[pwd,'/Model/MuscleModel'];
addpath(genpath(pathmusclefunctions));
[~,mai] = MomentArmIndices_3D(muscle_names,...
    muscle_spanning_joint_INFO);

pathcasadi_functions=[pwd,'/VariousFunctions/'];
addpath(pathcasadi_functions);
NMuscle_pol=nmuscles;
ndofs=7; %ndofs to construct polynomials must be 7 here
CasADi_functions;


%% Formulation of optimal control problem
ndofs=sum(Options.dofs_to_track); %3 hip, 1 knee 3 ankle
dofs_in_ID=[9:15]; %hip flx, hip add, hip int, knee flex, ankle flex, ankle add, ankle int

%Define scaling
scaling.lM0=0.05;
scaling.lTs=0.05;
scaling.dFTtilde=100;
scaling.FTtilde=5;
scaling.vMtilde=1;
scaling.lMtilde=1;
scaling.tact=0.02;
scaling.tdact=0.02;
scaling.q=1;
scaling.qdot=100;
scaling.qd2dot=1000;
scaling.QsQdots=repmat([scaling.q scaling.qdot],1,ndofs);
scaling.T=10;

if Options.useOptimizedIG
end

W.lM0lit=0.00001;
W.min_maxa=0.1;
W.state_for_reg=1e-6;
W.mindstate=0.01; %0.1
W.qtrack=1; %10
W.qdottrack=1;
W.penalizeHighFTtilde=10;
W.Kstiff=10;
W.Ddamp=10;
W.theta0=10;
W.penalizeFTtot=1e-4;
W.penalizeoutoflMtilde1=0;

ParallelMode='thread';
NThreads=8;

%% Define bounds
bounds.FTtilde.lower=0*ones(nmuscles,1);
bounds.FTtilde.upper=5*ones(nmuscles,1)/scaling.FTtilde;
load('lMTmax.mat');
bounds.lTs.lower=(1e-4)*ones(nmuscles,1)/scaling.lTs;
bounds.lTs.upper=lMTmax'/scaling.lTs;
bounds.lM0.lower=(lMTmax'/4)/scaling.lM0;
bounds.lM0.upper=ones(nmuscles,1)*0.06/scaling.lM0;
bounds.dFTtilde.lower = -1*ones(nmuscles,1);
bounds.dFTtilde.upper = 1*ones(nmuscles,1);
bounds.fiberdamping.lower=1e-3;
bounds.fiberdamping.upper=10;
bounds.tendondamping.lower=1e-3;
bounds.tendondamping.upper=200;
% bounds.tact.lower=0.005/scaling.tact;
% bounds.tact.upper=0.2/scaling.tact;
% bounds.tdact.lower=0.005/scaling.tdact;
% bounds.tdact.upper=0.2/scaling.tdact;
bounds.q.lower=-6*ones(1,ndofs)/scaling.q; 
bounds.q.upper=6*ones(1,ndofs)/scaling.q;
bounds.qdot.lower=-100/scaling.qdot;
bounds.qdot.upper=100/scaling.qdot;
bounds.QsQdots.lower(1:2:ndofs*2)=bounds.q.lower;
bounds.QsQdots.lower(2:2:ndofs*2)=repmat(bounds.qdot.lower,ndofs,1);
bounds.QsQdots.upper(1:2:ndofs*2)=bounds.q.upper;
bounds.QsQdots.upper(2:2:ndofs*2)=bounds.qdot.upper;
bounds.qd2dot.lower=-1000/scaling.qd2dot;
bounds.qd2dot.upper=1000/scaling.qd2dot;
bounds.K.lower=0;
bounds.K.upper=10;
bounds.D.lower=0;
bounds.D.upper=10;
bounds.theta0.lower=[0.0135, -0.1651,-0.0637,-1.1571,-0.1768]'-5*pi/180;
bounds.theta0.upper=[0.0135, -0.1651,-0.0637,-1.1571,-0.1768]'+5*pi/180;

%% Define guesses
guess.FTtilde=zeros(nmuscles,1)/scaling.FTtilde;
guess.fiberdamping=1;
guess.tendondamping=1;
if Options.useOptimizedIG
    guess.lM0=IG.sol_val.lM0';
    guess.lTs=IG.sol_val.lTs';
    guess.lTs(31:32)=0.015/scaling.lTs;
else
    guess.lM0=MTparam(2,:)/scaling.lM0;
    fnames=fieldnames(expdata);
    lMT_guess_at0=ComputelMTguess(expdata.(fnames{1}).q,expdata.(fnames{1}).qdot,f_lMT_vMT_dM,Options);
    if Options.useRigidTendon
        initlTs=(lMT_guess_at0-MTparam(2,:).*sqrt(1-MTparam(4,:)))/2;
    else
        initlTs=(lMT_guess_at0-MTparam(2,:));
    end
    initlTs(initlTs<0)=bounds.lTs.lower(initlTs<0)*scaling.lTs;
    guess.lTs=initlTs/scaling.lTs;
end
guess.vMtilde=zeros(nmuscles,1)/scaling.vMtilde;
guess.Kstiff=zeros(3,1);
guess.Ddamp=zeros(3,1);
guess.theta0=[0.0135, -0.1651,-0.0637,-1.1571,-0.1768]'; % computed as the mean of all initial values, TO BE REcalculated

%populate names of the constraints
g1_names=[];
g2_names=[];
g2_names_ineq=[];
g3_names=[];

%% Start with an empty optimization
opti = casadi.Opti();
if Options.optimizelTs
    lTs=opti.variable(nmuscles,1);
    opti.subject_to(bounds.lTs.lower < lTs < bounds.lTs.upper);
    opti.set_initial(lTs,guess.lTs);
    g1_names=[g1_names; repmat({'lts_bounds'},nmuscles,1)];
    lTs_k=MX.sym('lTs_k',nmuscles,1);
else 
    lTs     = MTparam(3,:);
    lTs=lTs/scaling.lTs;
end
if Options.optimizelM0
    lM0=opti.variable(nmuscles,1);
    opti.subject_to(bounds.lM0.lower < lM0 < bounds.lM0.upper);
    opti.set_initial(lM0,guess.lM0);
    g1_names=[g1_names; repmat({'lM0_bounds'},nmuscles,1)];
    lM0_k=MX.sym('lM0_k',nmuscles,1);
else 
    lTs     = MTparam(2,:);
    lM0=lM0'/scaling.lM0;
end
if Options.optimizefiberdamping
    fiber_damping=opti.variable(1,1);
    opti.subject_to(bounds.fiberdamping.lower < fiber_damping < ...
        bounds.fiberdamping.upper);
    opti.set_initial(fiber_damping,guess.fiberdamping);
    g1_names=[g1_names; 'fiber_damping_bounds'];
else
    fiber_damping=0.01;
end
fiber_damping_k=MX.sym('fiber_damping_k',1,1);

if Options.useRigidTendon==1
else
    if Options.optimizetendondamping
        tendon_damping=opti.variable(1,1);
        opti.subject_to(bounds.tendondamping.lower < tendon_damping < ...
            bounds.tendondamping.upper);
        opti.set_initial(tendon_damping,guess.tendondamping);
        g1_names=[g1_names; 'tendon_damping_bounds'];
    else
        tendon_damping=0;
    end
end
tendon_damping_k=MX.sym('tendon_damping_k',1,1);

if Options.optimizePassiveJointEl
    Kstiff=opti.variable(3,1);
    opti.subject_to(bounds.K.lower<= Kstiff <= bounds.K.upper);
    opti.set_initial(Kstiff,guess.Kstiff);
    g1_names=[g1_names; repmat({'Kstiff_bounds'},3,1)];
    Ddamp=opti.variable(3,1);
    opti.subject_to(bounds.D.lower<=Ddamp <= bounds.D.upper);
    opti.set_initial(Ddamp,guess.Ddamp);
    g1_names=[g1_names; repmat({'Ddamp_bounds'},3,1)];
    theta0=opti.variable(5,1);
    opti.subject_to(bounds.theta0.lower<= theta0 <= bounds.theta0.upper);
    opti.set_initial(theta0,guess.theta0);
    g1_names=[g1_names; repmat({'theta0_bounds'},3,1)];
else
    Kstiff=zeros(3,1);
    Ddamp=zeros(3,1);
    theta0=zeros(5,1);
end
Kstiff_k=MX.sym('Kstiff_k',3,1);
Ddamp_k=MX.sym('Ddamp_k',3,1);
theta0_k=MX.sym('theta0_k',5,1);

J=0;

if Options.optimizelM0
    J=J+W.lM0lit*sum((((lM0*scaling.lM0)-MTparam(2,:)')./MTparam(2,:)').^2);
    J_lM0=W.lM0lit*sum((((lM0*scaling.lM0)-MTparam(2,:)')./MTparam(2,:)').^2);
end


if Options.optimizePassiveJointEl
    J=J+W.Kstiff*sum(Kstiff.^2);
    J=J+W.Ddamp*sum(Ddamp.^2);
    J=J+W.theta0*sum((theta0-guess.theta0).^2);
    J_passjoint=W.Kstiff*sum(Kstiff.^2)+W.Ddamp*sum(Ddamp.^2)+W.theta0*sum((theta0-guess.theta0).^2);
end

clear FTtilde
clear Hilldiff
clear FT
clear lMT
clear Hilldiff;
clear FT;
clear a;

t_col_grid=1:N*(d+1);
t_col_grid(1:4:end)=[];
J_i=0;

nametrials=fieldnames(expdata);

% for i=1:length(nametrials)
for i=[1 3 6:2:12]
    eq_constr={};
    ineq_constr={};
    J_d=0;
    J1{i}=0;
    J2{i}=0;
    J3{i}=0;
    J4{i}=0;
    J5{i}=0;
    
    
    t0=expdata.(nametrials{i}).q(1,1);
    tf=expdata.(nametrials{i}).q(end,1);
    
    tgrid{i}=[];
    tgrid{i}(1:(d+1):(N+1)*(d+1))=t0:((tf-t0)/N):tf;
    deltat=tgrid{i}(d+1+1)-tgrid{i}(1);
    tgrid{i}(2:(d+1):N*(d+1))=tgrid{i}(1:(d+1):N*(d+1))+tau_root(2)*deltat;
    tgrid{i}(3:(d+1):N*(d+1))=tgrid{i}(1:(d+1):N*(d+1))+tau_root(3)*deltat;
    tgrid{i}(4:(d+1):N*(d+1))=tgrid{i}(1:(d+1):N*(d+1))+tau_root(4)*deltat;

    h{i}=(tgrid{i}(end)-tgrid{i}(1))/N;
    
    
    FT=MX.zeros(1,nmuscles);
    FTtilde{i}=opti.variable(nmuscles,N+1);
    opti.subject_to(bounds.FTtilde.lower <= FTtilde{i} <= ...
        bounds.FTtilde.upper);
    g1_names=[g1_names; repmat({'FTtilde bounds k'},nmuscles*(N+1),1)];
    opti.set_initial(FTtilde{i},repmat(guess.FTtilde,1,N+1));
    FTtilde_k=MX.sym('FTtilde_k',nmuscles,1);
    
    FTtilde_col{i}=opti.variable(nmuscles,d*N);
    opti.subject_to(bounds.FTtilde.lower <= FTtilde_col{i} <= ...
        bounds.FTtilde.upper);
    g1_names=[g1_names; repmat({'FTtilde bounds j'},nmuscles*d*N,1)];
    opti.set_initial(FTtilde_col{i},repmat(guess.FTtilde,1,d*N));
    FTtilde_j=MX.sym('FTtilde_j',nmuscles,d);

    dFTtilde_col{i}=opti.variable(nmuscles,d*N);
    opti.subject_to(bounds.dFTtilde.lower <= dFTtilde_col{i} <= ...
        bounds.dFTtilde.upper);
    g1_names=[g1_names; repmat({'dFTtilde bounds j'},nmuscles*d*N,1)];
    opti.set_initial(dFTtilde_col{i},zeros(nmuscles,d*N));
    dFTtilde_j=MX.sym('dFTtilde_j',nmuscles,d);
    
    QsQdot_prescribed{i}(:,1:2:7*2)=expdata.(nametrials{i}).q(:,2:8); %pelvis dofs and sacroiliac_flx are prescribed
    QsQdot_prescribed{i}(:,2:2:7*2)=expdata.(nametrials{i}).qdot(:,2:8);
    Qd2dot_prescribed{i}=expdata.(nametrials{i}).qd2dot(:,2:8);
    guess.QsQdots{i}(:,1:2:ndofs*2)=expdata.(nametrials{i}).q(:,9:15)/scaling.q;
    guess.QsQdots{i}(:,2:2:ndofs*2)=expdata.(nametrials{i}).qdot(:,9:15)/scaling.qdot;
    guess.Qd2dots{i}=expdata.(nametrials{i}).qd2dot(:,9:15)/scaling.qd2dot;
    
    QsQdots{i}=opti.variable(ndofs*2,N+1);
    opti.subject_to(bounds.QsQdots.lower' < QsQdots{i} < bounds.QsQdots.upper');
    g1_names=[g1_names; repmat({'QsQdot bounds'},ndofs*2*(N+1),1)];
    opti.set_initial(QsQdots{i},guess.QsQdots{i}(1:(d+1):end,:)');
    QsQdots_k=MX.sym('QsQdots_k',ndofs*2,1);
    
    QsQdots_col{i}=opti.variable(ndofs*2,d*N);
    opti.subject_to(bounds.QsQdots.lower' < QsQdots_col{i} < bounds.QsQdots.upper');
    opti.set_initial(QsQdots_col{i},guess.QsQdots{i}(t_col_grid,:)');
    g1_names=[g1_names; repmat({'QsQdot j bounds'},ndofs*2*N*d,1)];
    QsQdots_j=MX.sym('QsQdots_j',ndofs*2,d);
    QsQdots_prescribed_j=MX.sym('QsQdots_prescribed_j',14*2-ndofs*2,d);
    QsQdots_totrack_j=MX.sym('QsQdots_to_track_j',ndofs*2,d);
     
    Qd2dots_col{i}=opti.variable(ndofs,d*N);
    opti.subject_to(bounds.qd2dot.lower < Qd2dots_col{i} < bounds.qd2dot.upper);
    opti.set_initial(Qd2dots_col{i},guess.Qd2dots{i}(t_col_grid,:)');
    g1_names=[g1_names; repmat({'Qd2dots j bounds'},ndofs*N*d,1)];
    Qd2dots_j=MX.sym('Qd2dots_j',ndofs,d);
    Qd2dots_prescribed_j=MX.sym('Qd2dots_prescribed_j',14-ndofs,d);
    
    forces_prescribed_j=MX.sym('forces_prescribed_j',9,d);
    forces_prescribed{i}=expdata.(nametrials{i}).f(t_col_grid,2:end);
    
    if Options.useRigidTendon
%         akj=MX.sym('akj',nmuscles,1);    
        akj=zeros(nmuscles,d);
        lMTk=MX.sym('lMTk',nmuscles,1);
        vMTk=MX.sym('vMTk',nmuscles,1);
    else
%         akj=MX.sym('akj',nmuscles,1);    
        akj=zeros(nmuscles,d);   
        lMTk=MX.sym('lMTk',nmuscles,d);
        vMTk=MX.sym('vMTk',nmuscles,d);
    end
    
    out_d{i}=MX.zeros(14,d);
    for j=1:d
        %Get moment arms and muscle-tendon lengths at that frame
        all_qsleg=QsQdots_j(1:2:end,j)'.*scaling.QsQdots(1:2:end);
        all_qdotsleg=QsQdots_j(2:2:end,j)'.*scaling.QsQdots(2:2:end);
        [lMTj{j},vMTj{j},MAj_aux] =  f_lMT_vMT_dM(all_qsleg,all_qdotsleg); 

        MAj{j}.hip_flex   =  MAj_aux(mai(1).mus',1);
        MAj{j}.hip_add    =  MAj_aux(mai(2).mus',2);
        MAj{j}.hip_int    =  MAj_aux(mai(3).mus',3);
        MAj{j}.knee_flex  =  MAj_aux(mai(4).mus',4);
        MAj{j}.ankle_flex =  MAj_aux(mai(5).mus',5);  
        MAj{j}.ankle_add  =  MAj_aux(mai(6).mus',6); 
        MAj{j}.ankle_int  =  MAj_aux(mai(7).mus',7); 
        
        [Hilldiff_j{j},FT_j{j},~,~,Fp_j{j},lMtilde_j{j},lTtilde_j{j}]=f_forceEquilibrium_FtildeState(...
                           akj(:,j),FTtilde_j(:,j)*scaling.FTtilde,dFTtilde_j(:,j)*scaling.dFTtilde,lMTj{j},vMTj{j},lTs_k*scaling.lTs,lM0_k*scaling.lM0,fiber_damping_k,tendon_damping_k);
        
        all_QsQdot=MX.zeros(1,28);    
        all_QsQdot(1,1:14)=QsQdots_prescribed_j(1:14,j);
        all_QsQdot(1,15:28)=QsQdots_j(:,j)'.*scaling.QsQdots;
        
        all_Qd2dot=MX.zeros(1,14);
        all_Qd2dot(1:7)=Qd2dots_prescribed_j(1:7,j);
        all_Qd2dot(8:14)=Qd2dots_j(:,j)'.*scaling.qd2dot;
        
        out=F([all_QsQdot';all_Qd2dot';forces_prescribed_j(:,j)]);
        
        out_d{i}(:,j)=out;
        %dynamic constraints
        FTtildep_nsc = [FTtilde_k FTtilde_j]*scaling.FTtilde*C(:,j+1); %FTtilde_sol already in the original scale
        eq_constr{end+1}=(h{i}*dFTtilde_j(:,j)*scaling.dFTtilde-FTtildep_nsc)/scaling.dFTtilde;
        g2_names=[g2_names; repmat({'FTtilde cons'},nmuscles,1)];

        Qs_nsc = [QsQdots_k(1:2:end,:) QsQdots_j(1:2:end,:)].*scaling.QsQdots(1:2:end)'*C(:,j+1); %QsQdots_sol already in the original scale
        eq_constr{end+1}=(h{i}*QsQdots_j(2:2:end,j).*scaling.QsQdots(2:2:end)'-Qs_nsc)/scaling.QsQdots(2);
        g2_names=[g2_names; repmat({'Qs cons'},ndofs,1)];

        Qdots_nsc = [QsQdots_k(2:2:end,:) QsQdots_j(2:2:end,:)].*scaling.QsQdots(2:2:end)'*C(:,j+1); %QsQdots_sol already in the original scale
        eq_constr{end+1}=(h{i}*Qd2dots_j(:,j).*scaling.qd2dot'-Qdots_nsc)/scaling.qd2dot(1);
        g2_names=[g2_names; repmat({'Qsdot cons'},ndofs,1)];

        %Path constraints
        % Hill equilibrium
        eq_constr{end+1}=Hilldiff_j{j};
        g2_names=[g2_names; repmat({'Hilldiff_eq'},nmuscles,1)];
              
        %Muscle force sharing
        %moment equilibrium
        if Options.dofs_to_track(1)
            FT_hip_flx=FT_j{j}(mai(1).mus);
            T_hip_flx=FT_hip_flx'*MAj{j}.hip_flex;
            if Options.optimizePassiveJointEl
                PassiveM_hip_flx=-Kstiff_k(1)*(QsQdots_j(1*2-1,j)*scaling.q-theta0_k(1))-Ddamp_k(1)*QsQdots_j(1*2,j)*scaling.qdot;
                eq_constr{end+1}=(out(8)-T_hip_flx-PassiveM_hip_flx)/scaling.T;
            else
                eq_constr{end+1}=(out(8)-T_hip_flx)/scaling.T;
            end
            g2_names=[g2_names; 'hip flex'];
        end
        if Options.dofs_to_track(2)
            %hip adduction
            FT_hip_add=FT_j{j}(mai(2).mus);
            T_hip_add=FT_hip_add'*MAj{j}.hip_add;
            if Options.optimizePassiveJointEl
                PassiveM_hip_add=-Kstiff_k(1)*(QsQdots_j(2*2-1,j)*scaling.q-theta0_k(2))-Ddamp_k(1)*QsQdots_j(2*2,j)*scaling.qdot;
                eq_constr{end+1}=(out(9)-T_hip_add-PassiveM_hip_add)/scaling.T;
            else
                eq_constr{end+1}=(out(9)-T_hip_add)/scaling.T;
            end
            g2_names=[g2_names; 'hip_add'];
        end
        if Options.dofs_to_track(3)
            %hip rotation
            FT_hip_int=FT_j{j}(mai(3).mus);
            T_hip_int=FT_hip_int'*MAj{j}.hip_int;
            if Options.optimizePassiveJointEl
                PassiveM_hip_int=-Kstiff_k(1)*(QsQdots_j(3*2-1,j)*scaling.q-theta0_k(3))-Ddamp_k(1)*QsQdots_j(3*2,j)*scaling.qdot;
                eq_constr{end+1}=(out(10)-T_hip_int-PassiveM_hip_int)/scaling.T;
            else
                eq_constr{end+1}=(out(10)-T_hip_int)/scaling.T;
            end
            g2_names=[g2_names; {'hip rot'}];
        end
        if Options.dofs_to_track(4)
            %knee flexion
            FT_knee_flx=FT_j{j}(mai(4).mus);
            T_knee_flx=FT_knee_flx'*MAj{j}.knee_flex;
            if Options.optimizePassiveJointEl
                PassiveM_knee_flx=-Kstiff_k(2)*(QsQdots_j(4*2-1,j)*scaling.q-theta0_k(4))-Ddamp_k(2)*QsQdots_j(4*2,j)*scaling.qdot;
                eq_constr{end+1}=(out(11)-T_knee_flx-PassiveM_knee_flx)/scaling.T;
            else
                eq_constr{end+1}=(out(11)-T_knee_flx)/scaling.T;
            end
            g2_names=[g2_names; {'knee flex'}];
        end
        if Options.dofs_to_track(5)
            %ankle flexion
            FT_ankle_flx=FT_j{j}(mai(5).mus);
            T_ankle_flx=FT_ankle_flx'*MAj{j}.ankle_flex;
            if Options.optimizePassiveJointEl
                PassiveM_ankle_flx=-Kstiff_k(3)*(QsQdots_j(5*2-1,j)*scaling.q-theta0_k(5))-Ddamp_k(3)*QsQdots_j(5*2,j)*scaling.qdot;
                eq_constr{end+1}=(out(12)-T_ankle_flx-PassiveM_ankle_flx)/scaling.T;
            else
                eq_constr{end+1}=(out(12)-T_ankle_flx)/scaling.T;
            end
            g2_names=[g2_names; {'ankle flex'}];
        end
        if Options.dofs_to_track(6)
            %ankle adduction
            FT_ankle_add=FT_j{j}(mai(6).mus);
            T_ankle_add=FT_ankle_add'*MAj{j}.ankle_add;
            if Options.optimizePassiveJointEl
                PassiveM_ankle_add=-Kstiff_k(3)*(QsQdots_j(6*2-1,j)*scaling.q-theta0_k(6))-Ddamp_k(3)*QsQdots_j(6*2,j)*scaling.qdot;
                eq_constr{end+1}=(out(13)-T_ankle_add-PassiveM_ankle_add)/scaling.T;
            else
                eq_constr{end+1}=(out(13)-T_ankle_add)/scaling.T;
            end
            g2_names=[g2_names; {'ankle add'}];
        end
        if Options.dofs_to_track(7)
            %ankle int
            FT_ankle_int=FT_j{j}(mai(7).mus);
            T_ankle_int=FT_ankle_int'*MAj{j}.ankle_int;
            if Options.optimizePassiveJointEl
                PassiveM_ankle_int=-Kstiff_k(3)*(QsQdots_j(7*2-1,j)*scaling.q-theta0_k(7))-Ddamp_k(3)*QsQdots_j(7*2,j)*scaling.qdot;
                eq_constr{end+1}=(out(14)-T_ankle_int-PassiveM_ankle_int)/scaling.T;
            else
                eq_constr{end+1}=(out(14)-T_ankle_int)/scaling.T;
            end
            g2_names=[g2_names; {'ankle int'}];
        end
        %Cost function
        J_d=J_d+B(j+1)*W.qtrack*sum((QsQdots_j(1:2:end,j)-QsQdots_totrack_j(1:2:end,j)).^2);
%         J_d=J_d+B(j+1)*W.qdottrack*sum((QsQdots_j(2:2:end,j)-QsQdots_totrack_j(2:2:end,j)).^2);
        J_d=J_d+B(j+1)*W.mindstate*sum(dFTtilde_j(:,j).^2);
        if Options.penalizeHighFTtilde
            J_d=J_d+B(j+1)*W.penalizeHighFTtilde*sum(FTtilde_j(:,j).^2);
        end
        if Options.penalizeFTtot
            J_d=J_d+B(j+1)*W.penalizeFTtot*sum(FT_j{j}.^2);
        end
        if Options.penalizeoutoflMtilde1
            J_d=J_d+B(j+1)*W.penalizeoutoflMtilde1*sum((lMtilde_j{j}-1).^2);
        end
        %for debug
        J1{i}=J1{i}+B(j+1)*W.qtrack*sum((QsQdots_j(1:2:end,j)-QsQdots_totrack_j(1:2:end,j)).^2); %track q
        J2{i}=J2{i}+B(j+1)*W.mindstate*sum(dFTtilde_j(:,j).^2);
        if Options.penalizeHighFTtilde
            J3{i}=J3{i}+B(j+1)*W.penalizeHighFTtilde*sum(FTtilde_j(:,j).^2);
        end
        if Options.penalizeFTtot
            J4{i}=J4{i}+B(j+1)*W.penalizeFTtot*sum(FT_j{j}.^2);
        end
        if Options.penalizeoutoflMtilde1
            J5{i}=J5{i}+B(j+1)*W.penalizeoutoflMtilde1*sum((lMtilde_j{j}-1).^2);
        end
    end
    eq_constr = vertcat(eq_constr{:});
    ineq_constr=vertcat(ineq_constr{:});
    
    if Options.useRigidTendon
        keyboard;
    else
        f_coll = Function('f_coll',{FTtilde_k,FTtilde_j,...
            dFTtilde_j, QsQdots_k,QsQdots_j,...
            QsQdots_totrack_j,QsQdots_prescribed_j, ...
            Qd2dots_j,Qd2dots_prescribed_j, forces_prescribed_j, ...
            lTs_k, lM0_k, fiber_damping_k,tendon_damping_k,Kstiff_k,Ddamp_k,theta0_k},{eq_constr,J_d,...
            [FT_j{1} FT_j{2} FT_j{3}],...
            [lMtilde_j{1} lMtilde_j{2} lMtilde_j{3}],...
            [lTtilde_j{1} lTtilde_j{2} lTtilde_j{3}]});
        
        f_coll_map = f_coll.map(N,ParallelMode,NThreads);

        [coll_eq_constr{i}, Jall{i},FT_all{i},...
            lMtilde_all{i},lTtilde_all{i}] = ...
            f_coll_map(FTtilde{i}(:,1:end-1),...
            FTtilde_col{i}, dFTtilde_col{i}, ...
            QsQdots{i}(:,1:end-1),QsQdots_col{i},...
            guess.QsQdots{i}(t_col_grid,:)',... %guess.QsQdots is already scaled (this is for QsQdots_totrack_j)
            QsQdot_prescribed{i}(t_col_grid,:)',... %QsQdot_prescribed is unscaled
            Qd2dots_col{i},...
            Qd2dot_prescribed{i}(t_col_grid,:)',forces_prescribed{i}',...
            repmat(lTs,1,N),repmat(lM0,1,N),repmat(fiber_damping,1,N),...
            repmat(tendon_damping,1,N),repmat(Kstiff,1,N),repmat(Ddamp,1,N),repmat(theta0,1,N));
        
        %for debugging
        f_coll2 = Function('f_coll2',{FTtilde_k,FTtilde_j,...
            dFTtilde_j, QsQdots_k,QsQdots_j,...
            QsQdots_totrack_j,QsQdots_prescribed_j, ...
            Qd2dots_j,Qd2dots_prescribed_j, forces_prescribed_j, ...
            lTs_k, lM0_k, fiber_damping_k,tendon_damping_k,Kstiff_k,Ddamp_k,theta0_k},{J1{i},J2{i},J3{i},J4{i},J5{i},out_d{i}});
        f_coll_map2=f_coll2.map(N,ParallelMode,NThreads);
        [J1all{i},J2all{i},J3all{i},J4all{i},J5all{i},out_dall{i}]=f_coll_map2(FTtilde{i}(:,1:end-1),...
            FTtilde_col{i}, dFTtilde_col{i}, ...
            QsQdots{i}(:,1:end-1),QsQdots_col{i},...
            guess.QsQdots{i}(t_col_grid,:)',... %guess.QsQdots is already scaled (this is for QsQdots_totrack_j)
            QsQdot_prescribed{i}(t_col_grid,:)',... %QsQdot_prescribed is unscaled
            Qd2dots_col{i},...
            Qd2dot_prescribed{i}(t_col_grid,:)',forces_prescribed{i}',...
            repmat(lTs,1,N),repmat(lM0,1,N),repmat(fiber_damping,1,N),...
            repmat(tendon_damping,1,N),repmat(Kstiff,1,N),repmat(Ddamp,1,N),repmat(theta0,1,N));
    end
    opti.subject_to(coll_eq_constr{i}==0);
     % Add continuity constraints (next interval starts with end values of 
    % states from previous interval)
    if Options.useRigidTendon
        %no need for continuity constriants when using rigid
        %tendon
    else
        for k=1:N
            % Variables within current mesh interval
            % States    
            FTtilde_kj = [FTtilde{i}(:,k), FTtilde_col{i}(:,(k-1)*d+1:k*d)];
            opti.subject_to(FTtilde{i}(:,k+1) == FTtilde_kj*D); % scaled
            g3_names=[g3_names; repmat({'continuity constr FTtilde'},nmuscles,1)];

            QsQdots_kj=[QsQdots{i}(:,k), QsQdots_col{i}(:,(k-1)*d+1:k*d)];
            opti.subject_to(QsQdots{i}(:,k+1)==QsQdots_kj*D); %scaled
            g3_names=[g3_names; repmat({'continuity constr QsQdots'},ndofs*2,1)];
        end
    end  
    J_i=J_i+sum(Jall{i});
    
end

g_names=[g1_names; repmat([g2_names_ineq; g2_names],N,1); g3_names];  

J=J+J_i;

%% Solve OCP
opti.minimize(J);
options.ipopt.hessian_approximation = 'limited-memory'; %'exact'; %
options.ipopt.mu_strategy      = 'adaptive';
options.ipopt.max_iter = 2000;
options.ipopt.tol = 1e-4;   
opti.solver('ipopt', options);  
sol=opti.solve();


if opti.stats.success==1
    sol_val.stats=sol.stats;
    if Options.optimizelM0
        sol_val.lM0=sol.value(lM0);
        sol_val.lM0_unsc=sol_val.lM0*scaling.lM0;
    end
    if Options.optimizelTs
        sol_val.lTs=sol.value(lTs);
        sol_val.lTs_unsc=sol_val.lTs*scaling.lTs;
    end
    
else
    sol_val.stats=opti.stats;
    if Options.optimizelM0
        sol_val.lM0=opti.debug.value(lM0);
        sol_val.lM0_unsc=sol_val.lM0*scaling.lM0;
    end
    if Options.optimizelTs
        sol_val.lTs=opti.debug.value(lTs);
        sol_val.lTs_unsc=sol_val.lTs*scaling.lTs;
    end
end
sol_val.Options=Options;
sol_val.scaling=scaling;
sol_val.ndofs=ndofs;

sol_val.tgrid=tgrid;
sol_val.t2plot=t_col_grid;
sol_val.W=W;
sol_val.muscle_names=muscle_name;
sol_val.stats=opti.stats;

if opti.stats.success
    for i=1:size(QsQdots,2)
        sol_val.QsQdots{i}=sol.value(QsQdots{i});
        sol_val.QsQdots_col{i}=sol.value(QsQdots_col{i});
        sol_val.Qd2dot_col{i}=sol.value(Qd2dots_col{i});
        if ~isempty(sol_val.QsQdots{i})
            sol_val.QsQdots_unsc{i}=sol_val.QsQdots{i}.*scaling.QsQdots';
            sol_val.QsQdots_col_unsc{i}=sol_val.QsQdots_col{i}.*scaling.QsQdots';
            sol_val.Qd2dot_col_unsc{i}=sol_val.Qd2dot_col{i}*scaling.qd2dot;
        end
    end

    for i=1:size(FTtilde,2)
            sol_val.FTtilde{i}=[];
            sol_val.dFTtilde{i}=[];
            for k=1:size(FTtilde{i},2)
                if k<size(FTtilde{i},2)
                    sol_val.FTtilde{i}(:,(k-1)*(d+1)+1:(k*(d+1)))=[sol.value(FTtilde{i}(:,k)) sol.value(FTtilde_col{i}(:,(k-1)*d+1:k*d))];
                    sol_val.dFTtilde{i}=sol.value(dFTtilde_col{i});
                else
                    sol_val.FTtilde{i}(:,(k-1)*(d+1)+1)=sol.value(FTtilde{i}(:,k));
                end
            end
            sol_val.FTtilde_unsc{i}=sol_val.FTtilde{i}*scaling.FTtilde;
            sol_val.dFTtilde_unsc{i}=sol_val.dFTtilde{i}*scaling.dFTtilde;
            sol_val.FT_all{i}=sol.value(FT_all{i});
            sol_val.lMtilde_all{i}=sol.value(lMtilde_all{i});
            sol_val.lTtilde_all{i}=sol.value(lTtilde_all{i});
    end
    if Options.optimizetendondamping
        sol_val.tendon_damping=sol.value(tendon_damping);
    else
        sol_val.tendon_damping=0;
    end
    if Options.optimizefiberdamping
        sol_val.fiber_damping=sol.value(fiber_damping);
    else
        sol_val.fiber_damping=0;
    end   
else      
    for i=1:size(QsQdots,2)
        sol_val.QsQdots{i}=opti.debug.value(QsQdots{i});
        sol_val.QsQdots_col{i}=opti.debug.value(QsQdots_col{i});
        sol_val.Qd2dot_col{i}=opti.debug.value(Qd2dots_col{i});
        if ~isempty(sol_val.QsQdots{i})
            sol_val.QsQdots_unsc{i}=sol_val.QsQdots{i}.*scaling.QsQdots';
            sol_val.QsQdots_col_unsc{i}=sol_val.QsQdots_col{i}.*scaling.QsQdots';
            sol_val.Qd2dot_col_unsc{i}=sol_val.Qd2dot_col{i}*scaling.qd2dot;
        end
    end
    if Options.useRigidTendon==1
    else
        for i=1:size(FTtilde,2)
            sol_val.FTtilde{i}=[];
            sol_val.dFTtilde{i}=[];
            for k=1:size(FTtilde{i},2)
                if k<size(FTtilde{i},2)
                    sol_val.FTtilde{i}(:,(k-1)*(d+1)+1:(k*(d+1)))=[opti.debug.value(FTtilde{i}(:,k)) opti.debug.value(FTtilde_col{i}(:,(k-1)*d+1:k*d))];
                    sol_val.dFTtilde{i}=opti.debug.value(dFTtilde_col{i});
                else
                    sol_val.FTtilde{i}(:,(k-1)*(d+1)+1)=opti.debug.value(FTtilde{i}(:,k));
                end
            end
            
            sol_val.FT_all{i}=opti.debug.value(FT_all{i});
            sol_val.lMtilde_all{i}=opti.debug.value(lMtilde_all{i});
            sol_val.lTtilde_all{i}=opti.debug.value(lTtilde_all{i});
            
            if ~isempty(FTtilde{i})
                sol_val.FTtilde_unsc{i}=sol_val.FTtilde{i}*scaling.FTtilde;
                sol_val.dFTtilde_unsc{i}=sol_val.dFTtilde{i}*scaling.dFTtilde;
            end

        end

        if Options.optimizetendondamping
            sol_val.tendon_damping=opti.debug.value(tendon_damping);
        else
            sol_val.tendon_damping=0;
        end
    end
    if Options.optimizefiberdamping
        sol_val.fiber_damping=opti.debug.value(fiber_damping);
    else
        sol_val.fiber_damping=0;
    end

end
sol_val.guess=guess;
sol_val.QsQdot_prescribed=QsQdot_prescribed;
sol_val.Qd2dot_prescribed=Qd2dot_prescribed;
sol_val.force_prescribed=forces_prescribed;

% Recompute lMT
for i=1:size(FTtilde,2)
    if ~isempty(FTtilde{i})
        for k=1:size(FTtilde{i},2)-1
            all_qsleg=sol_val.QsQdots_col_unsc{i}(1:2:end,(k-1)*d+1:k*d)';
            all_qdotsleg=sol_val.QsQdots_col_unsc{i}(2:2:end,(k-1)*d+1:k*d)';
            for j=1:d
                %Get moment arms and muscle-tendon lengths at that frame
                [lMTj_aux,vMTj_aux,MAj_aux] =  f_lMT_vMT_dM(all_qsleg(j,:),all_qdotsleg(j,:));
                [Hilldiff_j_aux{j},FT_j_aux{j},~,~,Fp_j_aux{j},lMtilde_j_aux{j},lTtilde_j_aux{j}]=f_forceEquilibrium_FtildeState(...
                                       zeros(nmuscles,1),sol_val.FTtilde_unsc{i}(:,(k-1)*(d+1)+j+1),...
                                       sol_val.dFTtilde_unsc{i}(:,(k-1)*d+j),...
                                       full(lMTj_aux),full(vMTj_aux),sol_val.lTs_unsc,sol_val.lM0_unsc,sol_val.fiber_damping,sol_val.tendon_damping);
                eq_Hill{i}((k-1)*d+j,:)=full(Hilldiff_j_aux{j});
                sol_val.lMT{i}((k-1)*d+j,:)=full(lMTj_aux);

                MA_opt.hip_flex((k-1)*d+j,:)   =  full(MAj_aux(mai(1).mus',1));
                MA_opt.hip_add((k-1)*d+j,:)    =  full(MAj_aux(mai(2).mus',2));
                MA_opt.hip_int((k-1)*d+j,:)    =  full(MAj_aux(mai(3).mus',3));
                MA_opt.knee_flex((k-1)*d+j,:)  =  full(MAj_aux(mai(4).mus',4));
                MA_opt.ankle_flex((k-1)*d+j,:) =  full(MAj_aux(mai(5).mus',5));  
                MA_opt.ankle_add((k-1)*d+j,:)  =  full(MAj_aux(mai(6).mus',6)); 
                MA_opt.ankle_int((k-1)*d+j,:)  =  full(MAj_aux(mai(7).mus',7)); 

            end
        end
    end
end

%plot cost function terms
for i=1:length(J1)
    figure
    plot(opti.value(J1all{i}));
    hold all
    plot(opti.value(J2all{i}));
    plot(opti.value(J3all{i}));
    plot(opti.value(J4all{i}));
    plot(opti.value(J5all{i}));
    legend({'J1','J2','J3','J4','J5'});
end
if Options.optimizelM0
    opti.value(J_lM0)
end
if Options.optimizePassiveJointEl
    opti.value(J_passjoint)
end

%% write output kinematics
q.data(:,1)=sol_val.tgrid{1}(sol_val.t2plot);
q.data(:,2:8)=sol_val.QsQdot_prescribed{1}(sol_val.t2plot,1:2:end);
q.data(:,9:15)=sol_val.QsQdots_col_unsc{1}(1:2:end,:)';
q.data(:,[2:4 8:end])=q.data(:,[2:4 8:end])*180/pi;
q.labels={'time', 'sacrum_pitch','sacrum_roll','sacrum_yaw','sacrum_x','sacrum_y','sacrum_z','sacroiliac_flx', ...
    'hip_flx','hip_add','hip_int','knee_flx','ankle_flx','ankle_add','ankle_int'};
write_motionFile(q,'testoutmot.mot');

%% recompute constraint values
for i=1:length(sol_val.QsQdot_prescribed)
    if ~isempty(sol_val.QsQdot_prescribed{i})
        all_QsQdot_opt{i}(:,1:14)=sol_val.QsQdot_prescribed{i}(sol_val.t2plot,:);
        all_QsQdot_opt{i}(:,15:28)=sol_val.QsQdots_col_unsc{i}';

        all_Qd2dot_opt{i}(:,1:7)=sol_val.Qd2dot_prescribed{i}(sol_val.t2plot,:);
        all_Qd2dot_opt{i}(:,8:14)=sol_val.Qd2dot_col_unsc{i}';

        for k=1:size(all_Qd2dot_opt{i},1)
            out_opt{i}(k,:)=full(F([all_QsQdot_opt{i}(k,:)';all_Qd2dot_opt{i}(k,:)';forces_prescribed{i}(k,:)']));
        end
    end
end
%Muscle force sharing
%moment equilibrium
for i=1:length(sol_val.FT_all)
    if ~isempty(sol_val.FT_all{i})
        if Options.dofs_to_track(1)
            FT_hip_flx_opt=sol_val.FT_all{i}(mai(1).mus,:); 
            T_hip_flx_opt=FT_hip_flx_opt'.*MA_opt.hip_flex;
            if Options.optimizePassiveJointEl
                PassiveM_hip_flx_opt=-sol_val.Kstiff(1)*(sol_val.QsQdots_col_unsc{i}(1,:)'-sol_val.theta0(1))-sol_val.Ddamp(1)*sol_val.QsQdots_col_unsc{i}(2,:)';
                eq_constr_opt{i}(:,1)=(out_opt{i}(:,8)-sum(T_hip_flx_opt,2)-PassiveM_hip_flx_opt)/scaling.T;
            else
                eq_constr_opt{i}(:,1)=(out_opt{i}(:,8)-sum(T_hip_flx_opt,2))/scaling.T;
            end

        end
        if Options.dofs_to_track(2)
            %hip adduction
            FT_hip_add_opt=sol_val.FT_all{i}(mai(2).mus,:);
            T_hip_add_opt=FT_hip_add_opt'.*MA_opt.hip_add;
            if Options.optimizePassiveJointEl
                PassiveM_hip_add_opt=-sol_val.Kstiff(1)*(sol_val.QsQdots_col_unsc{i}(3,:)'-sol_val.theta0(2))-sol_val.Ddamp(1)*sol_val.QsQdots_col_unsc{i}(4,:)';
                eq_constr_opt{i}(:,2)=(out_opt{i}(:,9)-sum(T_hip_add_opt,2)-PassiveM_hip_add_opt)/scaling.T;
            else
                eq_constr_opt{i}(:,2)=(out_opt{i}(:,9)-sum(T_hip_add_opt,2))/scaling.T;
            end
        end
        if Options.dofs_to_track(3)
            %hip rotation
            FT_hip_int_opt=sol_val.FT_all{i}(mai(3).mus,:);
            T_hip_int_opt=FT_hip_int_opt'.*MA_opt.hip_int;
            if Options.optimizePassiveJointEl
                PassiveM_hip_int_opt=-sol_val.Kstiff(1)*(sol_val.QsQdots_col_unsc{i}(5,:)'-sol_val.theta0(3))-sol_val.Ddamp(1)*sol_val.QsQdots_col_unsc{i}(6,:)';
                eq_constr_opt{i}(:,3)=(out_opt{i}(:,10)-sum(T_hip_int_opt,2)-PassiveM_hip_int_opt)/scaling.T;
            else
                eq_constr_opt{i}(:,3)=(out_opt{i}(:,10)-sum(T_hip_int_opt,2))/scaling.T;
            end
        end
        if Options.dofs_to_track(4)
            %knee flexion
            FT_knee_flx_opt=sol_val.FT_all{i}(mai(4).mus,:);
            T_knee_flx_opt=FT_knee_flx_opt'.*MA_opt.knee_flex;
            if Options.optimizePassiveJointEl
                PassiveM_knee_flx_opt=-sol_val.Kstiff(2)*(sol_val.QsQdots_col_unsc{i}(7,:)'-sol_val.theta0(4))-sol_val.Ddamp(2)*sol_val.QsQdots_col_unsc{i}(8,:)';
                eq_constr_opt{i}(:,4)=(out_opt{i}(:,11)-sum(T_knee_flx_opt,2)-PassiveM_knee_flx_opt)/scaling.T;
            else
                eq_constr_opt{i}(:,4)=(out_opt{i}(:,11)-sum(T_knee_flx_opt,2))/scaling.T;
            end
        end
        if Options.dofs_to_track(5)
            %ankle flexion
            FT_ankle_flx_opt=sol_val.FT_all{i}(mai(5).mus,:);
            T_ankle_flx_opt=FT_ankle_flx_opt'.*MA_opt.ankle_flex;
            if Options.optimizePassiveJointEl
                PassiveM_ankle_flx_opt=-sol_val.Kstiff(3)*(sol_val.QsQdots_col_unsc{i}(9,:)'-sol_val.theta0(5))-sol_val.Ddamp(3)*sol_val.QsQdots_col_unsc(10,:)';
                eq_constr_opt{i}(:,5)=(out_opt{i}(:,12)-sum(T_ankle_flx_opt,2)-PassiveM_ankle_flx_opt)/scaling.T;
            else
                eq_constr_opt{i}(:,5)=(out_opt{i}(:,12)-sum(T_ankle_flx_opt,2))/scaling.T;
            end
        end
        if Options.dofs_to_track(6)
            %ankle adduction
            FT_ankle_add_opt=sol_val.FT_all{i}(mai(6).mus,:);
            T_ankle_add_opt=FT_ankle_add_opt'.*MA_opt.ankle_add;
            if Options.optimizePassiveJointEl
                PassiveM_ankle_add_opt=-sol_val.Kstiff(3)*(sol_val.QsQdots_col_unsc{i}(11,:)'-sol_val.theta0(6))-sol_val.Ddamp(3)*sol_val.QsQdots_col_unsc(12,:)';
                eq_constr_opt{i}(:,6)=(out_opt{i}(:,13)-sum(T_ankle_add_opt,2)-PassiveM_ankle_add_opt)/scaling.T;
            else
                eq_constr_opt{i}(:,6)=(out_opt{i}(:,13)-sum(T_ankle_add_opt,2))/scaling.T;
            end
        end
        if Options.dofs_to_track(7)
            %ankle int
            FT_ankle_int_opt=sol_val.FT_all{i}(mai(7).mus,:);
            T_ankle_int_opt=FT_ankle_int_opt'.*MA_opt.ankle_int;
            if Options.optimizePassiveJointEl
                PassiveM_ankle_int_opt=-sol_val.Kstiff(3)*(sol_val.QsQdots_col_unsc{i}(13,:)'-sol_val.theta0(7))-sol_val.Ddamp(3)*sol_val.QsQdots_col_unsc(14,:)';
                eq_constr_opt{i}(:,7)=(out_opt{i}(:,14)-sum(T_ankle_int_opt,2)-PassiveM_ankle_int_opt)/scaling.T;
            else
                eq_constr_opt{i}(:,7)=(out_opt{i}(:,14)-sum(T_ankle_int_opt,2))/scaling.T;
            end
        end
    end
end

%dynamic constraints
for i=1:length(sol_val.FTtilde_unsc)
    if ~isempty(sol_val.FTtilde_unsc{i})
        for k=1:N
            for j=1:d
                FTtildep_opt_nsc{i}(k,:) = [sol_val.FTtilde_unsc{i}(:,(k-1)*(d+1)+1:k*(d+1))]*C(:,j+1); %FTtilde_sol already in the original scale
                eq_constr_dyn_FTtilde{i}(k,:)=(h{i}*sol_val.dFTtilde_unsc{i}(:,(k-1)*d+j)-FTtildep_opt_nsc{i}(k,:)')/scaling.dFTtilde;

                Qs_opt_nsc{i}(k,:) = [sol_val.QsQdots_unsc{i}(1:2:end,k) sol_val.QsQdots_col_unsc{i}(1:2:end,(k-1)*d+1:(k-1)*d+d)]*C(:,j+1); %QsQdots_sol already in the original scale
                eq_constr_dyn_Qsopt{i}(k,:)=(h{i}*sol_val.QsQdots_col_unsc{i}(2:2:end,(k-1)*d+j)-Qs_opt_nsc{i}(k,:)')/scaling.QsQdots(2);

                Qdots_opt_nsc{i}(k,:) = [sol_val.QsQdots_unsc{i}(2:2:end,k) sol_val.QsQdots_col_unsc{i}(2:2:end,(k-1)*d+1:(k-1)*d+d)]*C(:,j+1); %QsQdots_sol already in the original scale
                eq_constr_dyn_Qdot{i}(k,:)=(h{i}*sol_val.Qd2dot_col_unsc{i}(:,(k-1)*d+j)-Qdots_opt_nsc{i}(k,:)')/scaling.qd2dot(1);
            end
        end
    end
end


function  expdata=LoadData(N,d,tau_root,trials)
    current_folder=pwd;
    for triali=1:length(trials)
        kinfiles=dir(['DataMay/kinematics/' trials{triali} '/*.mot']);
        forcefiles=dir(['DataMay/motion&force/' trials{triali} '/*.mot']);
        for i=1:length(kinfiles);
           kinfilename=kinfiles(i).name; 
           kindata=importdata([kinfiles(i).folder '/' kinfiles(i).name]);
           kindata=ProcessKinematics(kindata);
           C=strrep(kinfilename,'perturb','');
           sufix=strrep(C,'.mot','');
            
           trial_name=[trials{triali} '_' strrep(kinfilename,'.mot','')];
           %parameterize with splines
           t=kindata.data(1,1):0.002:kindata.data(end,1);
           expdata.(trial_name).kinematics(:,1)=t;
           expdata.(trial_name).kinematics_v(:,1)=t;
           expdata.(trial_name).kinematics_a(:,1)=t;
           for j=2:size(kindata.data,2)
                intdata=interp1(kindata.data(:,1),kindata.data(:,j),t,'spline');
                smoothed_kin=smooth(t,intdata,0.3,'rloess');

                kindata_spline(j-1)=spline(t,smoothed_kin);
                expdata.(trial_name).kinematics(:,j)=ppval(kindata_spline(j-1),t);
                expdata.(trial_name).kinematics_v(:,j)=ppval(fnder(kindata_spline(j-1),1),t);
                expdata.(trial_name).kinematics_a(:,j)=ppval(fnder(kindata_spline(j-1),2),t);
           end
           expdata.(trial_name).kinematics_labels=kindata.colheaders;
%            forcefilename2match=['high_posture' sufix '_ForwardResp_at_'];
%            found=false;
%            j=1;
%            while (~found)&(j<=length(forcefiles))
%                found=contains(forcefiles(j).name,forcefilename2match);
%                if found
%                else
%                 j=j+1;
%                end
%            end
           if strcmp(forcefiles(i).name,kinfiles(i).name)
           else
               keyboard;
           end
           forcedata=importdata([forcefiles(i).folder '/' forcefiles(i).name]);
           expdata.(trial_name).forces=forcedata.data;
           expdata.(trial_name).forces_labels=forcedata.colheaders;

           t0=expdata.(trial_name).kinematics(1,1);
           tf=expdata.(trial_name).kinematics(end,1);
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
           expdata.(trial_name).f=interp1(forcedata.data(:,1),forcedata.data,tgrid);
        end


    end
    
    
end
function kindata_out=ProcessKinematics_old(kindata);
%to avoid outliers
    nframes=size(kindata.data,1);
    v=diff(kindata.data(:,2:end))./diff(kindata.data(:,1));
    [I,J]=find(abs(v)>600);
    kindata_out=kindata;
    finished=false;
    I_done=[];
    J_done=[];
    while ~finished
        clear vals;
        for i=1:length(I)
            vals(i,:)=abs(v(I(i),J(i)));
        end
        [vals_sort,II]=sort(vals,'descend');
        [I_done,J_done,II]=IJ_done(I,J,I_done,J_done,II);
        kindata_out.data(I(II(1)),J(II(1))+1)=interp1(kindata.data([I(II(1))-1 I(II(1))+1],1),kindata_out.data([I(II(1))-1 I(II(1))+1],J(II(1))+1),kindata.data(I(II(1)),1),'spline');
        v=diff(kindata_out.data(:,2:end))./diff(kindata_out.data(:,1));
        [I,J]=find(abs(v)>600);
        finished=length(II)==1;
        length(II)
    end
    
end

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
function [I_done,J_done,II]=IJ_done(I,J,I_done,J_done,II)

isthis=false;
i=1;
while i<=length(II)&(~isthis)
    
    if ~isempty(I_done)
        found=any(sum([I(II(i)) J(II(i))]==[I_done J_done],2)==2);
        if found
            II(i)=[];
        else
            isthis=true;
            I_done=[I_done; I(II(i))];
            J_done=[J_done; J(II(i))];
        end
    else
        isthis=true;
        I_done=[I_done; I(II(i))];
        J_done=[J_done; J(II(i))];
    end
end

        
end

function lMT_guess_at0=ComputelMTguess(q,qdot,f_lMT_vMT_dM,Options)
    %Get moment arms and muscle-tendon lengths at that frame
    for i=1:size(q,1)
        [lMTj,vMTj,MAj_aux]=f_lMT_vMT_dM(q(i,9:15),qdot(i,9:15));
        lMT(i,:)=full(lMTj);
    end
    lMT_guess_at0=mean(lMT,1);

end
clear all;
import casadi.*
method='legendre';

%% Options
Options.useRigidTendon=0;
Options.optimizelM0=1;
Options.optimizelTs=1;
Options.optimizefiberdamping=1;
Options.optimizetendondamping=0;
Options.dofs_to_track=[1 0 0 1 1 0 0]; %1. hip flex 2. hip add 3. hip int 4. knee flex 5. ankle flex 6. ankle add 7. ankle int
Options.useOptimizedIG=0;
Options.trialstotrack='all'; %if 'all', all trials are taken
Options.penalizeHighFTtilde=1;
Options.penalizeFTtot=0;
Options.penalizeoutoflMtilde1=0;
Options.optInertiaParam=1;
Options.trackqd2dot=1;
Options.minPassiveProp=1;
Options.pendevPassive=1; 
Options.startfromprevioussolution=0;
    Options.prev_sol=load('sol_val_optKinematics_optJointPassive_optInertia_DataMarch2025AchillescutForward&Backward_tolem5_WRONGwithoutcontconstr.mat');
Options.optInertiapassiveParam=1;

Options.optimizeMuscleProp=0;
Options.optimizePassiveJointEl=1;
Options.individualpassiveprop=0;

Options.secondfolder=1;
main_folder='DataMay2025\baseline_ForwardOnly_rel2Dangles\';
main_folder2='DataMay2025\baseline_BackwardOnly_rel2Dangles\';

Options.main_folder=main_folder;
%% 
N=40;
Options.N=N;
d=3;
[tau_root,C,D,B] = CollocationScheme(d,method);


% trials={'h8'}; main_folder='DataMay';
% trials='all';main_folder='DataSeptember';
% trials={'data6'};
expdata=LoadData(N,d,tau_root,main_folder);

if Options.secondfolder==1;
    expdata2=LoadData(N,d,tau_root,main_folder2);

    fieldnames2=fieldnames(expdata2);
    for i=1:length(fieldnames2)
        expdata.([fieldnames2{i} '_back'])=expdata2.([fieldnames2{i}]);
    end
end

current_folder=pwd;
if Options.optimizeMuscleProp
    nmuscles=38;
    muscle_names={'BFa','BFp','CF','STp','STa','SM','QF','TFL','GMa','GMe',...
        'GMi','Pir','GP','GA','AL','AM','AB','Pec','IP','OE','OI','GS','GI',...
        'RF','VL','VI','VM','MG','LG','Pla','Sol','TA','EDL','TP','FDL','FHL',...
        'Per','Pop'};
    muscles_present=ones(1,nmuscles);
    if contains(main_folder,'achillescut')
        Imuscles_cut=28:31;
        muscles_present(Imuscles_cut)=0;

        nmuscles_cut=4;
    else
        nmuscles_cut=0;
        Imuscles_cut=[];
    end
    
    %% Create CasADi function for Hill-muscle relations
    % Function for Hill-equilibrium
    if Options.useRigidTendon
        FM          = SX(nmuscles-nmuscles_cut,1); %muscle fiber forces
        lMtilde     = SX.sym('lMtilde',nmuscles-nmuscles_cut); % Normalized fiber lengths
        FMactFL     = SX.sym('FMactFL',nmuscles-nmuscles_cut); % Normalized force-length term
        FMactFV     = SX.sym('FMactFV',nmuscles-nmuscles_cut); % Normalized force-velocity term
    else
        FTtilde     = SX.sym('FTtilde',nmuscles-nmuscles_cut); % Normalized tendon forces
        dFTtilde    = SX.sym('dFTtilde',nmuscles-nmuscles_cut); % Time derivative tendon forces
        % tension_SX  = SX.sym('tension',NMuscle); % Tensions not used here
        vMmax       = SX(nmuscles-nmuscles_cut,1); % Maximum contraction velocities
        Fpetilde    = SX(nmuscles-nmuscles_cut,1); % Normalized passive forces
        lMtilde     = SX(nmuscles-nmuscles_cut,1); % Normalized fiber lengths
    
    end
    a           = SX.sym('a',nmuscles-nmuscles_cut); % Muscle activations
    lMT         = SX.sym('lMT',nmuscles-nmuscles_cut); % Muscle-tendon lengths
    vMT         = SX.sym('vMT',nmuscles-nmuscles_cut); % Muscle-tendon velocities
    Hilldiff    = SX(nmuscles-nmuscles_cut,1); % Hill-equilibrium   
    Fce         = SX(nmuscles-nmuscles_cut,1); % Contractile element forces
    Fiso        = SX(nmuscles-nmuscles_cut,1); % Normalized forces from force-length curve 
    FT          = SX(nmuscles-nmuscles_cut,1); % Tendon forces    
    lTtilde     = SX(nmuscles-nmuscles_cut,1); % Normalized tendon lengths
    
    % Parameters of force-length-velocity curves
    load Fvparam
    load Fpparam
    load Faparam
    load MTparam %Load parameter values from osim model (FM0, lM0, alphao, vMax)
    
    MTparam(:,Imuscles_cut)=[];

    lTs     = SX.sym('lTs',1,nmuscles-nmuscles_cut); %tendon slack length
    if Options.optimizelM0
        lM0     = SX.sym('lM0',1,nmuscles-nmuscles_cut); %optimal fiber length
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
        for m = 1:nmuscles-nmuscles_cut
            [Hilldiff(m),FT(m),Fce(m),Fiso(m),vMmax(m),Fpetilde(m),lMtilde(m),lTtilde(m)] = ...
                ForceEquilibrium_FtildeState(a(m),FTtilde(m),dFTtilde(m),...
                lMT(m),vMT(m),MTparam(1,m),lM0(m),lTs(m),MTparam(4,m),MTparam(5,m),Fvparam,Fpparam,Faparam,fiber_damping,tendon_damping);
                %tension_SX(m));
        end
    end
end

cd(current_folder);
if Options.optimizeMuscleProp

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
end

%% Load external function
if Options.optInertiaParam
    F = external('F','RightRatHindlimb_Zhong_InertiaVar.dll')
else
    F = external('F','RightRatHindlimb_Zhong.dll');   
end

ndofs=7; %ndofs to construct polynomials must be 7 here
if Options.optimizeMuscleProp
    %% Prepare polynomials
    % Indices of the muscles actuating the different joints for later use
    pathpolynomial = [pwd,'/Polynomials'];
    addpath(genpath(pathpolynomial));
    load([pathpolynomial,'/muscle_spanning_joint_INFO_subject.mat']);
    load([pathpolynomial,'/MuscleInfo_subject.mat']);
    pathmusclefunctions=[pwd,'\MuscleModel'];
    rmpath(genpath(pathmusclefunctions));
    addpath(genpath(pathmusclefunctions),'-begin');
    if contains(main_folder,'achillescut')
        muscle_names(Imuscles_cut)=[];
        muscle_spanning_joint_INFO(Imuscles_cut,:)=[];
        MuscleInfo.muscle(Imuscles_cut)=[];
    end
    [~,mai] = MomentArmIndices_3D(muscle_names,...
        muscle_spanning_joint_INFO);
    pathcasadi_functions=[pwd,'/VariousFunctions/'];
    addpath(pathcasadi_functions);
    NMuscle_pol=nmuscles-nmuscles_cut;
    CasADi_functions;
end


%% Formulation of optimal control problem
% ndofs=sum(Options.dofs_to_track); %3 hip, 1 knee 3 ankle
% dofs_in_ID=[9:15]; %hip flx, hip add, hip int, knee flex, ankle flex, ankle add, ankle int
name_dofs={'hip flx', 'hip add', 'hip int', 'knee flex', 'ankle flex', 'ankle add', 'ankle int'};

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
scaling.T=1;
scaling.Kstiff=0.1;
scaling.Ddamp=0.1;
scaling.theta=0.1;
scaling.inertiaparam=repmat([0.05 1e-5 0.05],1,3);
scaling.Inertiapassparam=0.01;

if Options.useOptimizedIG
end

W.lM0lit=0.00001;
W.min_maxa=0.1;
W.state_for_reg=1e-6;
W.mindstate=0.01; %0.1
W.qtrack=100; %10
W.qdottrack=10; %0.1
W.penalizeHighFTtilde=10;
W.Kstiff=1e-5; %minimization
W.Ddamp=1e-5;
W.theta0=1e-5;
W.Inertiapassparam=1e-5;
W.Kstiffpen=1e-3; %penalization from mean values
W.Ddamppen=1e-3;
W.theta0pen=1e-3;
W.Inertiapassparampen=1e-3;
W.penalizeFTtot=1e-4;
W.penalizeoutoflMtilde1=0;
W.InertiaParam=1e-2;

ParallelMode='thread';
NThreads=8;

%% Define bounds
if Options.optimizeMuscleProp
    bounds.FTtilde.lower=0*ones(nmuscles-nmuscles_cut,1);
    bounds.FTtilde.upper=5*ones(nmuscles-nmuscles_cut,1)/scaling.FTtilde;
    load('lMTmax.mat');
    lMTmax=lMTmax(find(muscles_present));
    bounds.lTs.lower=(1e-4)*ones(nmuscles-nmuscles_cut,1)/scaling.lTs;
    bounds.lTs.upper=lMTmax'/scaling.lTs;
    bounds.lM0.lower=(lMTmax'/4)/scaling.lM0;
    bounds.lM0.upper=ones(nmuscles-nmuscles_cut,1)*0.06/scaling.lM0;
    bounds.dFTtilde.lower = -1*ones(nmuscles-nmuscles_cut,1);
    bounds.dFTtilde.upper = 1*ones(nmuscles-nmuscles_cut,1);
    bounds.fiberdamping.lower=1e-3;
    bounds.fiberdamping.upper=10;
    bounds.tendondamping.lower=1e-3;
    bounds.tendondamping.upper=200;
end
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
bounds.K.lower=1e-5/scaling.Kstiff;
bounds.K.upper=10/scaling.Kstiff;
bounds.D.lower=1e-5/scaling.Ddamp;
bounds.D.upper=10/scaling.Ddamp;
bounds.Inertiapassparam.lower=0;
bounds.Inertiapassparam.upper=0.1/scaling.Inertiapassparam;
bounds_theta0.lower=[[0.1174,-0.2571,0.1223,-0.9237,0.2676,0.3305,0.0909]'-180*pi/180]/scaling.theta;
bounds.theta0.lower=bounds_theta0.lower(find(Options.dofs_to_track));
bounds_theta0.upper=[[0.1174,-0.2571,0.1223,-0.9237,0.2676,0.3305,0.0909]'+180*pi/180]/scaling.theta;
bounds.theta0.upper=bounds_theta0.upper(find(Options.dofs_to_track));
bounds.inertiaParam.lower=[0.005 1.e-7 -0.020 0.002 1e-8 0.010 0.0005 1e-9 -0.001]./scaling.inertiaparam; %mfem Izfem yfem mtib Iztib ytib mfoot Izfoot yfoot 
bounds.inertiaParam.upper=[0.025 1.e-5 -0.005 0.015 1e-6 0.020 0.0050 1e-7      0]./scaling.inertiaparam; %mfem Izfem yfem mtib Iztib ytib mfoot Izfoot yfoot 

%% Define guesses
if Options.optimizeMuscleProp
    if Options.startfromprevioussolution
        keyboard;
    end
    guess.FTtilde=zeros(nmuscles-nmuscles_cut,1)/scaling.FTtilde;
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
    guess.vMtilde=zeros(nmuscles-nmuscles_cut,1)/scaling.vMtilde;
end
if Options.optimizePassiveJointEl
    guess.Kstiff=1e-5*ones(3,1)/scaling.Kstiff;
    guess.Ddamp=1e-5*ones(3,1)/scaling.Ddamp;
    guess.Inertiapassparam=zeros(3,1)/scaling.Inertiapassparam;
end
% guess.theta0=[0.0135, -0.1651,-0.0637,-1.1571,-0.1768]'; % computed as the mean of all initial values, TO BE REcalculated
guess_theta0=[0.1174,-0.2571,0.1223,-0.9237,0.2676,0.3305,0.0909]'/scaling.theta; % computed as the mean of all initial values for h8 perturb1
guess.theta0=guess_theta0(find(Options.dofs_to_track)); %take only the dofs that are tracked

if Options.startfromprevioussolution
    if any(contains(fieldnames(Options.prev_sol.sol_val),'inertiaParam_unsc'))
        guess.inertiaParam=Options.prev_sol.sol_val.inertiaParam_unsc'./scaling.inertiaparam;
    else
        guess.inertiaParam=Options.prev_sol.sol_val.inertiaParam'./scaling.inertiaparam;
    end
else
    guess.inertiaParam=[0.01351 1.086e-06 -0.014936 0.00538 8.204e-07 0.0152275 0.00193 3.727e-08 -0.00546174]./scaling.inertiaparam;
end

%populate names of the constraints
g1_names=[];
g2_names=[];
g2_names_ineq=[];
g3_names=[];

%% Start with an empty optimization
opti = casadi.Opti();
if Options.optimizeMuscleProp
    if Options.optimizelTs
        lTs=opti.variable(nmuscles-nmuscles_cut,1);
        opti.subject_to(bounds.lTs.lower < lTs < bounds.lTs.upper);
        opti.set_initial(lTs,guess.lTs);
        g1_names=[g1_names; repmat({'lts_bounds'},nmuscles-nmuscles_cut,1)];
        lTs_k=MX.sym('lTs_k',nmuscles-nmuscles_cut,1);
    else 
        lTs     = MTparam(3,:);
        lTs=lTs/scaling.lTs;
    end
    if Options.optimizelM0
        lM0=opti.variable(nmuscles-nmuscles_cut,1);
        opti.subject_to(bounds.lM0.lower < lM0 < bounds.lM0.upper);
        opti.set_initial(lM0,guess.lM0);
        g1_names=[g1_names; repmat({'lM0_bounds'},nmuscles-nmuscles_cut,1)];
        lM0_k=MX.sym('lM0_k',nmuscles-nmuscles_cut,1);
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
end

J=0;

if Options.optInertiaParam
    inertiaParam=opti.variable(9,1);
    opti.subject_to(bounds.inertiaParam.lower'<=inertiaParam<=bounds.inertiaParam.upper');
    opti.set_initial(inertiaParam,guess.inertiaParam');
    g1_names=[g1_names; repmat({'inertiaParam_bounds'},9,1)];
    inertiaParam_k=MX.sym('inertiaParam_k',9,1);

    J=J+W.InertiaParam*sum(...
        (inertiaParam(1:3:end)-guess.inertiaParam(1:3:end)').^2+...
        ((inertiaParam(2:3:end)-guess.inertiaParam(2:3:end)')).^2+...
        (inertiaParam(3:3:end)-guess.inertiaParam(3:3:end)').^2);
    J_inertiaParam=W.InertiaParam*sum(...
        (inertiaParam(1:3:end)-guess.inertiaParam(1:3:end)').^2+...
        ((inertiaParam(2:3:end)-guess.inertiaParam(2:3:end)')).^2+...
        (inertiaParam(3:3:end)-guess.inertiaParam(3:3:end)').^2);
end

if Options.optimizeMuscleProp&Options.optimizelM0
    J=J+W.lM0lit*sum((((lM0*scaling.lM0)-MTparam(2,:)')./MTparam(2,:)').^2);
    J_lM0=W.lM0lit*sum((((lM0*scaling.lM0)-MTparam(2,:)')./MTparam(2,:)').^2);
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

if Options.optimizePassiveJointEl
    if Options.individualpassiveprop==1
        ntrials4passprop=size(fieldnames(expdata),1);
    elseif contains(main_folder,'May2025')
        ntrials4passprop=4;
        list4passprop=[4 3 2 1 1 2 3 4 4 3 2 1 1 2 3 4]; %perturbations at 30 20 10 0 0 10...
    else
        %define how many different passive prop are considered
        keyboard;
    end
    if Options.startfromprevioussolution
        %provided that the perturbations are at the same locations
        keyboard;
        if any(any(contains(fieldnames(Options.prev_sol.sol_val),'Kstiff_unsc')))
            guess.Kstiff=Options.prev_sol.sol_val.Kstiff_unsc{i}/scaling.Kstiff;
            guess.Ddamp=Options.prev_sol.sol_val.Ddamp_unsc{i}/scaling.Ddamp;
            guess.theta0=Options.prev_sol.sol_val.theta0_unsc{i}/scaling.theta;
        else
            guess.Kstiff=Options.prev_sol.sol_val.Kstiff{i}/scaling.Kstiff;
            guess.Ddamp=Options.prev_sol.sol_val.Ddamp{i}/scaling.Ddamp;
            guess.theta0=Options.prev_sol.sol_val.theta0{i}/scaling.theta;
        end

    end
    Kstiff=opti.variable(3,ntrials4passprop);
    opti.subject_to(bounds.K.lower<= Kstiff <= bounds.K.upper);
    opti.set_initial(Kstiff,repmat(guess.Kstiff,1,ntrials4passprop));
    g1_names=[g1_names; repmat({'Kstiff_bounds'},3*ntrials4passprop,1)];
    Ddamp=opti.variable(3,ntrials4passprop);
    opti.subject_to(bounds.D.lower<=Ddamp <= bounds.D.upper);
    opti.set_initial(Ddamp,repmat(guess.Ddamp,1,ntrials4passprop));
    g1_names=[g1_names; repmat({'Ddamp_bounds'},3*ntrials4passprop,1)];
    theta0=opti.variable(sum(Options.dofs_to_track),ntrials4passprop);
    opti.subject_to(bounds.theta0.lower<= theta0 <= bounds.theta0.upper);
    opti.set_initial(theta0,repmat(guess.theta0,1,ntrials4passprop));
    g1_names=[g1_names; repmat({'theta0_bounds'},length(theta0(:)),1)];
    Kstiff_k=MX.sym('Kstiff_k',3,ntrials4passprop);
    Ddamp_k=MX.sym('Ddamp_k',3,ntrials4passprop);
    theta0_k=MX.sym('theta0_k',sum(Options.dofs_to_track),ntrials4passprop); %only needed for orderPassiveJoint==1 and for collmap when the order is 1
    if Options.optInertiapassiveParam
        InertiapassiveParam=opti.variable(3,ntrials4passprop);
        opti.subject_to(bounds.Inertiapassparam.lower<= InertiapassiveParam <= bounds.Inertiapassparam.upper);
        opti.set_initial(InertiapassiveParam,repmat(guess.Inertiapassparam,1,ntrials4passprop));
        g1_names=[g1_names; repmat({'Inertiapassiveparam_bounds'},length(InertiapassiveParam(:)),1)];
        InertiapassiveParam_k=MX.sym('InertiapassiveParam_k',3,ntrials4passprop);
    else
        InertiapassiveParam_k=zeros(3,ntrials4passprop);
    end
    
else
    Kstiff_k=zeros(3,ntrials4passprop);
    Ddamp_k=zeros(3,ntrials4passprop);
    theta0_k=zeros(sum(Options.dofs_to_track),ntrials4passprop);
    InertiapassiveParam_k=zeros(3,ntrials4passprop);
end
Kstiff_k_unsc=Kstiff_k*scaling.Kstiff;
Ddamp_k_unsc=Ddamp_k*scaling.Ddamp;
theta0_k_unsc=theta0_k*scaling.theta;
InertiapassiveParam_k_unsc=InertiapassiveParam_k*scaling.Inertiapassparam;

if Options.optimizePassiveJointEl&(Options.minPassiveProp==1)
    J=J+W.Kstiff*sum(Kstiff(:).^2);
    J=J+W.Ddamp*sum(Ddamp(:).^2);
    J_passjoint=W.Kstiff*sum(Kstiff(:).^2)+W.Ddamp*sum(Ddamp(:).^2);
    if Options.optInertiapassiveParam
        J=J+W.Inertiapassparam*sum(InertiapassiveParam(:).^2);
        J_passjoint=W.Inertiapassparam*sum(InertiapassiveParam(:).^2);
    end
    if Options.pendevPassive
        %do not include twice the penalization of theta0
    else
        difftheta0=theta0-repmat(guess.theta0,1,ntrials4passprop);
        J=J+W.theta0*sum(difftheta0(:).^2);
        J_passjoint=J_passjoint+W.theta0*sum(difftheta0(:).^2);
    end
end

% for i=1:length(nametrials)
% for i=[1 3 6:2:12]
for i=1:length(nametrials)
    eq_constr={};
    ineq_constr={};
    J_d=0;
    Jq{i}=0;
    Jqdot{i}=0;
    Jminqd2dot{i}=0;
    JmindFTtilde{i}=0;
    JpenhighFTtilde{i}=0;
    JpenhighFTtot{i}=0;
    JpenlMttilde{i}=0;
    
    
    t0=expdata.(nametrials{i}).q(1,1);
    tf=expdata.(nametrials{i}).q(end,1);
    
    tgrid{i}=[];
    tgrid{i}(1:(d+1):(N+1)*(d+1))=t0:((tf-t0)/N):tf;
    deltat=tgrid{i}(d+1+1)-tgrid{i}(1);
    tgrid{i}(2:(d+1):N*(d+1))=tgrid{i}(1:(d+1):N*(d+1))+tau_root(2)*deltat;
    tgrid{i}(3:(d+1):N*(d+1))=tgrid{i}(1:(d+1):N*(d+1))+tau_root(3)*deltat;
    tgrid{i}(4:(d+1):N*(d+1))=tgrid{i}(1:(d+1):N*(d+1))+tau_root(4)*deltat;

    h{i}=(tgrid{i}(end)-tgrid{i}(1))/N;
    
    if Options.optimizeMuscleProp
        if Options.startfromprevioussolution
            keyboard;
        end
        FT=MX.zeros(1,nmuscles-nmuscles_cut);
        FTtilde{i}=opti.variable(nmuscles-nmuscles_cut,N+1);
        opti.subject_to(bounds.FTtilde.lower <= FTtilde{i} <= ...
            bounds.FTtilde.upper);
        g1_names=[g1_names; repmat({'FTtilde bounds k'},(nmuscles-nmuscles_cut)*(N+1),1)];
        opti.set_initial(FTtilde{i},repmat(guess.FTtilde,1,N+1));
        FTtilde_k=MX.sym('FTtilde_k',nmuscles-nmuscles_cut,1);
        
        FTtilde_col{i}=opti.variable(nmuscles-nmuscles_cut,d*N);
        opti.subject_to(bounds.FTtilde.lower <= FTtilde_col{i} <= ...
            bounds.FTtilde.upper);
        g1_names=[g1_names; repmat({'FTtilde bounds j'},(nmuscles-nmuscles_cut)*d*N,1)];
        opti.set_initial(FTtilde_col{i},repmat(guess.FTtilde,1,d*N));
        FTtilde_j=MX.sym('FTtilde_j',nmuscles-nmuscles_cut,d);
    
        dFTtilde_col{i}=opti.variable(nmuscles-nmuscles_cut,d*N);
        opti.subject_to(bounds.dFTtilde.lower <= dFTtilde_col{i} <= ...
            bounds.dFTtilde.upper);
        g1_names=[g1_names; repmat({'dFTtilde bounds j'},(nmuscles-nmuscles_cut)*d*N,1)];
        opti.set_initial(dFTtilde_col{i},zeros(nmuscles-nmuscles_cut,d*N));
        dFTtilde_j=MX.sym('dFTtilde_j',nmuscles-nmuscles_cut,d);
    end

    QsQdot_prescribed{i}(:,1:2:7*2)=expdata.(nametrials{i}).q(:,2:8); %pelvis dofs and sacroiliac_flx are prescribed
    QsQdot_prescribed{i}(:,2:2:7*2)=expdata.(nametrials{i}).qdot(:,2:8);
    Qd2dot_prescribed{i}=expdata.(nametrials{i}).qd2dot(:,2:8);
    
    % if Options.startfromprevioussolution
    %     guess.QsQdots{i}(1:(d+1):(d+1)*N+1,:)=Options.prev_sol.sol_val.QsQdots_unsc{i}'./scaling.QsQdots;
    %     guess.QsQdots{i}(t_col_grid,:)=Options.prev_sol.sol_val.QsQdots_col_unsc{i}'./scaling.QsQdots;
    %     guess.Qd2dots{i}(t_col_grid,:)=Options.prev_sol.sol_val.Qd2dot_col_unsc{i}'/scaling.qd2dot;
    % else
        guess.QsQdots{i}(:,1:2:ndofs*2)=expdata.(nametrials{i}).q(:,9:15)/scaling.q;
        guess.QsQdots{i}(:,2:2:ndofs*2)=expdata.(nametrials{i}).qdot(:,9:15)/scaling.qdot;
        guess.Qd2dots{i}=expdata.(nametrials{i}).qd2dot(:,9:15)/scaling.qd2dot;
    % end
    
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
    Qd2dots_totrack_j=MX.sym('Qd2dots_totrack_j',ndofs,d);

    forces_prescribed_j=MX.sym('forces_prescribed_j',9,d);
    forces_prescribed{i}=expdata.(nametrials{i}).f(t_col_grid,2:end);
    
    if Options.optimizeMuscleProp
        akj=zeros(nmuscles-nmuscles_cut,d);
        lMTk=MX.sym('lMTk',nmuscles-nmuscles_cut,1);
        vMTk=MX.sym('vMTk',nmuscles-nmuscles_cut,1);
    end
    
    out_d{i}=MX.zeros(14,d);
    for j=1:d
        if Options.optimizeMuscleProp
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
        end
        all_QsQdot=MX.zeros(1,28);    
        all_QsQdot(1,1:14)=QsQdots_prescribed_j(1:14,j);
        all_QsQdot(1,15:28)=QsQdots_j(:,j)'.*scaling.QsQdots;
        
        all_Qd2dot=MX.zeros(1,14);
        all_Qd2dot(1:7)=Qd2dots_prescribed_j(1:7,j);
        all_Qd2dot(8:14)=Qd2dots_j(:,j)'.*scaling.qd2dot;
        
        if Options.optInertiaParam
            out=F([all_QsQdot';all_Qd2dot';forces_prescribed_j(:,j);inertiaParam_k.*scaling.inertiaparam']);
        else
            out=F([all_QsQdot';all_Qd2dot';forces_prescribed_j(:,j)]);
        end
        out_d{i}(:,j)=out;
        
        if Options.optimizeMuscleProp
            %dynamic constraints
            FTtildep_nsc = [FTtilde_k FTtilde_j]*scaling.FTtilde*C(:,j+1); %FTtilde_sol already in the original scale
            eq_constr{end+1}=(h{i}*dFTtilde_j(:,j)*scaling.dFTtilde-FTtildep_nsc)/scaling.dFTtilde;
            g2_names=[g2_names; repmat({'FTtilde cons'},nmuscles-nmuscles_cut,1)];
            
            %Path constraints
            % Hill equilibrium
            eq_constr{end+1}=Hilldiff_j{j};
            g2_names=[g2_names; repmat({'Hilldiff_eq'},nmuscles-nmuscles_cut,1)];
        end

        Qs_nsc = [QsQdots_k(1:2:end,:) QsQdots_j(1:2:end,:)].*scaling.QsQdots(1:2:end)'*C(:,j+1); %QsQdots_sol already in the original scale
        eq_constr{end+1}=(h{i}*QsQdots_j(2:2:end,j).*scaling.QsQdots(2:2:end)'-Qs_nsc)/scaling.QsQdots(2);
        g2_names=[g2_names; repmat({'Qs cons'},ndofs,1)];

        Qdots_nsc = [QsQdots_k(2:2:end,:) QsQdots_j(2:2:end,:)].*scaling.QsQdots(2:2:end)'*C(:,j+1); %QsQdots_sol already in the original scale
        eq_constr{end+1}=(h{i}*Qd2dots_j(:,j).*scaling.qd2dot'-Qdots_nsc)/scaling.qd2dot(1);
        g2_names=[g2_names; repmat({'Qsdot cons'},ndofs,1)];

        %Muscle force sharing
        %moment equilibrium
        if Options.individualpassiveprop==0;
            Itrial4passprop=list4passprop(i);
        else
            Itrial4passprop=i;
        end
        %hip flexion
        if Options.dofs_to_track(1)
            if Options.optimizeMuscleProp
                FT_hip_flx=FT_j{j}(mai(1).mus);
                T_hip_flx=FT_hip_flx'*MAj{j}.hip_flex;
            else
                T_hip_flx=0;
            end
            I=sum(Options.dofs_to_track(1:1));
            if Options.optimizePassiveJointEl
                PassiveM_hip_flx=-Kstiff_k_unsc(1,Itrial4passprop)*(QsQdots_j(1*2-1,j)*scaling.q-theta0_k_unsc(I,Itrial4passprop))-Ddamp_k_unsc(1,Itrial4passprop)*QsQdots_j(1*2,j)*scaling.qdot;
                if Options.optInertiapassiveParam
                    PassiveM_hip_flx=PassiveM_hip_flx-InertiapassiveParam_k_unsc(1,Itrial4passprop)*Qd2dots_j(1,j)*scaling.qd2dot;
                end
            else
                PassiveM_hip_flx=0;
            end
            eq_constr{end+1}=(out(8)-T_hip_flx-PassiveM_hip_flx)/scaling.T;
            g2_names=[g2_names; 'hip flex'];
        end

        %hip adduction
        if Options.dofs_to_track(2)
            if Options.optimizeMuscleProp
                FT_hip_add=FT_j{j}(mai(2).mus);
                T_hip_add=FT_hip_add'*MAj{j}.hip_add;
            else
                T_hip_add=0;
            end
            I=sum(Options.dofs_to_track(1:2));
            if Options.optimizePassiveJointEl
                PassiveM_hip_add=-Kstiff_k_unsc(1,Itrial4passprop)*(QsQdots_j(2*2-1,j)*scaling.q-theta0_k_unsc(I,Itrial4passprop))-Ddamp_k_unsc(1,Itrial4passprop)*QsQdots_j(2*2,j)*scaling.qdot;
                if Options.optInertiapassiveParam
                    PassiveM_hip_add=PassiveM_hip_add-InertiapassiveParam_k_unsc(1,Itrial4passprop)*Qd2dots_j(2,j)*scaling.qd2dot;
                end
            else
                PassiveM_hip_add=0;
            end
            eq_constr{end+1}=(out(9)-T_hip_add-PassiveM_hip_add)/scaling.T;
            g2_names=[g2_names; 'hip_add'];
        end

        %hip rotation
        if Options.dofs_to_track(3)
            if Options.optimizeMuscleProp
                FT_hip_int=FT_j{j}(mai(3).mus);
                T_hip_int=FT_hip_int'*MAj{j}.hip_int;
            else
                T_hip_int=0;
            end
            I=sum(Options.dofs_to_track(1:3));
            if Options.optimizePassiveJointEl
                PassiveM_hip_int=-Kstiff_k_unsc(1,Itrial4passprop)*(QsQdots_j(3*2-1,j)*scaling.q-theta0_k_unsc(I,Itrial4passprop))-Ddamp_k_unsc(1,Itrial4passprop)*QsQdots_j(3*2,j)*scaling.qdot;
                if Options.optInertiapassiveParam
                    PassiveM_hip_int=PassiveM_hip_int-InertiapassiveParam_k_unsc(1,Itrial4passprop)*Qd2dots_j(3,j)*scaling.qd2dot;
                end
            else
                PassiveM_hip_int=0;
            end
            eq_constr{end+1}=(out(10)-T_hip_int-PassiveM_hip_int)/scaling.T;
            g2_names=[g2_names; {'hip rot'}];
        end

        %knee flexion
        if Options.dofs_to_track(4)
            if Options.optimizeMuscleProp
                FT_knee_flx=FT_j{j}(mai(4).mus);
                T_knee_flx=FT_knee_flx'*MAj{j}.knee_flex;
            else
                T_knee_flx=0;
            end
            I=sum(Options.dofs_to_track(1:4));
            if Options.optimizePassiveJointEl
                PassiveM_knee_flx=-Kstiff_k_unsc(2,Itrial4passprop)*(QsQdots_j(4*2-1,j)*scaling.q-theta0_k_unsc(I,Itrial4passprop))-Ddamp_k_unsc(2,Itrial4passprop)*QsQdots_j(4*2,j)*scaling.qdot;
                if Options.optInertiapassiveParam
                    PassiveM_knee_flx=PassiveM_knee_flx-InertiapassiveParam_k_unsc(2,Itrial4passprop)*Qd2dots_j(4,j)*scaling.qd2dot;
                end
            else
                PassiveM_knee_flx=0;
            end
            eq_constr{end+1}=(out(11)-T_knee_flx-PassiveM_knee_flx)/scaling.T;
            g2_names=[g2_names; {'knee flex'}];
        end

        %ankle flexion
        if Options.dofs_to_track(5)
            if Options.optimizeMuscleProp
                FT_ankle_flx=FT_j{j}(mai(5).mus);
                T_ankle_flx=FT_ankle_flx'*MAj{j}.ankle_flex;
            else
                T_ankle_flx=0;
            end
            I=sum(Options.dofs_to_track(1:5));
            if Options.optimizePassiveJointEl
                PassiveM_ankle_flx=-Kstiff_k_unsc(3,Itrial4passprop)*(QsQdots_j(5*2-1,j)*scaling.q-theta0_k_unsc(I,Itrial4passprop))-Ddamp_k_unsc(3,Itrial4passprop)*QsQdots_j(5*2,j)*scaling.qdot;
                if Options.optInertiapassiveParam
                    PassiveM_ankle_flx=PassiveM_ankle_flx-InertiapassiveParam_k_unsc(3,Itrial4passprop)*Qd2dots_j(5,j)*scaling.qd2dot;
                end
            else
                PassiveM_ankle_flx=0;
            end
            eq_constr{end+1}=(out(12)-T_ankle_flx-PassiveM_ankle_flx)/scaling.T;
            g2_names=[g2_names; {'ankle flex'}];
        end

        %ankle adduction
        if Options.dofs_to_track(6)
            if Options.optimizeMuscleProp
                FT_ankle_add=FT_j{j}(mai(6).mus);
                T_ankle_add=FT_ankle_add'*MAj{j}.ankle_add;
            else
                T_ankle_add=0;
            end
            I=sum(Options.dofs_to_track(1:6));
            if Options.optimizePassiveJointEl
                PassiveM_ankle_add=-Kstiff_k_unsc(3,Itrial4passprop)*(QsQdots_j(6*2-1,j)*scaling.q-theta0_k_unsc(I,Itrial4passprop))-Ddamp_k_unsc(3,Itrial4passprop)*QsQdots_j(6*2,j)*scaling.qdot;
                if Options.optInertiapassiveParam
                    PassiveM_ankle_add=PassiveM_ankle_add-InertiapassiveParam_k_unsc(3,Itrial4passprop)*Qd2dots_j(6,j)*scaling.qd2dot;
                end
            else
                PassiveM_ankle_add  =0;
            end
            eq_constr{end+1}=(out(13)-T_ankle_add-PassiveM_ankle_add)/scaling.T;
            g2_names=[g2_names; {'ankle add'}];
        end

        %ankle int
        if Options.dofs_to_track(7)
            if Options.optimizeMuscleProp
                FT_ankle_int=FT_j{j}(mai(7).mus);
                T_ankle_int=FT_ankle_int'*MAj{j}.ankle_int;
            else
                T_ankle_int=0;
            end
            I=sum(Options.dofs_to_track(1:7));
            if Options.optimizePassiveJointEl
                PassiveM_ankle_int=-Kstiff_k_unsc(3,Itrial4passprop)*(QsQdots_j(7*2-1,j)*scaling.q-theta0_k_unsc(I,Itrial4passprop))-Ddamp_k_unsc(3,Itrial4passprop)*QsQdots_j(7*2,j)*scaling.qdot;
                if Options.optInertiapassiveParam
                    PassiveM_ankle_int=PassiveM_ankle_int-InertiapassiveParam_k_unsc(3,Itrial4passprop)*Qd2dots_j(7,j)*scaling.qd2dot;
                end
            else
                PassiveM_ankle_int=0;
            end
            eq_constr{end+1}=(out(14)-T_ankle_int-PassiveM_ankle_int)/scaling.T;
            g2_names=[g2_names; {'ankle int'}];
        end

        %Cost function
        J_d=J_d+B(j+1)*W.qtrack*sum((QsQdots_j(1:2:end,j)-QsQdots_totrack_j(1:2:end,j)).^2)*h{i};
        J_d=J_d+B(j+1)*W.qdottrack*sum((QsQdots_j(2:2:end,j)-QsQdots_totrack_j(2:2:end,j)).^2)*h{i};

        if Options.trackqd2dot
            J_d=J_d+B(j+1)*W.mindstate*sum((Qd2dots_j(:,j)-Qd2dots_totrack_j(:,j)).^2)*h{i};
        else
            J_d=J_d+B(j+1)*W.mindstate*sum(Qd2dots_j(:,j).^2)*h{i};         
        end
        if Options.optimizeMuscleProp
            J_d=J_d+B(j+1)*W.mindstate*sum(dFTtilde_j(:,j).^2)*h{i};
            
            if Options.penalizeHighFTtilde
                J_d=J_d+B(j+1)*W.penalizeHighFTtilde*sum(FTtilde_j(:,j).^2)*h{i};
            end
            if Options.penalizeFTtot
                J_d=J_d+B(j+1)*W.penalizeFTtot*sum(FT_j{j}.^2)*h{i};
            end
            if Options.penalizeoutoflMtilde1
                J_d=J_d+B(j+1)*W.penalizeoutoflMtilde1*sum((lMtilde_j{j}-1).^2)*h{i};
            end
        end
        %for debug
        Jq{i}=Jq{i}+B(j+1)*W.qtrack*sum((QsQdots_j(1:2:end,j)-QsQdots_totrack_j(1:2:end,j)).^2)*h{i}; %track q
        Jqdot{i}=Jqdot{i}+B(j+1)*W.qdottrack*sum((QsQdots_j(2:2:end,j)-QsQdots_totrack_j(2:2:end,j)).^2)*h{i}; %track qdot
        if Options.trackqd2dot
            Jminqd2dot{i}=Jminqd2dot{i}+B(j+1)*W.mindstate*sum((Qd2dots_j(:,j)-Qd2dots_totrack_j(:,j)).^2)*h{i};
        else
            Jminqd2dot{i}=Jminqd2dot{i}+B(j+1)*W.mindstate*sum(Qd2dots_j(:,j).^2)*h{i};
        end

        if Options.optimizeMuscleProp
            JmindFTtilde{i}=JmindFTtilde{i}+B(j+1)*W.mindstate*sum(dFTtilde_j(:,j).^2)*h{i};
            if Options.penalizeHighFTtilde
                JpenhighFTtilde{i}=JpenhighFTtilde{i}+B(j+1)*W.penalizeHighFTtilde*sum(FTtilde_j(:,j).^2)*h{i};
            end
            if Options.penalizeFTtot
                JpenhighFTtot{i}=JpenhighFTtot{i}+B(j+1)*W.penalizeFTtot*sum(FT_j{j}.^2)*h{i};
            end
            if Options.penalizeoutoflMtilde1
                JpenlMttilde{i}=JpenlMttilde{i}+B(j+1)*W.penalizeoutoflMtilde1*sum((lMtilde_j{j}-1).^2)*h{i};
            end
        end
    end
    eq_constr = vertcat(eq_constr{:});
    ineq_constr=vertcat(ineq_constr{:});
    
    if Options.optInertiaParam
        inputs_f_coll={inertiaParam_k};
        inputs_f_coll_map={repmat(inertiaParam,1,N)};
        inputs_f_coll2={inertiaParam_k};
        inputs_f_coll_map2={repmat(inertiaParam,1,N)};
    else
        inputs_f_coll=[];
        inputs_f_coll_map=[];
        inputs_f_coll2=[];
        inputs_f_coll_map2=[];
    end
    outputs_f_coll=[];

    inputs_f_coll=[inputs_f_coll {QsQdots_k,QsQdots_j,QsQdots_totrack_j,...
        QsQdots_prescribed_j, Qd2dots_j,Qd2dots_prescribed_j, ...
        forces_prescribed_j}];
    outputs_f_coll=[outputs_f_coll {eq_constr,J_d}];
    inputs_f_coll_map=[inputs_f_coll_map ...
        {QsQdots{i}(:,1:end-1),QsQdots_col{i},...
        guess.QsQdots{i}(t_col_grid,:)',... %guess.QsQdots is already scaled (this is for QsQdots_totrack_j)
        QsQdot_prescribed{i}(t_col_grid,:)',... %QsQdot_prescribed is unscaled
        Qd2dots_col{i},...
        Qd2dot_prescribed{i}(t_col_grid,:)',forces_prescribed{i}'}];

    outputs_f_coll_map = cell(1, 2); % Preallocate a cell array for two outputs
    outputs_f_coll2={Jq{i},Jqdot{i},Jminqd2dot{i},JmindFTtilde{i},JpenhighFTtilde{i},JpenhighFTtot{i},JpenlMttilde{i},out_d{i}};
    outputs_f_coll_map2=cell(1,8);

    if Options.trackqd2dot
        inputs_f_coll=[inputs_f_coll {Qd2dots_totrack_j}];
        outputs_f_coll=[outputs_f_coll];
        inputs_f_coll_map=[inputs_f_coll_map {guess.Qd2dots{i}(t_col_grid,:)'}];

        outputs_f_coll_map=outputs_f_coll_map;
    end


    if Options.optimizeMuscleProp
        if Options.useRigidTendon
            keyboard;
        else
            inputs_f_coll=[inputs_f_coll {FTtilde_k,FTtilde_j,...
                dFTtilde_j, lTs_k, lM0_k}];
            outputs_f_coll=[outputs_f_coll {[FT_j{1} FT_j{2} FT_j{3}],...
                [lMtilde_j{1} lMtilde_j{2} lMtilde_j{3}],...
                [lTtilde_j{1} lTtilde_j{2} lTtilde_j{3}]}];
            inputs_f_coll_map=[inputs_f_coll_map {FTtilde{i}(:,1:end-1),...
                FTtilde_col{i}, dFTtilde_col{i}, ...
                repmat(lTs,1,N),repmat(lM0,1,N)}];

            inputs_f_coll=[inputs_f_coll {fiber_damping_k}];
            if Options.optimizefiberdamping
                inputs_f_coll_map=[inputs_f_coll_map {repmat(fiber_damping,1,N)}];
            else
                inputs_f_coll_map=[inputs_f_coll_map {zeros(1,N)}];
            end
            inputs_f_coll=[inputs_f_coll {tendon_damping_k}];
            if Options.optimizetendondamping
                inputs_f_coll_map=[inputs_f_coll_map {repmat(tendon_damping,1,N)}];  
            else
                inputs_f_coll_map=[inputs_f_coll_map {zeros(1,N)}];  
            end
            outputs_f_coll_map = [outputs_f_coll_map cell(1, 3)]; % Preallocate a cell array for five outputs
            
            %for debugging
            inputs_f_coll2=inputs_f_coll;
            inputs_f_coll_map2=inputs_f_coll_map;
        end
    end

    if Options.optimizePassiveJointEl
        inputs_f_coll=[inputs_f_coll {Kstiff_k,Ddamp_k,theta0_k}];
        inputs_f_coll_map=[inputs_f_coll_map {repmat(Kstiff,1,N),repmat(Ddamp,1,N),repmat(theta0,1,N)}];
        if Options.optInertiapassiveParam
            inputs_f_coll=[inputs_f_coll {InertiapassiveParam_k}];
            inputs_f_coll_map=[inputs_f_coll_map {repmat(InertiapassiveParam,1,N)}];
        end
        
        %for debugging
        inputs_f_coll2=inputs_f_coll;
        inputs_f_coll_map2=inputs_f_coll_map;
    end

    f_coll = Function('f_coll',inputs_f_coll,outputs_f_coll);
    f_coll_map = f_coll.map(N,ParallelMode,NThreads);
    [outputs_f_coll_map{:}] = f_coll_map(inputs_f_coll_map{:});


    f_coll2 = Function('f_coll2',inputs_f_coll2,outputs_f_coll2);
    f_coll_map2=f_coll2.map(N,ParallelMode,NThreads);
    [outputs_f_coll_map2{:}] = f_coll_map2(inputs_f_coll_map2{:});
    

    [Jqall{i},Jqdotall{i},Jminqd2dotall{i},JmindFTtildeall{i},JpenhighFTtildeall{i},JpenhighFTtotall{i},JpenlMttildeall{i},out_dall{i}] = outputs_f_coll_map2{:};
    if Options.optimizeMuscleProp
        if Options.useRigidTendon
            keyboard;
        else
        end
        [coll_eq_constr{i}, Jall{i}, FT_all{i}, lMtilde_all{i}, lTtilde_all{i}] = outputs_f_coll_map{:};
    else
        [coll_eq_constr{i}, Jall{i}] = outputs_f_coll_map{:};
    end

    opti.subject_to(coll_eq_constr{i}==0);
    % Add continuity constraints (next interval starts with end values of 
    % states from previous interval)
    if Options.optimizeMuscleProp
        if Options.useRigidTendon
            %no need for continuity constraints when using rigid
            %tendon
        else
            for k=1:N
                % Variables within current mesh interval
                % States    
                FTtilde_kj = [FTtilde{i}(:,k), FTtilde_col{i}(:,(k-1)*d+1:k*d)];
                opti.subject_to(FTtilde{i}(:,k+1) == FTtilde_kj*D); % scaled
                g3_names=[g3_names; repmat({'continuity constr FTtilde'},nmuscles-nmuscles_cut,1)];
            end
        end  
    end
    for k=1:N
        QsQdots_kj=[QsQdots{i}(:,k), QsQdots_col{i}(:,(k-1)*d+1:k*d)];
        opti.subject_to((QsQdots{i}(:,k+1))/10==(QsQdots_kj/10)*D); %scaled
        g3_names=[g3_names; repmat({'continuity constr QsQdots'},ndofs*2,1)];
    end
    J_i=J_i+sum(Jall{i});
    
end

if Options.optimizePassiveJointEl&(Options.pendevPassive==1)

    Kstiff_sum=sum(Kstiff,2);
    Ddamp_sum=sum(Ddamp,2);

    Kstiffmean=Kstiff_sum/size(Kstiff,2);
    Ddampmean=Ddamp_sum/size(Ddamp,2);
    if Options.optInertiapassiveParam
        InertiapassiveParam_sum=sum(InertiapassiveParam,2);
        Inertiapassparammean=InertiapassiveParam_sum/size(InertiapassiveParam,2);
    end
    J_passjointpen=0;
    for i=1:size(Kstiff,2)
        J=J+W.Kstiffpen*sum((Kstiff(:,i)-Kstiffmean).^2);
        J=J+W.Ddamppen*sum((Ddamp(:,i)-Ddampmean).^2);
        J=J+W.theta0pen*sum((theta0(:,i)-guess.theta0).^2);
        J_passjointpen=J_passjointpen+W.Kstiffpen*sum((Kstiff(:,i)-Kstiffmean).^2)+W.Ddamppen*sum((Ddamp(:,i)-Ddampmean).^2)+W.theta0pen*sum((theta0(:,i)-guess.theta0).^2);
        if Options.optInertiapassiveParam
            J=J+W.Inertiapassparampen*sum((InertiapassiveParam(:,i)-Inertiapassparammean).^2);
            J_passjointpen=J_passjointpen+W.Inertiapassparampen*sum((InertiapassiveParam(:,i)-Inertiapassparammean).^2);
        end
    end
end


g_names=[g1_names; repmat([g2_names_ineq; g2_names],N,1); g3_names];  

J=J+J_i;

%% Solve OCP
opti.minimize(J);
options.ipopt.hessian_approximation = 'limited-memory'; %'exact'; %
options.ipopt.mu_strategy      = 'adaptive';
options.ipopt.max_iter = 2000;
options.ipopt.tol = 1e-5;   
opti.solver('ipopt', options);  
sol=opti.solve();

if opti.stats.success==1
    sol_val.opt_x=opti.value(opti.x);
    sol_val.lam_g=opti.value(opti.lam_g);
else
    sol_val.opt_x=opti.debug.value(opti.x);
    sol_val.lam_g=opti.debug.value(opti.lam_g);
end


if Options.optInertiaParam
    if opti.stats.success==1
        sol_val.inertiaParam=sol.value(inertiaParam);
    else
        sol_val.inertiaParam=opti.debug.value(inertiaParam);
    end
    sol_val.inertiaParam_unsc=sol_val.inertiaParam.*scaling.inertiaparam';
end

if Options.optimizeMuscleProp
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
end
sol_val.Options=Options;
sol_val.scaling=scaling;
sol_val.ndofs=ndofs;

sol_val.tgrid=tgrid;
% sol_val.t2plot=t_col_grid;
sol_val.W=W;
if Options.optimizeMuscleProp
    sol_val.muscle_names=muscle_name;
end
sol_val.stats=opti.stats;


if Options.optimizePassiveJointEl
    if opti.stats.success
        sol_val.Kstiff=sol.value(Kstiff);
        sol_val.Ddamp=sol.value(Ddamp);
        sol_val.theta0=sol.value(theta0);
        if Options.optInertiapassiveParam
            sol_val.InertiapassiveParam = sol.value(InertiapassiveParam);
        end
    else 
        sol_val.Kstiff=opti.debug.value(Kstiff);
        sol_val.Ddamp=opti.debug.value(Ddamp);
        sol_val.theta0=opti.debug.value(theta0);
        if Options.optInertiapassiveParam
            sol_val.InertiapassiveParam = opti.debug.value(InertiapassiveParam);
        end
    end
    sol_val.Kstiff_unsc=sol_val.Kstiff*scaling.Kstiff;
    sol_val.Ddamp_unsc=sol_val.Ddamp*scaling.Ddamp;
    sol_val.theta0_unsc=sol_val.theta0*scaling.theta;
    if Options.optInertiapassiveParam
        sol_val.InertiapassiveParam_unsc=sol_val.InertiapassiveParam*scaling.Inertiapassparam;
    end
end


if opti.stats.success
    
    for i=1:size(QsQdots,2)
        sol_val.tgrid{i}=tgrid{i};
        sol_val.t2plot{i}=t_col_grid;
        sol_val.QsQdots{i}=sol.value(QsQdots{i});
        sol_val.QsQdots_col{i}=sol.value(QsQdots_col{i});
        sol_val.Qd2dot_col{i}=sol.value(Qd2dots_col{i});
        if ~isempty(sol_val.QsQdots{i})
            sol_val.QsQdots_unsc{i}=sol_val.QsQdots{i}.*scaling.QsQdots';
            sol_val.QsQdots_col_unsc{i}=sol_val.QsQdots_col{i}.*scaling.QsQdots';
            sol_val.Qd2dot_col_unsc{i}=sol_val.Qd2dot_col{i}*scaling.qd2dot;
        end
        sol_val.tgrid_col{i}=tgrid{i}(t_col_grid);
    end
    
    if Options.optimizeMuscleProp
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
    end
else  
    
    for i=1:size(QsQdots,2)
        sol_val.tgrid{i}=tgrid{i};
        sol_val.t2plot{i}=t_col_grid;

        sol_val.QsQdots{i}=opti.debug.value(QsQdots{i});
        sol_val.QsQdots_col{i}=opti.debug.value(QsQdots_col{i});
        sol_val.Qd2dot_col{i}=opti.debug.value(Qd2dots_col{i});
        if ~isempty(sol_val.QsQdots{i})
            sol_val.QsQdots_unsc{i}=sol_val.QsQdots{i}.*scaling.QsQdots';
            sol_val.QsQdots_col_unsc{i}=sol_val.QsQdots_col{i}.*scaling.QsQdots';
            sol_val.Qd2dot_col_unsc{i}=sol_val.Qd2dot_col{i}*scaling.qd2dot;
        end
        sol_val.tgrid_col{i}=tgrid{i}(t_col_grid);
    end
    if Options.optimizeMuscleProp
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

end
sol_val.guess=guess;
sol_val.bounds=bounds;
sol_val.QsQdot_prescribed=QsQdot_prescribed;
sol_val.Qd2dot_prescribed=Qd2dot_prescribed;
sol_val.force_prescribed=forces_prescribed;
sol_val.name_dofs=name_dofs;
sol_val.nametrials=nametrials;
if Options.optimizeMuscleProp
    sol_val.MTparam=MTparam;
end

if Options.optimizeMuscleProp
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
                                           zeros(nmuscles-nmuscles_cut,1),sol_val.FTtilde_unsc{i}(:,(k-1)*(d+1)+j+1),...
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
end

J_opt=0;
J_opt_trackq=0;
J_opt_trackqdot=0;
J_opt_minqs2dot=0;
J_optpassjoint=0;
J_optpassjointpen=0;
Kstiff_sum_opt=0;
Ddamp_sum_opt=0;
%plot cost function terms
for i=1:length(Jq)
    figure
    plot(opti.value(Jqall{i}));
    hold all
    plot(opti.value(Jqdotall{i}));
    plot(opti.value(Jminqd2dotall{i}));
    plot(opti.value(JmindFTtildeall{i}));
    plot(opti.value(JpenhighFTtildeall{i}));
    plot(opti.value(JpenhighFTtotall{i}));
    plot(opti.value(JpenlMttildeall{i}));
    legend({'Jq','Jqdot','Jminqd2dot','JmindFTtilde','JpenhighFTtilde','JpenhighFTtot','JpenlMttilde'});
    J_opt=J_opt+sum(opti.value(Jqall{i}))+sum(opti.value(Jqdotall{i}))+...
        sum(opti.value(Jminqd2dotall{i}))+sum(opti.value(JmindFTtildeall{i}))+...
        sum(opti.value(JpenhighFTtildeall{i}))+sum(opti.value(JpenhighFTtotall{i}))+...
        sum(opti.value(JpenlMttildeall{i}));

    J_opt_trackq=J_opt_trackq+sum(opti.value(Jqall{i}));
    J_opt_trackqdot=J_opt_trackqdot+sum(opti.value(Jqdotall{i}));
    J_opt_minqs2dot=J_opt_minqs2dot+sum(opti.value(Jminqd2dotall{i}));
end
if Options.optimizePassiveJointEl
    if Options.minPassiveProp
        opti.value(J_passjoint{i});
        J_opt=J_opt+opti.value(J_passjoint);
        J_optpassjoint=J_optpassjoint+opti.value(J_passjoint);
    end
    if Options.pendevPassive
        J_opt=J_opt+opti.value(J_passjointpen);
        J_optpassjointpen=J_optpassjointpen+opti.value(J_passjointpen);
        %for debug
        Kstiff_sum_opt=Kstiff_sum_opt+sum(sol_val.Kstiff,2);
        Ddamp_sum_opt=Ddamp_sum_opt+sum(sol_val.Ddamp,2);
        %for debug
        Kstiffmean_opt=Kstiff_sum_opt/size(Kstiff,2);
        Ddampmean_opt=Ddamp_sum_opt/size(Ddamp,2);
        
        J_passjointpen_test=0;
        for i=1:size(Kstiff,2)
            J_passjointpen_test=J_passjointpen_test+W.Kstiffpen*sum((sol_val.Kstiff-Kstiffmean_opt).^2)+W.Ddamppen*sum((sol_val.Ddamp(:,i)-Ddampmean_opt).^2)+W.theta0pen*sum((sol_val.theta0(:,i)-guess.theta0).^2);
        end
        if Options.optInertiapassiveParam
            InertiapassiveParammean_opt=mean(sol_val.InertiapassiveParam,2);
            J_passjointpen_test=J_passjointpen_test+W.Inertiapassparampen*sum((sol_val.InertiapassiveParam-InertiapassiveParammean_opt).^2);
        end
    end
end
if Options.optimizeMuscleProp&Options.optimizelM0
    opti.value(J_lM0);
    J_opt=J_opt+opti.value(J_lM0);
end
if Options.optInertiaParam
    J_opt=J_opt+opti.value(J_inertiaParam);
end
figure;
bar([J_opt_trackq J_opt_trackqdot J_opt_minqs2dot opti.value(J_inertiaParam) J_optpassjoint J_optpassjointpen]);
set(gca,'XTickLabels',{'qtrak','qdottrack','minqd2dot','penInertia','minPass','penPass'})

%% write output kinematics
q.data(:,1)=sol_val.tgrid{1}(sol_val.t2plot{1});
q.data(:,2:8)=sol_val.QsQdot_prescribed{1}(sol_val.t2plot{1},1:2:end);
q.data(:,9:15)=sol_val.QsQdots_col_unsc{1}(1:2:end,:)';
q.data(:,[2:4 8:end])=q.data(:,[2:4 8:end])*180/pi;
q.labels={'time', 'sacrum_pitch','sacrum_roll','sacrum_yaw','sacrum_x','sacrum_y','sacrum_z','sacroiliac_flx', ...
    'hip_flx','hip_add','hip_int','knee_flx','ankle_flx','ankle_add','ankle_int'};
write_motionFile(q,'testoutmot.mot');

%% recompute constraint values
for i=1:length(sol_val.QsQdot_prescribed)
    if ~isempty(sol_val.QsQdot_prescribed{i})
        all_QsQdot_opt{i}(:,1:14)=sol_val.QsQdot_prescribed{i}(sol_val.t2plot{i},:);
        all_QsQdot_opt{i}(:,15:28)=sol_val.QsQdots_col_unsc{i}';

        all_Qd2dot_opt{i}(:,1:7)=sol_val.Qd2dot_prescribed{i}(sol_val.t2plot{i},:);
        all_Qd2dot_opt{i}(:,8:14)=sol_val.Qd2dot_col_unsc{i}';

        for k=1:size(all_Qd2dot_opt{i},1)
            if Options.optInertiaParam
                out_opt{i}(k,:)=full(F([all_QsQdot_opt{i}(k,:)';all_Qd2dot_opt{i}(k,:)';forces_prescribed{i}(k,:)';sol_val.inertiaParam_unsc]));
            else
                out_opt{i}(k,:)=full(F([all_QsQdot_opt{i}(k,:)';all_Qd2dot_opt{i}(k,:)';forces_prescribed{i}(k,:)']));
            end
        end
    end
end
sol_val.out_opt=out_opt;

%Muscle force sharing
%moment equilibrium
if Options.optimizePassiveJointEl
    if Options.individualpassiveprop==0;
        Itrial4passprop=list4passprop(i);
    else
        Itrial4passprop=i;
    end
end
for i=1:length(sol_val.QsQdots_col_unsc)

    %hip flexion
    if Options.dofs_to_track(1)
        if Options.optimizeMuscleProp
            FT_hip_flx_opt=sol_val.FT_all{i}(mai(1).mus,:); 
            T_hip_flx_opt=FT_hip_flx_opt'.*MA_opt.hip_flex;
        else
            T_hip_flx_opt=zeros(N*d,1);
        end
        I=sum(Options.dofs_to_track(1:1));
        if Options.optimizePassiveJointEl
            PassiveM_hip_flx_opt=-sol_val.Kstiff_unsc(1,Itrial4passprop)*(sol_val.QsQdots_col_unsc{i}(1,:)'-sol_val.theta0_unsc(I,Itrial4passprop))-sol_val.Ddamp_unsc(1,Itrial4passprop)*sol_val.QsQdots_col_unsc{i}(2,:)';
            if Options.optInertiapassiveParam
                PassiveM_hip_flx_opt=PassiveM_hip_flx_opt-sol_val.InertiapassiveParam_unsc(1,Itrial4passprop)*(sol_val.Qd2dot_col_unsc{i}(:,1));
            end
        else
            PassiveM_hip_flx_opt=zeros(d*N,1);   
        end
        sol_val.PassiveM_hip_flx_opt{i}=PassiveM_hip_flx_opt;
        eq_constr_opt{i}(:,1)=(out_opt{i}(:,8)-sum(T_hip_flx_opt,2)-PassiveM_hip_flx_opt)/scaling.T;

    end
    
    %hip adduction
    if Options.dofs_to_track(2)
        if Options.optimizeMuscleProp
            FT_hip_add_opt=sol_val.FT_all{i}(mai(2).mus,:);
            T_hip_add_opt=FT_hip_add_opt'.*MA_opt.hip_add;
        else
            T_hip_add_opt=zeros(N*d,1);
        end
        I=sum(Options.dofs_to_track(1:2));
        if Options.optimizePassiveJointEl
            PassiveM_hip_add_opt=-sol_val.Kstiff_unsc(1,Itrial4passprop)*(sol_val.QsQdots_col_unsc{i}(3,:)'-sol_val.theta0_unsc(I,Itrial4passprop))-sol_val.Ddamp_unsc(1,Itrial4passprop)*sol_val.QsQdots_col_unsc{i}(4,:)';
            if Options.optInertiapassiveParam
                PassiveM_hip_add_opt=PassiveM_hip_add_opt-sol_val.InertiapassiveParam_unsc(1,Itrial4passprop)*sol_val.Qd2dot_col_unsc{i}(:,2);
            end
        else
            PassiveM_hip_add_opt=zeros(d*N,1);
        end
        sol_val.PassiveM_hip_add_opt{i}=PassiveM_hip_add_opt;
        eq_constr_opt{i}(:,2)=(out_opt{i}(:,9)-sum(T_hip_add_opt,2)-PassiveM_hip_add_opt)/scaling.T;
    end

    %hip rotation
    if Options.dofs_to_track(3)
        if Options.optimizeMuscleProp
            FT_hip_int_opt=sol_val.FT_all{i}(mai(3).mus,:);
            T_hip_int_opt=FT_hip_int_opt'.*MA_opt.hip_int;
        else
            T_hip_int_opt=zeros(d*N,1);
        end
        I=sum(Options.dofs_to_track(1:3));
        if Options.optimizePassiveJointEl
            PassiveM_hip_int_opt=-sol_val.Kstiff_unsc(1,Itrial4passprop)*(sol_val.QsQdots_col_unsc{i}(5,:)'-sol_val.theta0_unsc(I,Itrial4passprop))-sol_val.Ddamp_unsc(1,Itrial4passprop)*sol_val.QsQdots_col_unsc{i}(6,:)';
            if Options.optInertiapassiveParam
                PassiveM_hip_int_opt=PassiveM_hip_int_opt-sol_val.InertiapassiveParam_unsc(1,Itrial4passprop)*sol_val.Qd2dot_col_unsc{i}(:,3);
            end
        else
            PassiveM_hip_int_opt=zeros(d*N,1);
        end
        sol_val.PassiveM_hip_int_opt{i}=PassiveM_hip_int_opt;
        eq_constr_opt{i}(:,3)=(out_opt{i}(:,10)-sum(T_hip_int_opt,2)-PassiveM_hip_int_opt)/scaling.T;

    end

    %knee flexion
    if Options.dofs_to_track(4)
        if Options.optimizeMuscleProp
            FT_knee_flx_opt=sol_val.FT_all{i}(mai(4).mus,:);
            T_knee_flx_opt=FT_knee_flx_opt'.*MA_opt.knee_flex;
        else
            T_knee_flx_opt=zeros(d*N,1);
        end
        I=sum(Options.dofs_to_track(1:4));
        if Options.optimizePassiveJointEl
            PassiveM_knee_flx_opt=-sol_val.Kstiff_unsc(2,Itrial4passprop)*(sol_val.QsQdots_col_unsc{i}(7,:)'-sol_val.theta0_unsc(I,Itrial4passprop))-sol_val.Ddamp_unsc(2,Itrial4passprop)*sol_val.QsQdots_col_unsc{i}(8,:)';
            if Options.optInertiapassiveParam
                PassiveM_knee_flx_opt=PassiveM_knee_flx_opt-sol_val.InertiapassiveParam_unsc(2,Itrial4passprop)*sol_val.Qd2dot_col_unsc{i}(:,4);
            end
        else
            PassiveM_knee_flx_opt=zeros(d*N,1);
        end
        sol_val.PassiveM_knee_flx_opt{i}=PassiveM_knee_flx_opt;
        eq_constr_opt{i}(:,4)=(out_opt{i}(:,11)-sum(T_knee_flx_opt,2)-PassiveM_knee_flx_opt)/scaling.T;

    end

    %ankle flexion
    if Options.dofs_to_track(5)
        if Options.optimizeMuscleProp
            FT_ankle_flx_opt=sol_val.FT_all{i}(mai(5).mus,:);
            T_ankle_flx_opt=FT_ankle_flx_opt'.*MA_opt.ankle_flex;
        else
            T_ankle_flx_opt=zeros(d*N,1);
        end
        I=sum(Options.dofs_to_track(1:5));
        if Options.optimizePassiveJointEl
            PassiveM_ankle_flx_opt=-sol_val.Kstiff_unsc(3,Itrial4passprop)*(sol_val.QsQdots_col_unsc{i}(9,:)'-sol_val.theta0_unsc(I,Itrial4passprop))-sol_val.Ddamp_unsc(3,Itrial4passprop)*sol_val.QsQdots_col_unsc{i}(10,:)';
            if Options.optInertiapassiveParam
                PassiveM_ankle_flx_opt=PassiveM_ankle_flx_opt-sol_val.InertiapassiveParam_unsc(3,Itrial4passprop)*sol_val.Qd2dot_col_unsc{i}(:,5);
            end
        else
            PassiveM_ankle_flx_opt=zeros(d*N,1);       
        end
        sol_val.PassiveM_ankle_flx_opt{i}=PassiveM_ankle_flx_opt;
        eq_constr_opt{i}(:,5)=(out_opt{i}(:,12)-sum(T_ankle_flx_opt,2)-PassiveM_ankle_flx_opt)/scaling.T;
    end

    %ankle adduction
    if Options.dofs_to_track(6)
        if Options.optimizeMuscleProp
            FT_ankle_add_opt=sol_val.FT_all{i}(mai(6).mus,:);
            T_ankle_add_opt=FT_ankle_add_opt'.*MA_opt.ankle_add;
        else
            T_ankle_add_opt=zeros(d*N,1);
        end
        I=sum(Options.dofs_to_track(1:6));
        if Options.optimizePassiveJointEl
            PassiveM_ankle_add_opt=-sol_val.Kstiff_unsc(3,Itrial4passprop)*(sol_val.QsQdots_col_unsc{i}(11,:)'-sol_val.theta0_unsc(I,Itrial4passprop))-sol_val.Ddamp_unsc(3,Itrial4passprop)*sol_val.QsQdots_col_unsc{i}(12,:)';
            if Options.optInertiapassiveParam
                PassiveM_ankle_add_opt=PassiveM_ankle_add_opt-sol_val.InertiapassiveParam_unsc(3,Itrial4passprop)*sol_val.Qd2dot_col_unsc{i}(:,6);
            end
        else
            PassiveM_ankle_add_opt=zeros(d*N,1);  
        end
        sol_val.PassiveM_ankle_add_opt{i}=PassiveM_ankle_add_opt;
        eq_constr_opt{i}(:,6)=(out_opt{i}(:,13)-sum(T_ankle_add_opt,2)-PassiveM_ankle_add_opt)/scaling.T;
    end

    %ankle int
    if Options.dofs_to_track(7)
        if Options.optimizeMuscleProp   
            FT_ankle_int_opt=sol_val.FT_all{i}(mai(7).mus,:);
            T_ankle_int_opt=FT_ankle_int_opt'.*MA_opt.ankle_int;
        else
             T_ankle_int_opt=zeros(d*N,1);
        end
        I=sum(Options.dofs_to_track(1:7));
        if Options.optimizePassiveJointEl
            PassiveM_ankle_int_opt=-sol_val.Kstiff_unsc(3,Itrial4passprop)*(sol_val.QsQdots_col_unsc{i}(13,:)'-sol_val.theta0_unsc(I,Itrial4passprop))-sol_val.Ddamp_unsc(3,Itrial4passprop)*sol_val.QsQdots_col_unsc(14,:)';
            if Options.optInertiapassiveParam
                PassiveM_ankle_int_opt=PassiveM_ankle_int_opt-sol_val.InertiapassiveParam_unsc(3,Itrial4passprop)*sol_val.Qd2dot_col_unsc{i}(:,7);
            end
        else
            PassiveM_ankle_int_opt=zeros(d*N,1);
        end
        sol_val.PassiveM_ankle_int_opt{i}=PassiveM_ankle_int_opt;
        eq_constr_opt{i}(:,7)=(out_opt{i}(:,14)-sum(T_ankle_int_opt,2)-PassiveM_ankle_int_opt)/scaling.T;
    end
end

%dynamic and continuity constraints
for i=1:length(sol_val.QsQdots_unsc)
    if ~isempty(sol_val.QsQdots_unsc{i})
        for k=1:N
            for j=1:d
                if Options.optimizeMuscleProp
                    FTtildep_opt_nsc{i}(k,:) = [sol_val.FTtilde_unsc{i}(:,(k-1)*(d+1)+1:k*(d+1))]*C(:,j+1); %FTtilde_sol already in the original scale
                    eq_constr_dyn_FTtilde{i}(k,:)=(h{i}*sol_val.dFTtilde_unsc{i}(:,(k-1)*d+j)-FTtildep_opt_nsc{i}(k,:)')/scaling.dFTtilde;
                end
                Qs_opt_nsc{i}(k,:) = [sol_val.QsQdots_unsc{i}(1:2:end,k) sol_val.QsQdots_col_unsc{i}(1:2:end,(k-1)*d+1:(k-1)*d+d)]*C(:,j+1); %QsQdots_sol already in the original scale
                eq_constr_dyn_Qsopt{i}((k-1)*d+j,:)=(h{i}*sol_val.QsQdots_col_unsc{i}(2:2:end,(k-1)*d+j)-Qs_opt_nsc{i}(k,:)')/scaling.QsQdots(2);

                Qdots_opt_nsc{i}(k,:) = [sol_val.QsQdots_unsc{i}(2:2:end,k) sol_val.QsQdots_col_unsc{i}(2:2:end,(k-1)*d+1:(k-1)*d+d)]*C(:,j+1); %QsQdots_sol already in the original scale
                eq_constr_dyn_Qdot{i}((k-1)*d+j,:)=(h{i}*sol_val.Qd2dot_col_unsc{i}(:,(k-1)*d+j)-Qdots_opt_nsc{i}(k,:)')/scaling.qd2dot(1);
            end


            QsQdots_kj_opt=[sol_val.QsQdots{i}(:,k), sol_val.QsQdots_col{i}(:,(k-1)*d+1:k*d)];
            eq_constr_cont_QsQdots{i}(k,:)=(sol_val.QsQdots{i}(:,k+1)-QsQdots_kj_opt*D)/10; %scaled

        end
    end
end


function  expdata=LoadData(N,d,tau_root,main_folder)
    current_folder=pwd;
    kinfiles=dir([main_folder '/kinematics/' '/*.mot']);
    kinfiles=kinfiles(~contains({kinfiles.name}, '_2Dangles'));
    forcefiles=dir([main_folder '/perturbation/' '/*.mot']);
    for i=1:length(kinfiles);
        kinfilename=kinfiles(i).name; 
        kindata=importdata([kinfiles(i).folder '/' kinfiles(i).name]);
        if strcmp(main_folder,'DataSeptember')
            kindata.data(:,1)=kindata.data(:,1);
        end
       kindata=ProcessKinematics(kindata);
       C=strrep(kinfilename,'perturb','');
       sufix=strrep(C,'.mot','');
        
       if strcmp(main_folder,'DataSeptember')
            trial_name=[strrep(kinfilename,'.mot','')];
       elseif strcmp(main_folder,'DataNovember')||strcmp(main_folder,'DataDecember')||contains(main_folder,'DataMarch2025')||contains(main_folder,'DataMay2025')
            trial_name=[strrep(strrep(kinfilename,'.mot',''),'.','_')];
            trial_name=[strrep(trial_name,'-','m')];
       else
        keyboard; %need to update
        trial_name=[trials{triali} '_' strrep(kinfilename,'.mot','')];
       end
       %parameterize with splines
       t=kindata.data(1,1):0.0002:kindata.data(end,1);
       expdata.(trial_name).kinematics(:,1)=t;
       expdata.(trial_name).kinematics_v(:,1)=t;
       expdata.(trial_name).kinematics_a(:,1)=t;
       [B,A]=butter(3,100/(5000/2));
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
       if contains(main_folder,'November')||contains(main_folder,'December')||contains(main_folder,'DataMarch2025')
            forcedata=importdata([forcefiles(i).folder '/' strrep(kinfiles(i).name,'kinematics','motor')]);
       elseif contains(main_folder,'September')
            forcedata=importdata([forcefiles(i).folder '/' strrep(kinfiles(i).name,'_kin','')]);
       else
            forcedata=importdata([forcefiles(i).folder '/' forcefiles(i).name]);
       end
       expdata.(trial_name).forces=forcedata.data;
       expdata.(trial_name).forces_labels=forcedata.colheaders;

       t0=max(forcedata.data(1,1),expdata.(trial_name).kinematics(1,1));
       tf=min(forcedata.data(end,1),expdata.(trial_name).kinematics(end,1));
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
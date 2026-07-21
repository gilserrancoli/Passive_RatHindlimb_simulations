clear all;
filepath=fileparts(matlab.desktop.editor.getActiveFilename);
cd(filepath);

import casadi.*
method='legendre';

%% Options
Options.useRigidTendon=0;
Options.optimizelM0=1;
Options.optimizelTs=1;
Options.optimizefiberdamping=0;
Options.optimizetendondamping=0;
Options.dofs_to_track=[1 0 0 1 1 0 0]; %1. hip flex 2. hip add 3. hip int 4. knee flex 5. ankle flex 6. ankle add 7. ankle int
Options.useOptimizedIG=0;
% Options.trialstotrack='all'; %if 'all', all trials are taken
Options.penalizeHighFTtilde=1;
Options.penalizeFTtot=0;
Options.penalizeoutoflMtilde1=0;
Options.optInertiaParam=0;
Options.minPassiveProp=1;
Options.pendevPassive=1; 
Options.useRestingMoments=1;
Options.optInertiapassiveParam=0;

Options.optimizeMuscleProp=0;
Options.optimizePassiveJointEl=1;
    Options.orderPassiveJoint=1; %either 1 or 3
    Options.KDwithinteractionterms=1;
        Options.nointeraction_hipankle=1;
Options.individualpassiveprop=0; %consider different stiffness and damper parameters for each case (each perturbation)
    Options.samepassiveprop=1; %only if previous is 0, all stiffness and damping parameters have the same values (ideal for testing variable inertia parameters for each perturbation)

Options.removefirstpertRat22=1;
Options.normalizeCostFunction=0;

Options.secondfolder=1;
main_folder='Datarat25\baseline_forward_2mm\';
main_folder2='Datarat25\baseline_backward_2mm\';
if Options.normalizeCostFunction
    ref_solution=load('sol_val_optJointPassive_Datarat25baselineFB_tolem5_Jminres777_JpenKDT1e-3_JminKDT1e-5_JminInertiaP1e-2_nokneemarker2D_pert2_wKDinteractionsnoanklehip_ParamID.mat');
else
    ref_solution=[];
end
Options.main_folder=main_folder;

%% 
N=40; %just used to discretize data
d=3;
[tau_root,C,D,B] = CollocationScheme(d,method);
Options.N=N;
Options.d=d;

% Load expdata in .mot
% trials={'h8'};
% expdata=LoadData_May(N,d,tau_root,trials);
expdata=LoadData(N,d,tau_root,main_folder);
if Options.secondfolder==1;
    expdata2=LoadData(N,d,tau_root,main_folder2);

    fieldnames2=fieldnames(expdata2);
    for i=1:length(fieldnames2)
        expdata.([fieldnames2{i} '_back'])=expdata2.([fieldnames2{i}]);
    end
end

if Options.removefirstpertRat22
    fieldnames_struct=fieldnames(expdata);
    I=find(contains(fieldnames_struct,'m20'));
    expdata=rmfield(expdata,fieldnames_struct(I));
end

%% load .trc files
current_folder=pwd;
cd([main_folder '\kinematics']);
names_trc=dir('*.trc');
for i=1:length(names_trc)
    name_trc=strrep(names_trc(i).name,'-','m');
    name_trc=strrep(name_trc,'.trc','');
    name_trc=strrep(name_trc,'.','_');
    all_trc_data.(name_trc)=readtable(names_trc(i).name,'FileType','text');
end
cd(current_folder);
cd([main_folder2 '\kinematics']);
if Options.secondfolder==1;
    names_trc=dir('*.trc');
    for i=1:length(names_trc)
        name_trc=strrep(names_trc(i).name,'-','m');
        name_trc=strrep(name_trc,'.trc','');
        name_trc=strrep(name_trc,'.','_');
        all_trc_data.([name_trc '_back'])=readtable(names_trc(i).name,'FileType','text');
    end
end
cd(current_folder);

% Prepare muscle-related data
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
        FM          = SX(nmuscles,1); %muscle fiber forces
        lMtilde     = SX.sym('lMtilde',nmuscles); % Normalized fiber lengths
        FMactFL     = SX.sym('FMactFL',nmuscles); % Normalized force-length term
        FMactFV     = SX.sym('FMactFV',nmuscles); % Normalized force-velocity term
    else
        %not possible to use compliant-tendon model without using Direct
        %Collocation here
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
    
    MTparam(:,Imuscles_cut)=[];

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
    
    
    cd('MuscleModel');
    if Options.useRigidTendon
        [FT, FM, lMtilde, FMactFL, FMactFV, FMpas, cos_alpha] = ...
                HillModel_RigidTendon(a',lMT',vMT',MTparam,lM0,lTs,fiber_damping);
        %
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
        % if (Options.optimizelM0)&&(Options.optimizelTs)
        %     f_forceEquilibrium_FtildeState = ...
        %         Function('f_forceEquilibrium_FtildeState',{a,FTtilde,dFTtilde,...
        %         lMT,vMT,lTs,lM0,fiber_damping,tendon_damping},{Hilldiff,FT,Fce,Fiso,vMmax,lMtilde,lTtilde});
        % elseif (~Options.optimizelM0)&&(Options.optimizelTs)
        %     f_forceEquilibrium_FtildeState = ...
        %         Function('f_forceEquilibrium_FtildeState',{a,FTtilde,dFTtilde,...
        %         lMT,vMT,lTs,fiber_damping,tendon_damping},{Hilldiff,FT,Fce,Fiso,vMmax,lMtilde,lTtilde});
        % else
        % end
    end
end

%% Load external function
% if Options.optInertiaParam
    if contains(main_folder,'rat21')
        F = external('F','RightRatHindlimb_Zhong_InertiaVar_rat21.dll')
    elseif contains(main_folder,'rat22')
        F = external('F','RightRatHindlimb_Zhong_InertiaVar_rat22.dll')
    elseif contains(main_folder,'rat23')
        F = external('F','RightRatHindlimb_Zhong_InertiaVar_rat23.dll')
    elseif contains(main_folder,'rat24')
        F = external('F','RightRatHindlimb_Zhong_InertiaVar_rat24.dll')
    elseif contains(main_folder,'rat25')
        F = external('F','RightRatHindlimb_Zhong_InertiaVar_rat25.dll')
    else
        keyboard;
    end
% else
%     keyboard;
%     F = external('F','RightRatHindlimb_Zhong.dll');   
% end
nPoints = Options.N * Options.d;
Fmap=F.map(nPoints,'serial');

ndofs=7; %ndofs to construct polynomials must be 7 here
if Options.optimizeMuscleProp
    %% Prepare polynomials
    % Indices of the muscles actuating the different joints for later use
    pathpolynomial = [pwd,'/Polynomials'];
    addpath(genpath(pathpolynomial));
    load([pathpolynomial,'/muscle_spanning_joint_INFO_subject.mat']);
    load([pathpolynomial,'/MuscleInfo_subject.mat']);
    pathmusclefunctions=[pwd,'/MuscleModel'];
    addpath(genpath(pathmusclefunctions));
    [~,mai] = MomentArmIndices_3D(muscle_names,...
        muscle_spanning_joint_INFO);


    pathcasadi_functions=[pwd,'/VariousFunctions/'];
    addpath(pathcasadi_functions);
    NMuscle_pol=nmuscles;
    CasADi_functions;
end

%% Formulation of the optimization problem
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
scaling.res=0.02;
scaling.Kstiff=0.1;
scaling.Ddamp=0.1;
scaling.theta=0.1;
scaling.inertiaparam=repmat([0.05 1e-5 0.05],1,3);
scaling.Inertiapassparam=0.01;
Options.scaling=scaling;

W.lM0lit=0.00001;
W.min_maxa=0.1;
W.state_for_reg=1e-6;
W.mindstate=0.01; %0.1
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
W.penalizeFTtot=1e-4;
W.penalizeoutoflMtilde1=0.01;
W.minresidual=50;
Options.W=W;


ParallelMode='thread';
NThreads=8;

%% Define bounds
if Options.optimizeMuscleProp
    load('lMTmax.mat');
    bounds.lTs.lower=(1e-4)*ones(nmuscles,1)/scaling.lTs;
    bounds.lTs.upper=lMTmax'/scaling.lTs;
    bounds.lM0.lower=(lMTmax'/4)/scaling.lM0;
    bounds.lM0.upper=ones(nmuscles,1)*0.06/scaling.lM0;
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
if Options.orderPassiveJoint==1
    if Options.KDwithinteractionterms
        if Options.nointeraction_hipankle
            bounds.K.lower=[0 0 0 -inf -inf]';
            bounds.D.lower=[0 0 0 -inf -inf]';
        else
            bounds.K.lower=[0 0 0 -inf -inf -inf]';
            bounds.D.lower=[0 0 0 -inf -inf -inf]';
        end
    else
        bounds.K.lower=[0 0 0];
        bounds.D.lower=[0 0 0];
    end
    bounds.K.upper=10;
    bounds.D.upper=10;
    bounds.Inertiapassparam.lower=0;
    bounds.Inertiapassparam.upper=10;
    if Options.useRestingMoments
        bounds.K0.lower=-inf*ones(3,1);
        bounds.K0.upper=inf*ones(3,1);
    else
        bounds_theta0.lower=[0.1174,-0.2571,0.1223,-0.9237,0.2676,0.3305,0.0909]'-180*pi/180;
        bounds.theta0.lower=bounds_theta0.lower(find(Options.dofs_to_track));
        bounds_theta0.upper=[0.1174,-0.2571,0.1223,-0.9237,0.2676,0.3305,0.0909]'+180*pi/180;
        bounds.theta0.upper=bounds_theta0.upper(find(Options.dofs_to_track));
    end
end
bounds.inertiaParam.lower=[0.005 1.e-7 -0.020 0.002 1e-8 0.010 0.0005 1e-9 -0.001]./scaling.inertiaparam; %mfem Izfem yfem mtib Iztib ytib mfoot Izfoot yfoot 
bounds.inertiaParam.upper=[0.025 1.e-5 -0.005 0.015 1e-6 0.020 0.0050 1e-7      0]./scaling.inertiaparam; %mfem Izfem yfem mtib Iztib ytib mfoot Izfoot yfoot 
bounds.hip_ankleqs.lower=[-10;  -20]*pi/180;
bounds.hip_ankleqs.upper=[ 10;   20]*pi/180;

%% Define guesses
if Options.optimizeMuscleProp
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
end
if Options.optimizePassiveJointEl
    if Options.orderPassiveJoint==1
        if Options.KDwithinteractionterms
            if Options.nointeraction_hipankle==1
                n_interactions=2;
            else
                n_interactions=3;
            end
            guess.Kstiff=[1e-5*ones(3,1); zeros(n_interactions,1)]/scaling.Kstiff;
            guess.Ddamp=[1e-5*ones(3,1); zeros(n_interactions,1)]/scaling.Ddamp;
            guess.Inertiapassparam=zeros(3+n_interactions,1)/scaling.Inertiapassparam;
        else
            guess.Kstiff=1e-5*ones(3,1)/scaling.Kstiff;
            guess.Ddamp=1e-5*ones(3,1)/scaling.Ddamp;
            guess.Inertiapassparam=zeros(3,1)/scaling.Inertiapassparam;
        end
    elseif Options.orderPassiveJoint==3
        guess.Kstiff=zeros(4*sum(Options.dofs_to_track),1);
        guess.Ddamp=zeros(3*sum(Options.dofs_to_track),1);
    end
end
if Options.useRestingMoments
    guess.K0=zeros(3,1);
else
    % guess.theta0=[0.0135, -0.1651,-0.0637,-1.1571,-0.1768]'; % computed as the mean of all initial values, TO BE REcalculated
    guess_theta0=[0.1174,-0.2571,0.1223,-0.9237,0.2676,0.3305,0.0909]'/scaling.theta;; % computed as the mean of all initial values for h8 perturb1
    guess.theta0=guess_theta0(find(Options.dofs_to_track)); %take only the dofs that are tracked
end
guess.inertiaParam=[0.01351 1.086e-06 -0.014936 0.00538 8.204e-07 0.0152275 0.00193 3.727e-08 -0.00546174]./scaling.inertiaparam;
guess.hip_ankleqs=[0;0]; %hip add hip rot ankle rot



%% Start with an empty optimization
x0=[];
lb=[];
ub=[];
varnames=[];
J=0;

if Options.optimizeMuscleProp
    %to be defined, optimize lT, lM0, fiber_damping, tendon_damping...
    keyboard;
end

if Options.optInertiaParam
    x0=[x0; guess.inertiaParam'];
    lb=[lb; bounds.inertiaParam.lower'];
    ub=[ub; bounds.inertiaParam.upper'];
    varnames=[varnames; repmat({'inertiaParam'},size(guess.inertiaParam'))];
else
    if contains(main_folder,'rat21')
        sol_val=load('sol_val_optJointPassive_Datarat21baselineFB_tolem5_Jminres777_JpenKDT1e-3_JminKDT1e-5_JminInertiaP1e-2_nokneemarker2D_pert2_wKDinteractionsnoanklehip_ParamID.mat');
    elseif contains(main_folder,'rat22')
        sol_val=load('sol_val_optJointPassive_Datarat22baselineFB_tolem5_Jminres777_JpenKDT1e-3_JminKDT1e-5_JminInertiaP1e-2_nokneemarker2D_pert2_wKDinteractionsnoanklehip_ParamID.mat');
    elseif contains(main_folder,'rat23')
        sol_val=load('sol_val_optJointPassive_Datarat23baselineFB_tolem5_Jminres777_JpenKDT1e-3_JminKDT1e-5_JminInertiaP1e-2_nokneemarker2D_pert2_wKDinteractionsnoanklehip_ParamID.mat');
    elseif contains(main_folder,'rat24')
        sol_val=load('sol_val_optJointPassive_Datarat24baselineFB_tolem5_Jminres777_JpenKDT1e-3_JminKDT1e-5_JminInertiaP1e-2_nokneemarker2D_pert2_wKDinteractionsnoanklehip_ParamID.mat');
    elseif contains(main_folder,'rat25')
        sol_val=load('sol_val_optJointPassive_Datarat25baselineFB_tolem5_Jminres777_JpenKDT1e-3_JminKDT1e-5_JminInertiaP1e-2_nokneemarker2D_pert2_wKDinteractionsnoanklehip_ParamID.mat');
    else
    end
    inertiaParam_unsc=sol_val.sol_val.inertiaParam_unsc;
    Options.inertiaParam_unsc=inertiaParam_unsc;
end

if Options.optimizePassiveJointEl
    if Options.individualpassiveprop==1
        ntrials4passprop=size(fieldnames(expdata),1);
    elseif (Options.individualpassiveprop==0)&&(Options.samepassiveprop==1)
        ntrials4passprop=size(fieldnames(expdata),1);
        list4passprop=ones(1,size(fieldnames(expdata),1));
    elseif contains(main_folder,'May2025')&&Options.secondfolder==0
        ntrials4passprop=4;
        list4passprop=[4 3 2 1 1 2 3 4]; %perturbations at 30 20 10 0 0 10...
    elseif contains(main_folder,'May2025')&&Options.secondfolder==1
        ntrials4passprop=4;
        list4passprop=[4 3 2 1 1 2 3 4 4 3 2 1 1 2 3 4]; %perturbations at 30 20 10 0 0 10...
    elseif contains(main_folder,'rat22')&&Options.secondfolder==0
        if Options.removefirstpertRat22
            ntrials4passprop=5;
            list4passprop=[5 4 3 2 1 1 2 3 4 5]; %perturbations at 30 20 10 0 -10 10 0 10...
        else
            ntrials4passprop=6;
            list4passprop=[6 5 4 3 2 1 1 2 3 4 5 6]; %perturbations at 30 20 10 0 -10 -20 -20 -10 0 10...
        end
    elseif contains(main_folder,'rat22')&&Options.secondfolder==1
        if Options.removefirstpertRat22
            ntrials4passprop=5;
            list4passprop=[5 4 3 2 1 1 2 3 4 5 5 4 3 2 1 1 2 3 4 5]; %perturbations at 30 20 10 0 -10 -20 -20 -10 0 10...
        else
            ntrials4passprop=6;
            list4passprop=[6 5 4 3 2 1 1 2 3 4 5 6 6 5 4 3 2 1 1 2 3 4 5 6]; %perturbations at 30 20 10 0 -10 -20 -20 -10 0 10...
        end
    elseif (contains(main_folder,'August2025')||contains(main_folder,'Datarat25')||contains(main_folder,'Datarat23'))&&Options.secondfolder==0
        ntrials4passprop=5;
        list4passprop=[5 4 3 2 1 1 2 3 4 5]; %perturbations at 30 20 10 0 -10 -10 0 10...
    elseif (contains(main_folder,'August2025')||contains(main_folder,'Datarat25')||contains(main_folder,'Datarat23'))&&Options.secondfolder==1
        ntrials4passprop=5;
        list4passprop=[5 4 3 2 1 1 2 3 4 5 5 4 3 2 1 1 2 3 4 5]; %perturbations at 30 20 10 0 -10 -10 0 10...
    else
        %define how many different passive prop are considered
        keyboard;
    end
    Options.ntrials4passprop=ntrials4passprop;
    Options.list4passprop=list4passprop;

    nvar_passparam=length(guess.Kstiff);
    Options.nvar_passparam=nvar_passparam;
    x0=[x0; repmat(guess.Kstiff,ntrials4passprop,1)];
    lb=[lb; repmat(bounds.K.lower,ntrials4passprop,1)];
    ub=[ub; repmat(bounds.K.upper,ntrials4passprop*nvar_passparam,1)];
    varnames=[varnames; repmat({'Ks'},ntrials4passprop*nvar_passparam,1)];
    
    x0=[x0; repmat(guess.Ddamp,ntrials4passprop,1)];
    lb=[lb; repmat(bounds.D.lower,ntrials4passprop,1)];
    ub=[ub; repmat(bounds.D.upper,ntrials4passprop*nvar_passparam,1)];
    varnames=[varnames; repmat({'Ds'},ntrials4passprop*nvar_passparam,1)];

    if Options.useRestingMoments
        x0=[x0; repmat(guess.K0,ntrials4passprop,1)];
        lb=[lb; repmat(bounds.K0.lower,ntrials4passprop,1)];
        ub=[ub; repmat(bounds.K0.upper,ntrials4passprop,1)];
        varnames=[varnames; repmat({'K0s'},ntrials4passprop*3,1)];
    else
        x0=[x0; repmat(guess.theta0,ntrials4passprop,1)];
        lb=[lb; repmat(bounds.theta0.lower,ntrials4passprop,1)];
        ub=[ub; repmat(bounds.theta0.upper,ntrials4passprop,1)];
        varnames=[varnames; repmat({'theta0s'},ntrials4passprop*3,1)];
    end

    if Options.optInertiapassiveParam
        x0=[x0; repmat(guess.Inertiapassparam,ntrials4passprop,1)];
        lb=[lb; repmat(bounds.Inertiapassparam.lower,ntrials4passprop*3,1)];
        ub=[ub; repmat(bounds.Inertiapassparam.upper,ntrials4passprop*3,1)];
        varnames=[varnames; repmat({'inertiapassprop'},ntrials4passprop*3,1)];
    end

end

x0=[x0; guess.hip_ankleqs];
lb=[lb; bounds.hip_ankleqs.lower];
ub=[ub; bounds.hip_ankleqs.upper];
varnames=[varnames; repmat({'hip_ankle_qs'},length(guess.hip_ankleqs),1)];

idx_var.inertiaParam = contains(varnames,'inertiaParam');
idx_var.Kstiff       = contains(varnames,'Ks');
idx_var.Ddamp        = contains(varnames,'Ds');
idx_var.K0           = contains(varnames,'K0s');
idx_var.theta0       = contains(varnames,'theta0s');
idx_var.inertiaPass  = contains(varnames,'inertiapassprop');
idx_var.hipAnkle     = contains(varnames,'hip_ankle_qs');

aux=costfun(x0,guess,idx_var, varnames,Options,expdata,Fmap,ref_solution,all_trc_data);
Jpattern_name=['Jpattern' num2str(length(aux)) 'x' num2str(length(x0)) '.mat'];
if exist(Jpattern_name)
    load(Jpattern_name);
else
    Jpattern = estimateJacobianPattern( ...
        @(x)costfun(x,guess,idx_var, varnames,Options,expdata,Fmap,ref_solution,all_trc_data), ...
        x0(:), ...
        lb(:), ...
        ub(:));
    save(Jpattern_name,'Jpattern');
end

A=[];
b=[];
Aeq=[];
beq=[];
% options=optimset('Display','iter','UseParallel',true,...
%     'Jacobian', 'off', 'JacobPattern', sparse(double(Jpattern)), ...
%     'FinDiffType','forward');
options = optimoptions('lsqnonlin', ...
'Display','iter-detailed', ...
'Algorithm','interior-point', ...
'UseParallel',true, ...
'SpecifyObjectiveGradient',false, ...
'JacobPattern',sparse(double(Jpattern)), ...
'FiniteDifferenceType','forward');
fprintf('Start optimization');
[x,resnorm,residual,exitflag,output] = lsqnonlin(@(x)costfun(x, guess, idx_var, varnames, Options, expdata,Fmap,ref_solution,all_trc_data),x0,lb,ub,A,b,Aeq,beq,@(x)nonlcon(x,varnames,Options),options);
fprintf('Optimization finished');

if Options.optInertiaParam
    I=contains(varnames,'inertiaParam');
    inertiaParam=x(I);
    inertiaParam_unsc=inertiaParam.*scaling.inertiaparam';
else
    inertiaParam_unsc=Options.inertiaParam_unsc;
end
if Options.optimizeMuscleProp&Options.optimizelM0
    %need to write, lM0 deviation
end
if Options.optimizePassiveJointEl
    I=contains(varnames,'Ks');
    Kstiff=x(I);
    I=contains(varnames,'Ds');
    Ddamp=x(I);
    if Options.useRestingMoments
        I=contains(varnames,'K0s');
        K0=x(I);
    else
        I=contains(varnames,'theta0s');
        theta0=x(I);
    end
    if Options.optInertiapassiveParam
        I=contains(varnames,'inertiapassprop');
        InertiapassiveParam=x(I);
    end

    Kstiff_unsc=Kstiff*scaling.Kstiff;
    Ddamp_unsc=Ddamp*scaling.Ddamp;
    if Options.useRestingMoments
        K0_unsc=K0*scaling.Kstiff;
    else
        theta0_unsc=theta0*scaling.theta;
    end
    if Options.optInertiapassiveParam
        InertiapassiveParam_unsc=InertiapassiveParam*scaling.Inertiapassparam;
    end
end
I=contains(varnames,'hip_ankle_qs');
hip_ankle_qs=x(I);

nametrials=fieldnames(expdata);
t_col_grid=1:N*(d+1);
t_col_grid(1:4:end)=[];
for i=1:size(nametrials,1)

    t0=expdata.(nametrials{i}).q(1,1);
    tf=expdata.(nametrials{i}).q(end,1);
    
    tgrid{i}=[];
    tgrid{i}(1:(d+1):(N+1)*(d+1))=t0:((tf-t0)/N):tf;
    deltat=tgrid{i}(d+1+1)-tgrid{i}(1);
    tgrid{i}(2:(d+1):N*(d+1))=tgrid{i}(1:(d+1):N*(d+1))+tau_root(2)*deltat;
    tgrid{i}(3:(d+1):N*(d+1))=tgrid{i}(1:(d+1):N*(d+1))+tau_root(3)*deltat;
    tgrid{i}(4:(d+1):N*(d+1))=tgrid{i}(1:(d+1):N*(d+1))+tau_root(4)*deltat;
    sol_val.t2plot{i}=t_col_grid;
    
    QsQdot_prescribed{i}(:,1:2:7*2)=expdata.(nametrials{i}).q(:,2:8); %pelvis dofs and sacroiliac_flx are prescribed
    QsQdot_prescribed{i}(:,2:2:7*2)=expdata.(nametrials{i}).qdot(:,2:8);
    Qd2dot_prescribed{i}=expdata.(nametrials{i}).qd2dot(:,2:8);
    QsQdots{i}(:,1:2:ndofs*2)=expdata.(nametrials{i}).q(1:(d+1):end,9:15);
    QsQdots_col{i}(:,1:2:ndofs*2)=expdata.(nametrials{i}).q(t_col_grid,9:15);
    QsQdots{i}(:,2:2:ndofs*2)=expdata.(nametrials{i}).qdot(1:(d+1):end,9:15);
    QsQdots_col{i}(:,2:2:ndofs*2)=expdata.(nametrials{i}).qdot(t_col_grid,9:15);
    Qd2dot_col{i}=expdata.(nametrials{i}).qd2dot(t_col_grid,9:15);
    forces_prescribed{i}=expdata.(nametrials{i}).f(t_col_grid,2:end);


    if Options.optimizeMuscleProp
        %define a, lM0, lTs...
    end

    for k=1:N
        for j=1:d
            if Options.optimizeMuscleProp
                %Get moment arms and muscle-tendon lengths at that frame
                keyboard;
            end
            if Options.KDwithinteractionterms
                Itrial4passprop=list4passprop(i);
                Kstiff_unsci6=Kstiff_unsc((Itrial4passprop-1)*nvar_passparam+1:Itrial4passprop*nvar_passparam);
                Ddamp_unsci6=Ddamp_unsc((Itrial4passprop-1)*nvar_passparam+1:Itrial4passprop*nvar_passparam);
                if Options.nointeraction_hipankle
                    Kstiff_unsci=[Kstiff_unsci6(1) Kstiff_unsci6(4) 0; Kstiff_unsci6(4) Kstiff_unsci6(2) Kstiff_unsci6(5); 0 Kstiff_unsci6(5) Kstiff_unsci6(3)];
                    Ddamp_unsci=[  Ddamp_unsci6(1)  Ddamp_unsci6(4) 0;  Ddamp_unsci6(4)  Ddamp_unsci6(2)  Ddamp_unsci6(5); 0  Ddamp_unsci6(5)  Ddamp_unsci6(3)];
                else
                    Kstiff_unsci=[Kstiff_unsci6(1) Kstiff_unsci6(4) Kstiff_unsci6(5); Kstiff_unsci6(4) Kstiff_unsci6(2) Kstiff_unsci6(6); Kstiff_unsci6(5) Kstiff_unsci6(6) Kstiff_unsci6(3)];
                    Ddamp_unsci=[Ddamp_unsci6(1) Ddamp_unsci6(4) Ddamp_unsci6(5); Ddamp_unsci6(4) Ddamp_unsci6(2) Ddamp_unsci6(6); Ddamp_unsci6(5) Ddamp_unsci6(6) Ddamp_unsci6(3)];
                end
                if Options.useRestingMoments
                    K0_unsci=K0_unsc((Itrial4passprop-1)*3+1:Itrial4passprop*3);
                else
                    theta0_unsci=theta0_unsc((Itrial4passprop-1)*3+1:Itrial4passprop*3);
                end
            end
            all_QsQdot(1,1:14)=QsQdot_prescribed{i}((k-1)*(d+1)+j+1,:);
            all_QsQdot(1,15:28)=QsQdots_col{i}((k-1)*d+j,:);

            all_Qd2dot(1:7)=Qd2dot_prescribed{i}((k-1)*(d+1)+j+1,:);
            all_Qd2dot(8:14)=Qd2dot_col{i}((k-1)*d+j,:);

            forces_prescribed_j=forces_prescribed{i}((k-1)*d+j,:);

            % if Options.optInertiaParam
                out{i}((k-1)*d+j,:)=full(F([all_QsQdot';all_Qd2dot';forces_prescribed_j';inertiaParam_unsc]));
            % else
            %     out{i}((k-1)*d+j,:)=full(F([all_QsQdot';all_Qd2dot';forces_prescribed_j']));
            % end

            if Options.optimizeMuscleProp
                %Hill equilibrium
                keyboard;
            end
        
            %moment equilibrium
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
                    if Options.individualpassiveprop==0;
                        Itrial4passprop=list4passprop(i);
                    else
                        Itrial4passprop=i;
                    end
                    if Options.KDwithinteractionterms
                        %this is valid only if using sagittal plane angles
                        %and moments
                        if Options.useRestingMoments
                            PassiveM_hip_flx=K0_unsci(1)-Kstiff_unsci(1,1:3)*(QsQdots_col{i}((k-1)*d+j,[1 7 9])')-Ddamp_unsci(1,1:3)*QsQdots_col{i}((k-1)*d+j,[2 8 10])';
                        else
                            PassiveM_hip_flx=-Kstiff_unsci(1,1:3)*(QsQdots_col{i}((k-1)*d+j,[1 7 9])'-theta0_unsci)-Ddamp_unsci(1,1:3)*QsQdots_col{i}((k-1)*d+j,[2 8 10])';
                        end
                    else
                        if Options.useRestingMoments
                            PassiveM_hip_flx=K0_unsci(1)-Kstiff_unsc((Itrial4passprop-1)*3+1)*(QsQdots_col{i}((k-1)*d+j,1))-Ddamp_unsc((Itrial4passprop-1)*3+1)*QsQdots_col{i}((k-1)*d+j,2);
                        else
                            PassiveM_hip_flx=-Kstiff_unsc((Itrial4passprop-1)*3+1)*(QsQdots_col{i}((k-1)*d+j,1)-theta0_unsc((Itrial4passprop-1)*3+1))-Ddamp_unsc((Itrial4passprop-1)*3+1)*QsQdots_col{i}((k-1)*d+j,2);
                        end
                    end
                    if Options.optInertiapassiveParam
                        PassiveM_hip_flx=PassiveM_hip_flx-InertiapassiveParam_unsc((Itrial4passprop-1)*3+1)*Qd2dot_col{i}((k-1)*d+j,1);
                    end
                    
                else
                        PassiveM_hip_flx=0;
                end
                out_model{i}((k-1)*d+j,8)= T_hip_flx+PassiveM_hip_flx;
            end

            %hip adduction
            if Options.dofs_to_track(2)
                if Options.KDwithinteractionterms
                    keyboard;
                end
                if Options.optimizeMuscleProp
                    FT_hip_add=FT_j{j}(mai(2).mus);
                    T_hip_add=FT_hip_add'*MAj{j}.hip_add;
                else
                    T_hip_add=0;
                end
                I=sum(Options.dofs_to_track(1:2));
                if Options.optimizePassiveJointEl
                    if Options.individualpassiveprop==0;
                        Itrial4passprop=list4passprop(i);
                    else
                        Itrial4passprop=i;
                    end
                    PassiveM_hip_add=-Kstiff_unsc((Itrial4passprop-1)*3+1)*(QsQdots_col{i}((k-1)*d+j,3)-theta0_unsc((Itrial4passprop-1)*3+1))-Ddamp_unsc((Itrial4passprop-1)*3+1)*QsQdots_col{i}((k-1)*d+j,4);
                    if Options.optInertiapassiveParam
                        PassiveM_hip_add=PassiveM_hip_add-InertiapassiveParam_unsc((Itrial4passprop-1)*3+1)*Qd2dot_col{i}((k-1)*d+j,2);
                    end
                    
                else
                        PassiveM_hip_add=0;
                end
                out_model{i}((k-1)*d+j,9)=T_hip_add+PassiveM_hip_add;
            end

            %hip rotation
            if Options.dofs_to_track(3)
                if Options.KDwithinteractionterms
                    keyboard;
                end
                if Options.optimizeMuscleProp
                    FT_hip_rot=FT_j{j}(mai(3).mus);
                    T_hip_rot=FT_hip_rot'*MAj{j}.hip_rot;
                else
                    T_hip_rot=0;
                end
                I=sum(Options.dofs_to_track(1:3));
                if Options.optimizePassiveJointEl
                    if Options.individualpassiveprop==0;
                        Itrial4passprop=list4passprop(i);
                    else
                        Itrial4passprop=i;
                    end
                    PassiveM_hip_rot=-Kstiff_unsc((Itrial4passprop-1)*3+1)*(QsQdots_col{i}((k-1)*d+j,5)-theta0_unsc((Itrial4passprop-1)*3+1))-Ddamp_unsc((Itrial4passprop-1)*3+1)*QsQdots_col{i}((k-1)*d+j,6);
                    if Options.optInertiapassiveParam
                        PassiveM_hip_rot=PassiveM_hip_rot-InertiapassiveParam_unsc((Itrial4passprop-1)*3+1)*Qd2dot_col{i}((k-1)*d+j,3);
                    end
                    
                else
                        PassiveM_hip_rot=0;
                end
                out_model{i}((k-1)*d+j,10)=T_hip_rot+PassiveM_hip_rot;
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
                    if Options.individualpassiveprop==0;
                        Itrial4passprop=list4passprop(i);
                    else
                        Itrial4passprop=i;
                    end
                    if Options.KDwithinteractionterms
                        %this is valid only if using sagittal plane angles
                        %and moments
                        if Options.useRestingMoments
                            PassiveM_knee_flx=K0_unsci(2)-Kstiff_unsci(2,1:3)*(QsQdots_col{i}((k-1)*d+j,[1 7 9])')-Ddamp_unsci(2,1:3)*QsQdots_col{i}((k-1)*d+j,[2 8 10])';
                        else
                            PassiveM_knee_flx=-Kstiff_unsci(2,1:3)*(QsQdots_col{i}((k-1)*d+j,[1 7 9])'-theta0_unsci)-Ddamp_unsci(2,1:3)*QsQdots_col{i}((k-1)*d+j,[2 8 10])';
                        end
                    else
                        if Options.useRestingMoments
                            PassiveM_knee_flx=K0_unsci(2)-Kstiff_unsc((Itrial4passprop-1)*3+2)*(QsQdots_col{i}((k-1)*d+j,7))-Ddamp_unsc((Itrial4passprop-1)*3+2)*QsQdots_col{i}((k-1)*d+j,8);
                        else
                            PassiveM_knee_flx=-Kstiff_unsc((Itrial4passprop-1)*3+2)*(QsQdots_col{i}((k-1)*d+j,7)-theta0_unsc((Itrial4passprop-1)*3+2))-Ddamp_unsc((Itrial4passprop-1)*3+2)*QsQdots_col{i}((k-1)*d+j,8);
                        end
                    end
                    if Options.optInertiapassiveParam
                        PassiveM_knee_flx=PassiveM_knee_flx-InertiapassiveParam_unsc((Itrial4passprop-1)*3+2)*Qd2dot_col{i}((k-1)*d+j,4);
                    end
                    
                else
                        PassiveM_knee_flx=0;
                end
                out_model{i}((k-1)*d+j,11)=T_knee_flx+PassiveM_knee_flx;
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
                    if Options.individualpassiveprop==0;
                        Itrial4passprop=list4passprop(i);
                    else
                        Itrial4passprop=i;
                    end
                    if Options.KDwithinteractionterms
                        %this is valid only if using sagittal plane angles
                        %and moments
                        if Options.useRestingMoments
                            PassiveM_ankle_flx=K0_unsci(3)-Kstiff_unsci(3,1:3)*(QsQdots_col{i}((k-1)*d+j,[1 7 9])')-Ddamp_unsci(3,1:3)*QsQdots_col{i}((k-1)*d+j,[2 8 10])';
                        else
                            PassiveM_ankle_flx=-Kstiff_unsci(3,1:3)*(QsQdots_col{i}((k-1)*d+j,[1 7 9])'-theta0_unsci)-Ddamp_unsci(3,1:3)*QsQdots_col{i}((k-1)*d+j,[2 8 10])';
                        end
                    else
                        if Options.useRestingMoments
                            PassiveM_ankle_flx=K0_unsci(3)-Kstiff_unsc((Itrial4passprop-1)*3+3)*(QsQdots_col{i}((k-1)*d+j,9))-Ddamp_unsc((Itrial4passprop-1)*3+3)*QsQdots_col{i}((k-1)*d+j,10);
                        else
                            PassiveM_ankle_flx=-Kstiff_unsc((Itrial4passprop-1)*3+3)*(QsQdots_col{i}((k-1)*d+j,9)-theta0_unsc((Itrial4passprop-1)*3+3))-Ddamp_unsc((Itrial4passprop-1)*3+3)*QsQdots_col{i}((k-1)*d+j,10);
                        end
                    end
                    if Options.optInertiapassiveParam
                        PassiveM_ankle_flx=PassiveM_ankle_flx-InertiapassiveParam_unsc((Itrial4passprop-1)*3+3)*Qd2dot_col{i}((k-1)*d+j,5);
                    end
                    
                else
                        PassiveM_ankle_flx=0;
                end
                out_model{i}((k-1)*d+j,12)=T_ankle_flx+PassiveM_ankle_flx;
            end

            %ankle adduction
            if Options.dofs_to_track(6)
                if Options.KDwithinteractionterms
                    keyboard;
                end
                if Options.optimizeMuscleProp
                    FT_ankle_add=FT_j{j}(mai(6).mus);
                    T_ankle_add=FT_ankle_add'*MAj{j}.ankle_add;
                else
                    T_ankle_add=0;
                end
                I=sum(Options.dofs_to_track(1:6));
                if Options.optimizePassiveJointEl
                    if Options.individualpassiveprop==0;
                        Itrial4passprop=list4passprop(i);
                    else
                        Itrial4passprop=i;
                    end
                    PassiveM_ankle_add=-Kstiff_unsc((Itrial4passprop-1)*3+3)*(QsQdots_col{i}((k-1)*d+j,11)-theta0_unsc((Itrial4passprop-1)*3+3))-Ddamp_unsc((Itrial4passprop-1)*3+3)*QsQdots_col{i}((k-1)*d+j,12);
                    if Options.optInertiapassiveParam
                        PassiveM_ankle_add=PassiveM_ankle_add-InertiapassiveParam_unsc((Itrial4passprop-1)*3+3)*Qd2dot_col{i}((k-1)*d+j,6);
                    end
                    
                else
                        PassiveM_ankle_add=0;
                end
                out_model{i}((k-1)*d+j,13)=T_ankle_add+PassiveM_ankle_add;
            end

            %ankle rotation
            if Options.dofs_to_track(7)
                if Options.KDwithinteractionterms
                    keyboard;
                end
                if Options.optimizeMuscleProp
                    FT_ankle_int=FT_j{j}(mai(7).mus);
                    T_ankle_int=FT_ankle_int'*MAj{j}.ankle_int;
                else
                    T_ankle_int=0;
                end
                I=sum(Options.dofs_to_track(1:7));
                if Options.optimizePassiveJointEl
                    if Options.individualpassiveprop==0;
                        Itrial4passprop=list4passprop(i);
                    else
                        Itrial4passprop=i;
                    end
                    PassiveM_ankle_int=-Kstiff_unsc((Itrial4passprop-1)*3+3)*(QsQdots_col{i}((k-1)*d+j,13)-theta0_unsc((Itrial4passprop-1)*3+3))-Ddamp_unsc((Itrial4passprop-1)*3+3)*QsQdots_col{i}((k-1)*d+j,14);
                    if Options.optInertiapassiveParam
                        PassiveM_ankle_int=PassiveM_ankle_int-InertiapassiveParam_unsc((Itrial4passprop-1)*3+3)*Qd2dot_col{i}((k-1)*d+j,7);
                    end
                    
                else
                        PassiveM_ankle_int=0;
                end
                out_model{i}((k-1)*d+j,14)=T_ankle_int+PassiveM_ankle_int;
            end

            if Options.optimizeMuscleProp
                %penalize high tendon forces and deviation of lMtilde...
            end

        end
    end

end

if Options.optInertiapassiveParam
    sol_val.InertiapassiveParam=InertiapassiveParam;
    sol_val.InertiapassiveParam_unsc=InertiapassiveParam_unsc;
end
sol_val.Kstiff=Kstiff;
sol_val.Kstiff_unsc=Kstiff_unsc;
sol_val.Ddamp=Ddamp;
sol_val.Ddamp_unsc=Ddamp_unsc;
if Options.useRestingMoments
    sol_val.K0=K0;
    sol_val.K0_unsc=K0_unsc;
else
    sol_val.theta0=theta0;
    sol_val.theta0_unsc=theta0_unsc;
end
sol_val.out=out;
sol_val.out_model=out_model;
sol_val.Options=Options;
sol_val.guess=guess;
sol_val.bounds=bounds;
sol_val.scaling=scaling;
sol_val.W=W;
if Options.optInertiaParam
    sol_val.inertiaParam=inertiaParam;
end
sol_val.inertiaParam_unsc=inertiaParam_unsc;
sol_val.tgrid=tgrid;
sol_val.nametrials=nametrials;
sol_val.forces_prescribed=forces_prescribed;
sol_val.main_folder=main_folder;
sol_val.QsQdot_prescribed=QsQdot_prescribed;
sol_val.QsQdots_col_unsc=QsQdots_col;
sol_val.Qd2dot_col_unsc=Qd2dot_col;
sol_val.Qd2dot_prescribed=Qd2dot_prescribed;
sol_val.hip_ankle_qs=hip_ankle_qs;

save('solution.mat','sol_val','x','Options','resnorm','residual','exitflag','output');


function f=costfun(x,guess,idx_var,varnames,Options,expdata,Fmap,ref_solution,all_trc_data)
f=[];

ndofs=7;
scaling=Options.scaling;
W=Options.W;
ntrials4passprop=Options.ntrials4passprop;
d=Options.d;
N=Options.N;
t_col_grid=1:N*(d+1);
t_col_grid(1:4:end)=[];
list4passprop=Options.list4passprop;
nvar_passparam=length(guess.Kstiff);

nTrials = numel(fieldnames(expdata));
nTrackedDofs = nnz(Options.dofs_to_track);

nDynamicTerms = ...
    nTrials * Options.N * Options.d * nTrackedDofs;

nRegularizationTerms = 0;

if Options.optInertiaParam
    nRegularizationTerms = nRegularizationTerms + ...
        nnz(contains(varnames,'inertiaParam'));
end

if Options.optimizePassiveJointEl && Options.minPassiveProp == 1
    nRegularizationTerms = nRegularizationTerms + ...
        nnz(contains(varnames,'Ks'));

    nRegularizationTerms = nRegularizationTerms + ...
        nnz(contains(varnames,'Ds'));

    if Options.optInertiapassiveParam
        nRegularizationTerms = nRegularizationTerms + ...
            nnz(contains(varnames,'inertiapassprop'));
    end

    if ~Options.pendevPassive && ~Options.useRestingMoments
        nRegularizationTerms = nRegularizationTerms + ...
            numel(guess.theta0) * Options.ntrials4passprop;
    end

    if Options.useRestingMoments
        nRegularizationTerms = nRegularizationTerms + ...
            nnz(contains(varnames,'K0s'));
    end
end

nHipAnkleTerms = nnz(contains(varnames,'hip_ankle_qs'));

nResiduals = ...
    nRegularizationTerms + ...
    nHipAnkleTerms + ...
    nDynamicTerms;

f = zeros(nResiduals,1);
idx = 0;



if Options.optInertiaParam
    % I=contains(varnames,'inertiaParam');
    inertiaParam=x(idx_var.inertiaParam);

    % f=[f; inertiaParam-guess.inertiaParam'];
    nv=numel(inertiaParam);
    f(idx+1:idx+nv)=inertiaParam-guess.inertiaParam';
    idx=idx+nv;

    inertiaParam_unsc=inertiaParam.*scaling.inertiaparam';
else
    inertiaParam_unsc=Options.inertiaParam_unsc;
end
if Options.optimizeMuscleProp&Options.optimizelM0
    %need to write, lM0 deviation
end
if Options.optimizePassiveJointEl
    % I=contains(varnames,'Ks');
    Kstiff=x(idx_var.Kstiff);
    % I=contains(varnames,'Ds');
    Ddamp=x(idx_var.Ddamp);
    if Options.useRestingMoments
        % I=contains(varnames,'K0s');
        K0=x(idx_var.K0);
    else
        % I=contains(varnames,'theta0s');
        theta0=x(idx_var.theta0);
    end
    if Options.optInertiapassiveParam
        % I=contains(varnames,'inertiapassprop');
        InertiapassiveParam=x(idx_var.InertiaPass);
    end
    
    if Options.minPassiveProp==1
        % f=[f; W.Kstiff*Kstiff(:)];
        nv=numel(Kstiff(:));
        f(idx+1:idx+nv)=W.Kstiff*Kstiff(:);
        idx = idx +nv; 

        % f=[f; W.Ddamp*Ddamp(:)];
        nv=numel(Ddamp(:));
        f(idx+1:idx+nv)=W.Ddamp*Ddamp(:);
        idx=idx+nv;

        if Options.optInertiapassiveParam
            % f=[f; W.Inertiapassparam*InertiapassiveParam];
            nv=numel(InertiapassiveParam);
            f(idx+1:idx+nv)=W.Inertiapassparam*InertiapassiveParam;
            idx=idx+nv;
        end
        if Options.pendevPassive
            %do not include twice the penalization of theta0
        else
            difftheta0=theta0-repmat(guess.theta0,ntrials4passprop,1);
            % f=[f; W.theta0*difftheta0(:)];
            nv=numel(difftheta0(:));
            f(idx+1:idx+nv)=W.theta0*difftheta0(:);
            idx=idx+nv;
        end
        if Options.useRestingMoments
            % f=[f; W.Kstiff*K0(:)];
            nv=numel(K0(:));
            f(idx+1:idx+nv)=W.Kstiff*K0(:);
            idx=idx+nv;
        end
    end
    Kstiff_unsc=Kstiff*scaling.Kstiff;
    Ddamp_unsc=Ddamp*scaling.Ddamp;
    if Options.useRestingMoments
        K0_unsc=K0*scaling.Kstiff;
    else
        theta0_unsc=theta0*scaling.theta;
    end
    if Options.optInertiapassiveParam
        InertiapassiveParam_unsc=InertiapassiveParam*scaling.Inertiapassparam;
    end
end
% I=contains(varnames,'hip_ankle_qs');
hip_ankle_qs=x(idx_var.hipAnkle);
% f=[f; hip_ankle_qs];
nv=numel(hip_ankle_qs);
f(idx+1:idx+nv)=hip_ankle_qs;
idx=idx+nv;

nametrials=fieldnames(expdata);
trc_nametrials=fieldnames(all_trc_data);
for i=1:size(nametrials,1)
    trc_data=all_trc_data.(strrep(nametrials{i},'_2D',''));
    [QsQdots_new,Qd2dots_new]=ComputeIK(expdata.(nametrials{i}).q,hip_ankle_qs,trc_data,Options,nametrials{i},expdata.(nametrials{i}).q(:,1));
    QsQdot_prescribed{i}(:,1:2:7*2)=expdata.(nametrials{i}).q(:,2:8); %pelvis dofs and sacroiliac_flx are prescribed
    QsQdot_prescribed{i}(:,2:2:7*2)=expdata.(nametrials{i}).qdot(:,2:8);
    Qd2dot_prescribed{i}=expdata.(nametrials{i}).qd2dot(:,2:8);

    % QsQdots{i}(:,1:2:ndofs*2)=expdata.(nametrials{i}).q(1:(d+1):end,9:15);
    % QsQdots_col{i}(:,1:2:ndofs*2)=expdata.(nametrials{i}).q(t_col_grid,9:15); 
    % QsQdots{i}(:,2:2:ndofs*2)=expdata.(nametrials{i}).qdot(1:(d+1):end,9:15);
    % QsQdots_col{i}(:,2:2:ndofs*2)=expdata.(nametrials{i}).qdot(t_col_grid,9:15);
    % Qd2dot_col{i}=expdata.(nametrials{i}).qd2dot(t_col_grid,9:15);
    QsQdots{i}=QsQdots_new(1:(d+1):end,:);
    QsQdots_col{i}=QsQdots_new(t_col_grid,:);
    Qd2dot_col{i}=Qd2dots_new(t_col_grid,:);
    forces_prescribed{i}=expdata.(nametrials{i}).f(t_col_grid,2:end);

    % ---------------------------------------------------------------
    % Batch all CasADi inputs for this trial. Each column is one
    % collocation point, ordered as p = (k-1)*d + j. Fmap must be
    % created once before lsqnonlin using:
    %     Fmap = F.map(Options.N*Options.d,'serial');
    % ---------------------------------------------------------------
    nPoints = N*d;

    all_QsQdot_batch = [ ...
        QsQdot_prescribed{i}(t_col_grid,:), ...
        QsQdots_col{i}];                         % nPoints x 28

    all_Qd2dot_batch = [ ...
        Qd2dot_prescribed{i}(t_col_grid,:), ...
        Qd2dot_col{i}];                          % nPoints x 14

    inertiaParam_batch = repmat( ...
        inertiaParam_unsc(:),1,nPoints);

    Finput_batch = [ ...
        all_QsQdot_batch.'; ...
        all_Qd2dot_batch.'; ...
        forces_prescribed{i}.'; ...
        inertiaParam_batch];

    % One MATLAB-to-CasADi call for all collocation points in the trial.
    out_batch = full(Fmap(Finput_batch));

    if size(out_batch,2) ~= nPoints
        error('costfun_batched:UnexpectedCasadiOutput', ...
            ['Expected Fmap to return one output column per collocation ' ...
             'point (%d), but received %d columns.'], ...
            nPoints,size(out_batch,2));
    end

    if Options.optimizeMuscleProp
        %define a, lM0, lTs...
    end
    if Options.optimizePassiveJointEl
        if Options.KDwithinteractionterms
            Itrial4passprop=list4passprop(i);
            try
                Kstiff_unsci6=Kstiff_unsc((Itrial4passprop-1)*nvar_passparam+1:Itrial4passprop*nvar_passparam);
                Ddamp_unsci6=Ddamp_unsc((Itrial4passprop-1)*nvar_passparam+1:Itrial4passprop*nvar_passparam);
            catch
                keyboard;
            end
            if Options.nointeraction_hipankle
                Kstiff_unsci=[Kstiff_unsci6(1)  Kstiff_unsci6(4)                0; Kstiff_unsci6(4) Kstiff_unsci6(2) Kstiff_unsci6(5); 0                Kstiff_unsci6(5) Kstiff_unsci6(3)];
                Ddamp_unsci= [ Ddamp_unsci6(1)   Ddamp_unsci6(4)                0;  Ddamp_unsci6(4)  Ddamp_unsci6(2)  Ddamp_unsci6(5); 0                Ddamp_unsci6(5)  Ddamp_unsci6(3)];
            else
                Kstiff_unsci=[Kstiff_unsci6(1) Kstiff_unsci6(4) Kstiff_unsci6(5); Kstiff_unsci6(4) Kstiff_unsci6(2) Kstiff_unsci6(6); Kstiff_unsci6(5) Kstiff_unsci6(6) Kstiff_unsci6(3)];
                Ddamp_unsci=[Ddamp_unsci6(1) Ddamp_unsci6(4) Ddamp_unsci6(5); Ddamp_unsci6(4) Ddamp_unsci6(2) Ddamp_unsci6(6); Ddamp_unsci6(5) Ddamp_unsci6(6) Ddamp_unsci6(3)];
            end
            if Options.useRestingMoments
                K0_unsci=K0_unsc((Itrial4passprop-1)*3+1:Itrial4passprop*3);
            else
                theta0_unsci=theta0_unsc((Itrial4passprop-1)*3+1:Itrial4passprop*3);
            end
        end
    end

    for k=1:N
        for j=1:d
            p=(k-1)*d+j;
            if Options.optimizeMuscleProp
                %Get moment arms and muscle-tendon lengths at that frame
                keyboard;
            end
            all_QsQdot(1,1:14)=QsQdot_prescribed{i}((k-1)*(d+1)+j+1,:);
            all_QsQdot(1,15:28)=QsQdots_col{i}((k-1)*d+j,:);

            all_Qd2dot(1:7)=Qd2dot_prescribed{i}((k-1)*(d+1)+j+1,:);
            all_Qd2dot(8:14)=Qd2dot_col{i}((k-1)*d+j,:);

            forces_prescribed_j=forces_prescribed{i}((k-1)*d+j,:);

            % % if Options.optInertiaParam
                % out=full(F([all_QsQdot';all_Qd2dot';forces_prescribed_j';inertiaParam_unsc]));
            % % else
            % %     out=full(F([all_QsQdot';all_Qd2dot';forces_prescribed_j']));
            % % end
            
            % Retrieve the already-computed CasADi output for this point.
            out=out_batch(:,p);

            if Options.optimizeMuscleProp
                %Hill equilibrium
                keyboard;
            end
        
            %moment equilibrium
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
                    if Options.individualpassiveprop==0;
                        Itrial4passprop=list4passprop(i);
                    else
                        Itrial4passprop=i;
                    end
                    if Options.KDwithinteractionterms
                        %this is valid only if using sagittal plane angles
                        %and moments
                        if Options.useRestingMoments
                            PassiveM_hip_flx=K0_unsci(1)-Kstiff_unsci(1,1:3)*(QsQdots_col{i}((k-1)*d+j,[1 7 9])')-Ddamp_unsci(1,1:3)*QsQdots_col{i}((k-1)*d+j,[2 8 10])';
                        else
                            PassiveM_hip_flx=-Kstiff_unsci(1,1:3)*(QsQdots_col{i}((k-1)*d+j,[1 7 9])'-theta0_unsci)-Ddamp_unsci(1,1:3)*QsQdots_col{i}((k-1)*d+j,[2 8 10])';
                        end
                    else
                        if Options.useRestingMoments
                            PassiveM_hip_flx=K0_unsci(1)-Kstiff_unsc((Itrial4passprop-1)*3+1)*(QsQdots_col{i}((k-1)*d+j,1))-Ddamp_unsc((Itrial4passprop-1)*3+1)*QsQdots_col{i}((k-1)*d+j,2);
                        else
                            PassiveM_hip_flx=-Kstiff_unsc((Itrial4passprop-1)*3+1)*(QsQdots_col{i}((k-1)*d+j,1)-theta0_unsc((Itrial4passprop-1)*3+1))-Ddamp_unsc((Itrial4passprop-1)*3+1)*QsQdots_col{i}((k-1)*d+j,2);
                        end
                    end
                    if Options.optInertiapassiveParam
                        PassiveM_hip_flx=PassiveM_hip_flx-InertiapassiveParam_unsc((Itrial4passprop-1)*3+1)*Qd2dot_col{i}((k-1)*d+j,1);
                    end
                    
                else
                        PassiveM_hip_flx=0;
                end
                if Options.normalizeCostFunction
                    vari=var(ref_solution.sol_val.out{i}(:,8));
                    % f=[f; (out(8)-T_hip_flx-PassiveM_hip_flx)/(vari*1e5)];
                    nv=1;
                    f(idx+1:idx+nv)=(out(8)-T_hip_flx-PassiveM_hip_flx)/(vari*1e5);
                    idx=idx+nv;
                else
                    % f=[f; (out(8)-T_hip_flx-PassiveM_hip_flx)/scaling.T];
                    nv=1;
                    f(idx+1:idx+nv)=(out(8)-T_hip_flx-PassiveM_hip_flx)/scaling.T;
                    idx=idx+nv;
                end
            end

            %hip adduction
            if Options.dofs_to_track(2)
                if Options.KDwithinteractionterms
                    keyboard;
                end
                if Options.optimizeMuscleProp
                    FT_hip_add=FT_j{j}(mai(2).mus);
                    T_hip_add=FT_hip_add'*MAj{j}.hip_add;
                else
                    T_hip_add=0;
                end
                I=sum(Options.dofs_to_track(1:2));
                if Options.optimizePassiveJointEl
                    if Options.individualpassiveprop==0;
                        Itrial4passprop=list4passprop(i);
                    else
                        Itrial4passprop=i;
                    end
                    PassiveM_hip_add=-Kstiff_unsc((Itrial4passprop-1)*3+1)*(QsQdots_col{i}((k-1)*d+j,3)-theta0_unsc((Itrial4passprop-1)*3+1))-Ddamp_unsc((Itrial4passprop-1)*3+1)*QsQdots_col{i}((k-1)*d+j,4);
                    if Options.optInertiapassiveParam
                        PassiveM_hip_add=PassiveM_hip_add-InertiapassiveParam_unsc((Itrial4passprop-1)*3+1)*Qd2dot_col{i}((k-1)*d+j,2);
                    end
                    
                else
                        PassiveM_hip_add=0;
                end
                if Options.normalizeCostFunction
                    vari=var(ref_solution.sol_val.out{i}(:,9));
                    % f=[f; (out(9)-T_hip_add-PassiveM_hip_add)/(vari*1e5)];
                    nv=1;
                    f(idx+1:idx+nv)=(out(9)-T_hip_add-PassiveM_hip_add)/(vari*1e5);
                    idx=idx+nv;
                else
                    % f=[f; (out(9)-T_hip_add-PassiveM_hip_add)/scaling.T];
                    nv=1;
                    f(idx+1:idx+nv)=(out(9)-T_hip_add-PassiveM_hip_add)/scaling.T;
                    idx=idx+nv;
                end
            end

            %hip rotation
            if Options.dofs_to_track(3)
                if Options.KDwithinteractionterms
                    keyboard;
                end
                if Options.optimizeMuscleProp
                    FT_hip_rot=FT_j{j}(mai(3).mus);
                    T_hip_rot=FT_hip_rot'*MAj{j}.hip_rot;
                else
                    T_hip_rot=0;
                end
                I=sum(Options.dofs_to_track(1:3));
                if Options.optimizePassiveJointEl
                    if Options.individualpassiveprop==0;
                        Itrial4passprop=list4passprop(i);
                    else
                        Itrial4passprop=i;
                    end
                    PassiveM_hip_rot=-Kstiff_unsc((Itrial4passprop-1)*3+1)*(QsQdots_col{i}((k-1)*d+j,5)-theta0_unsc((Itrial4passprop-1)*3+1))-Ddamp_unsc((Itrial4passprop-1)*3+1)*QsQdots_col{i}((k-1)*d+j,6);
                    if Options.optInertiapassiveParam
                        PassiveM_hip_rot=PassiveM_hip_rot-InertiapassiveParam_unsc((Itrial4passprop-1)*3+1)*Qd2dot_col{i}((k-1)*d+j,3);
                    end
                    
                else
                        PassiveM_hip_rot=0;
                end
                if Options.normalizeCostFunction
                    vari=var(ref_solution.sol_val.out{i}(:,10));
                    % f=[f; (out(10)-T_hip_rot-PassiveM_hip_rot)/(vari*1e5)];
                    nv=1;
                    f(idx+1:idx+nv)=(out(10)-T_hip_rot-PassiveM_hip_rot)/(vari*1e5);
                    idx=idx+nv;
                else
                    % f=[f; (out(10)-T_hip_rot-PassiveM_hip_rot)/scaling.T];
                    nv=1;
                    f(idx+1:idx+nv)=(out(10)-T_hip_rot-PassiveM_hip_rot)/scaling.T;
                    idx=idx+nv;
                end
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
                    if Options.individualpassiveprop==0;
                        Itrial4passprop=list4passprop(i);
                    else
                        Itrial4passprop=i;
                    end
                    if Options.KDwithinteractionterms
                        %this is valid only if using sagittal plane angles
                        %and moments
                        if Options.useRestingMoments
                            PassiveM_knee_flx=K0_unsci(2)-Kstiff_unsci(2,1:3)*(QsQdots_col{i}((k-1)*d+j,[1 7 9])')-Ddamp_unsci(2,1:3)*QsQdots_col{i}((k-1)*d+j,[2 8 10])';
                        else
                            PassiveM_knee_flx=-Kstiff_unsci(2,1:3)*(QsQdots_col{i}((k-1)*d+j,[1 7 9])'-theta0_unsci)-Ddamp_unsci(2,1:3)*QsQdots_col{i}((k-1)*d+j,[2 8 10])';
                        end
                    else
                        if Options.useRestingMoments
                            PassiveM_knee_flx=K0_unsci(2)-Kstiff_unsc((Itrial4passprop-1)*3+2)*(QsQdots_col{i}((k-1)*d+j,7))-Ddamp_unsc((Itrial4passprop-1)*3+2)*QsQdots_col{i}((k-1)*d+j,8);
                        else
                            PassiveM_knee_flx=-Kstiff_unsc((Itrial4passprop-1)*3+2)*(QsQdots_col{i}((k-1)*d+j,7)-theta0_unsc((Itrial4passprop-1)*3+2))-Ddamp_unsc((Itrial4passprop-1)*3+2)*QsQdots_col{i}((k-1)*d+j,8);
                        end
                    end
                    if Options.optInertiapassiveParam
                        PassiveM_knee_flx=PassiveM_knee_flx-InertiapassiveParam_unsc((Itrial4passprop-1)*3+2)*Qd2dot_col{i}((k-1)*d+j,4);
                    end
                    
                else
                        PassiveM_knee_flx=0;
                end
                if Options.normalizeCostFunction
                    vari=var(ref_solution.sol_val.out{i}(:,11));
                    % f=[f; (out(11)-T_knee_flx-PassiveM_knee_flx)/(vari*1e5)];
                    nv=1;
                    f(idx+1:idx+nv)=(out(11)-T_knee_flx-PassiveM_knee_flx)/(vari*1e5);
                    idx=idx+nv;
                else
                    % f=[f; (out(11)-T_knee_flx-PassiveM_knee_flx)/scaling.T];
                    nv=1;
                    f(idx+1:idx+nv)=(out(11)-T_knee_flx-PassiveM_knee_flx)/scaling.T;
                    idx=idx+nv;
                end
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
                    if Options.individualpassiveprop==0;
                        Itrial4passprop=list4passprop(i);
                    else
                        Itrial4passprop=i;
                    end
                    if Options.KDwithinteractionterms
                        %this is valid only if using sagittal plane angles
                        %and moments
                        if Options.useRestingMoments
                            PassiveM_ankle_flx=K0_unsci(3)-Kstiff_unsci(3,1:3)*(QsQdots_col{i}((k-1)*d+j,[1 7 9])')-Ddamp_unsci(3,1:3)*QsQdots_col{i}((k-1)*d+j,[2 8 10])';
                        else
                            PassiveM_ankle_flx=-Kstiff_unsci(3,1:3)*(QsQdots_col{i}((k-1)*d+j,[1 7 9])'-theta0_unsci)-Ddamp_unsci(3,1:3)*QsQdots_col{i}((k-1)*d+j,[2 8 10])';
                        end
                    else
                        if Options.useRestingMoments
                            PassiveM_ankle_flx=K0_unsci(3)-Kstiff_unsc((Itrial4passprop-1)*3+3)*(QsQdots_col{i}((k-1)*d+j,9))-Ddamp_unsc((Itrial4passprop-1)*3+3)*QsQdots_col{i}((k-1)*d+j,10);
                        else
                            PassiveM_ankle_flx=-Kstiff_unsc((Itrial4passprop-1)*3+3)*(QsQdots_col{i}((k-1)*d+j,9)-theta0_unsc((Itrial4passprop-1)*3+3))-Ddamp_unsc((Itrial4passprop-1)*3+3)*QsQdots_col{i}((k-1)*d+j,10);
                        end
                    end
                    if Options.optInertiapassiveParam
                        PassiveM_ankle_flx=PassiveM_ankle_flx-InertiapassiveParam_unsc((Itrial4passprop-1)*3+3)*Qd2dot_col{i}((k-1)*d+j,5);
                    end
                    
                else
                        PassiveM_ankle_flx=0;
                end
                if Options.normalizeCostFunction
                    vari=var(ref_solution.sol_val.out{i}(:,12));
                    % f=[f; (out(12)-T_ankle_flx-PassiveM_ankle_flx)/(vari*1e5)];
                    nv=1;
                    f(idx+1:idx+nv)=(out(12)-T_ankle_flx-PassiveM_ankle_flx)/(vari*1e5);
                    idx=idx+nv;
                else
                    % f=[f; (out(12)-T_ankle_flx-PassiveM_ankle_flx)/scaling.T];
                    nv=1;
                    f(idx+1:idx+nv)=(out(12)-T_ankle_flx-PassiveM_ankle_flx)/scaling.T;
                    idx=idx+nv;
                end
            end

            %ankle adduction
            if Options.dofs_to_track(6)
                if Options.KDwithinteractionterms
                    keyboard;
                end
                if Options.optimizeMuscleProp
                    FT_ankle_add=FT_j{j}(mai(6).mus);
                    T_ankle_add=FT_ankle_add'*MAj{j}.ankle_add;
                else
                    T_ankle_add=0;
                end
                I=sum(Options.dofs_to_track(1:6));
                if Options.optimizePassiveJointEl
                    if Options.individualpassiveprop==0;
                        Itrial4passprop=list4passprop(i);
                    else
                        Itrial4passprop=i;
                    end
                    PassiveM_ankle_add=-Kstiff_unsc((Itrial4passprop-1)*3+3)*(QsQdots_col{i}((k-1)*d+j,11)-theta0_unsc((Itrial4passprop-1)*3+3))-Ddamp_unsc((Itrial4passprop-1)*3+3)*QsQdots_col{i}((k-1)*d+j,12);
                    if Options.optInertiapassiveParam
                        PassiveM_ankle_add=PassiveM_ankle_add-InertiapassiveParam_unsc((Itrial4passprop-1)*3+3)*Qd2dot_col{i}((k-1)*d+j,6);
                    end
                    
                else
                        PassiveM_ankle_add=0;
                end
                if Options.normalizeCostFunction
                    vari=var(ref_solution.sol_val.out{i}(:,13));
                    % f=[f; (out(13)-T_ankle_add-PassiveM_ankle_add)/(vari*1e5)];
                    nv=1;
                    f(idx+1:idx+nv)=(out(13)-T_ankle_add-PassiveM_ankle_add)/(vari*1e5);
                    idx=idx+nv;
                else
                    % f=[f; (out(13)-T_ankle_add-PassiveM_ankle_add)/scaling.T];
                    nv=1;
                    f(idx+1:idx+nv)=(out(13)-T_ankle_add-PassiveM_ankle_add)/scaling.T;
                    idx=idx+nv;
                end
            end

            %ankle rotation
            if Options.dofs_to_track(7)
                if Options.KDwithinteractionterms
                    keyboard;
                end
                if Options.optimizeMuscleProp
                    FT_ankle_int=FT_j{j}(mai(7).mus);
                    T_ankle_int=FT_ankle_int'*MAj{j}.ankle_int;
                else
                    T_ankle_int=0;
                end
                I=sum(Options.dofs_to_track(1:7));
                if Options.optimizePassiveJointEl
                    if Options.individualpassiveprop==0;
                        Itrial4passprop=list4passprop(i);
                    else
                        Itrial4passprop=i;
                    end
                    PassiveM_ankle_int=-Kstiff_unsc((Itrial4passprop-1)*3+3)*(QsQdots_col{i}((k-1)*d+j,13)-theta0_unsc((Itrial4passprop-1)*3+3))-Ddamp_unsc((Itrial4passprop-1)*3+3)*QsQdots_col{i}((k-1)*d+j,14);
                    if Options.optInertiapassiveParam
                        PassiveM_ankle_int=PassiveM_ankle_int-InertiapassiveParam_unsc((Itrial4passprop-1)*3+3)*Qd2dot_col{i}((k-1)*d+j,7);
                    end
                    
                else
                        PassiveM_ankle_int=0;
                end
                if Options.normalizeCostFunction
                    vari=var(ref_solution.sol_val.out{i}(:,14));
                    % f=[f; (out(14)-T_ankle_int-PassiveM_ankle_int)/(vari*1e5)];
                    nv=1;
                    f(idx+1*idx+nv)=(out(14)-T_ankle_int-PassiveM_ankle_int)/(vari*1e5);
                    idx=idx+nv;
                else
                    % f=[f; (out(14)-T_ankle_int-PassiveM_ankle_int)/scaling.T];
                    nv=1;
                    f(idx+1*idx+nv)=(out(14)-T_ankle_int-PassiveM_ankle_int)/scaling.T;
                    idx=idx+nv;
                end
            end

            if Options.optimizeMuscleProp
                %penalize high tendon forces and deviation of lMtilde...
            end

        end
    end

end




assert(idx == nResiduals, ...
    'Residual count mismatch: expected %d, wrote %d.', ...
    nResiduals,idx);    



end




function  expdata=LoadData(N,d,tau_root,main_folder)
    current_folder=pwd;
    kinfiles=dir([main_folder '/kinematics/' '/*.mot']);
    names = {kinfiles.name};
    nums = str2double(regexprep(names, '^kinematics_(\d+)_.*$', '$1'));
    [~, order] = sort(nums, 'ascend');
    kinfiles = kinfiles(order);

    if contains(main_folder,'nokneemarker')||contains(main_folder,'August2025')
        kinfiles=kinfiles(contains({kinfiles.name}, '_2D'));
    else
        kinfiles=kinfiles(~contains({kinfiles.name}, '_2Dangles'));
    end
    forcefiles=dir([main_folder '/perturbation/' '/*.mot']);
    names = {forcefiles.name};
    nums = str2double(regexprep(names, '^motor_(\d+)_.*$', '$1'));
    [~, order] = sort(nums, 'ascend');
    forcefiles = forcefiles(order);

        for i=1:length(kinfiles);
           kinfilename=kinfiles(i).name; 
           kindata=importdata([kinfiles(i).folder '/' kinfiles(i).name]);
           if strcmp(main_folder,'DataSeptember')
            trial_name=[strrep(kinfilename,'.mot','')];
           elseif strcmp(main_folder,'DataNovember')||strcmp(main_folder,'DataDecember')||contains(main_folder,'DataMarch2025')||contains(main_folder,'DataMay2025')||contains(main_folder,'August2025')||contains(main_folder,'Datarat25')||contains(main_folder,'Datarat23')||contains(main_folder,'Datarat22')
                trial_name=[strrep(strrep(kinfilename,'.mot',''),'.','_')];
                trial_name=[strrep(trial_name,'-','m')];
           else
            keyboard; %need to update
            trial_name=[trials{triali} '_' strrep(kinfilename,'.mot','')];
           end

           % kindata=ProcessKinematics(kindata);
           % C=strrep(kinfilename,'perturb','');
           % sufix=strrep(C,'.mot','');
            
           % trial_name=[strrep(strrep(kinfilename,'.mot',''),'.','_')];

           %parameterize with splines
           t=kindata.data(1,1):0.0002:kindata.data(end,1);
           expdata.(trial_name).kinematics(:,1)=t;
           expdata.(trial_name).kinematics_v(:,1)=t;
           expdata.(trial_name).kinematics_a(:,1)=t;
           [B,A]=butter(3,100/(5000/2));
           for j=2:size(kindata.data,2)
                intdata=interp1(kindata.data(:,1),kindata.data(:,j),t,'spline');
                % smoothed_kin=smooth(t,intdata,0.4,'rloess'); %test when using
                % IK from OpenSim
                smoothed_filt_kin=filtfilt(B,A,intdata);

                kindata_spline(j-1)=spline(t,smoothed_filt_kin);

                expdata.(trial_name).kinematics(:,j)=ppval(kindata_spline(j-1),t);
                expdata.(trial_name).kinematics_v(:,j)=ppval(fnder(kindata_spline(j-1),1),t);
                expdata.(trial_name).kinematics_a(:,j)=ppval(fnder(kindata_spline(j-1),2),t);

                % expdata.(trial_name).kinematics(:,j)=filtfilt(B,A,expdata.(trial_name).kinematics(:,j));
                % expdata.(trial_name).kinematics_v(:,j)=filtfilt(B,A,expdata.(trial_name).kinematics_v(:,j));
                % expdata.(trial_name).kinematics_a(:,j)=filtfilt(B,A,expdata.(trial_name).kinematics_a(:,j));

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
           % if strcmp(forcefiles(i).name,strrep(kinfiles(i).name,'_kin',''))
           % else
           %     keyboard;
           % end
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

function [c,ceq]=nonlcon(x,varnames,Options)

    ceq=[];
    c=[];
    ntrials4passprop=Options.ntrials4passprop;
    nvar_passparam =Options.nvar_passparam;
    I=contains(varnames,'Ks');
    Kstiff=x(I);
    I=contains(varnames,'Ds');
    Ddamp=x(I);
    for i=1:ntrials4passprop
        Kstiffi=Kstiff((i-1)*nvar_passparam+1:i*nvar_passparam);
        K11=Kstiffi(1);
        K22=Kstiffi(2);
        K33=Kstiffi(3);
        K12=Kstiffi(4);
        
        if Options.nointeraction_hipankle
            K13=0;
            K23=Kstiffi(5);
        else
            K13=Kstiffi(5);
            K23=Kstiffi(6);
        end
        c=[c; K12^2-K11*K22; K13^2-K11*K33; K23^2-K22*K33];
        c=[c; -(K11*K22*K33+2*K12*K13*K23-K11*(K23^2)-K22*(K13^2)-K33*(K12^2))]; %(detK >0)
    end

    

end


function [QsQdots,Qd2dots]=ComputeIK(expdata_q,hip_ankle_qs,trc_data,Options,nametrial,t)

    if contains(Options.main_folder,'DataMay2025') %rat 21
        load('MeanPelvisData_May2025.mat');
        lfemur=38; %mm
        ltibia=42; %mm
    elseif contains(Options.main_folder,'rat22') %rat 22
        load('MeanPelvisData_Rat22.mat');
        lfemur=36; %mm
        ltibia=43; %mm
    elseif contains(Options.main_folder,'rat23')
        load('MeanPelvisData_Rat23.mat');
        lfemur=36; %mm
        ltibia=45; %mm
    elseif contains(Options.main_folder,'August2025') %rat24
        load('MeanPelvisData_August2025');
        lfemur=38; %mm
        ltibia=42; %mm
    elseif contains(Options.main_folder,'rat25')
        load('MeanPelvisData_Rat25.mat');
        lfemur=36; %mm
        ltibia=46; %mm
    end
    
    %Obtain the number of perturbation
    strpos=strfind(nametrial,'pos');
    if nametrial(strpos+3)=='-'
        pert=nametrial(strpos+3:strpos+5);
    else
        pert=nametrial(strpos+3:strpos+4);
    end
    pert=strrep(pert,'_','');
    found=false;
    p=1;
    try
        while (p<=length(mdata_bypert))&(~found)
            found=str2num(pert)==mdata_bypert(p).perts;
            if found
                continue;
            end
            p=p+1;
        end
    catch
        keyboard;
    end

    I_ankle=find(contains(trc_data.Properties.VariableNames,'ankle'));
    I_mtp=find(contains(trc_data.Properties.VariableNames,'mtp'));
    I_toe=find(contains(trc_data.Properties.VariableNames,'toe'));

    nf=size(trc_data,1)-1;
    %position hip and ankle markers
    pos_pelvis_top=repmat(mdata_bypert(p).pelvis_top,nf,1);
    pos_pelvis_bottom=repmat(mdata_bypert(p).pelvis_bottom,nf,1);
    pos_hip=repmat(mdata_bypert(p).hip,nf,1);
    % pos_ankle_0=table2array(trc_data(2:end,I_ankle:(I_ankle+2)));
    pos_ankle_0i=mdata_bypert(p).ankle;
    pos_mtp_0i=mdata_bypert(p).mtp;
    pos_toe_0i=mdata_bypert(p).toe;
    [R_spine_ground, R_pelvis_ground, t_spine, spine_euler]=ComputePelvisConf(pos_pelvis_top,pos_pelvis_bottom,pos_hip,Options.main_folder);

    q.data(:,1)=0:0.0001:0.08;
    q.data(:,2:4)=repmat(spine_euler,size(q.data,1),1);
    q.data(:,5:7)=repmat(t_spine',size(q.data,1),1);
    q.data(:,8)=3.7;

    %create new ankle, mtp and toe positions
    path_split=split(Options.main_folder,'\');
    prefix=path_split{2};
    if contains(lower(prefix),'forward')
        perturbation_in=mdata_bypert(p).pert_forward;
    elseif contains(lower(prefix),'backward')
        perturbation_in=mdata_bypert(p).pert_backward;
    else
        keyboard;
    end
    v_ankle=mdata_bypert(p).v_ankle;

    if v_ankle(1)<0
        pos_ankle=pos_ankle_0i(1,:)-(perturbation_in-perturbation_in(1))*v_ankle';
    else
        pos_ankle=pos_ankle_0i(1,:)+(perturbation_in-perturbation_in(1))*v_ankle';
    end

    v_mtp=mdata_bypert(p).v_mtp;
    if v_mtp(1)<0
        pos_mtp=pos_mtp_0i(1,:)-(perturbation_in-perturbation_in(1))*v_mtp';
    else
        pos_mtp=pos_mtp_0i(1,:)+(perturbation_in-perturbation_in(1))*v_mtp';
    end

    v_toe=mdata_bypert(p).v_toe;
    if v_toe(1)<0
        pos_toe=pos_toe_0i(1,:)-(perturbation_in-perturbation_in(1))*v_toe';
    else
        pos_toe=pos_toe_0i(1,:)+(perturbation_in-perturbation_in(1))*v_toe';
    end

    solutionGuess=[];
    previousSolution=[];
    for i=1:nf
        % [hip_flexion(i,:), hip_introtation(i,:), pos_knee(i,:), Z_femur(i,:), solutionGuess] = computeHipFlexion_Kneepos_FromDistances(R_spine_ground, R_pelvis_ground, hip_ankle_qs, pos_hip(i,:), pos_pelvis_bottom(i,:), pos_pelvis_top(i,:), pos_ankle(i,:), pos_ankle(1,:), pos_ankle(end,:), lfemur, ltibia,solutionGuess);
        solution = solveHipAnglesAnalytical(R_pelvis_ground, ...
            hip_ankle_qs(1), pos_hip(i,:), pos_ankle(i,:), lfemur, ...
            ltibia, previousSolution);
        hip_flexion(i,1) = solution.hipFlexion;
        hip_introtation(i,1) = solution.hipInternalRotation;
        pos_knee(i,:) = solution.posKnee.';
        Z_femur(i,:) = solution.ZFemur.';
        previousSolution = [solution.hipFlexion; ...
            solution.hipInternalRotation];

        knee_angle_deg(i,:) = -computeKneeAngle(pos_hip(i,:), pos_knee(i,:), pos_ankle(i,:));
        
        ankle_flexion(i,:) = computeAnkleFlexionAngle(hip_ankle_qs, pos_hip(i,:), pos_knee(i,:), pos_ankle(i,:), pos_ankle(1,:), pos_ankle(end,:), pos_toe(i,:),Z_femur(i,:));

    end

    q.data(:,1)=0:0.0001:0.08;
    % q.data(:,2:4)=repmat(spine_euler,size(q.data,1),1);
    % q.data(:,5:7)=repmat(t_spine',size(q.data,1),1);
    % q.data(:,8)=3.7;
    % q.data(:,2:8)=expdata_q(:,2:8);
    % q.data(:,[9:15])=[hip_flexion repmat(hip_ankle_qs(1)*pi/180,size(q.data,1),1) hip_introtation*pi/180 knee_angle_deg ankle_flexion repmat([0 hip_ankle_qs(2)]*pi/180,size(q.data,1),1)];

    new_qs=[hip_flexion repmat(hip_ankle_qs(1),size(q.data,1),1) hip_introtation knee_angle_deg ankle_flexion repmat([0 hip_ankle_qs(2)],size(q.data,1),1)];

    time_trc=table2array(trc_data(2:end,2));
    [B,A]=butter(3,100/(5000/2));
    % smoothed_filt_kin=filtfilt(B,A,new_qs);
   for j=1:size(new_qs,2)
        % intdata=interp1(time_trc(:,1),new_qs(:,j),t,'spline');      

        kindata_spline(j)=spline(time_trc,new_qs(:,j));

        new_qs2.kinematics(:,j)=ppval(kindata_spline(j),t);
        new_qs2.kinematics_v(:,j)=ppval(fnder(kindata_spline(j),1),t);
        new_qs2.kinematics_a(:,j)=ppval(fnder(kindata_spline(j),2),t);
   end


    QsQdots=zeros(size(new_qs2.kinematics,1),size(new_qs2.kinematics,2)*2);
    QsQdots(:,1:2:end)=new_qs2.kinematics;
    QsQdots(:,2:2:end)=new_qs2.kinematics_v;
    Qd2dots=new_qs2.kinematics_a;
end

function [R_spine_ground, R_pelvis, t_spine, spine_euler]=ComputePelvisConf(pos_pelvis_top,pos_pelvis_bottom,pos_hip,in_folder)

    x_pelvis=pos_pelvis_top-pos_pelvis_bottom;
    x_pelvis=x_pelvis./sqrt(sum(x_pelvis.^2,2));

    vert_ground=[0 0 1];

    z_pelvis=cross(x_pelvis,repmat(vert_ground,size(x_pelvis,1),1));
    z_pelvis=z_pelvis./sqrt(sum(z_pelvis.^2,2));
    
    y_pelvis=cross(z_pelvis,x_pelvis);
    y_pelvis=y_pelvis./sqrt(sum(y_pelvis.^2,2));
    
R_pelvis=[x_pelvis(1,:)' y_pelvis(1,:)' z_pelvis(1,:)'];

R_sacroiliac=eul2rotm([3.7*pi/180,0,0],'ZXY');
if contains(in_folder,'May2025') %rat21
    pos_hip_in_pelvis_frame=[0 0 0.007]';
    t_sacroiliac=[0.00428703 -0.00257222 0.00857407]';
elseif contains(in_folder,'rat22')
    pos_hip_in_pelvis_frame=[0 0 0.007]';
    t_sacroiliac=[0.00465216 -0.0027913 0.00930433]';
elseif contains(in_folder,'rat23')
    pos_hip_in_pelvis_frame=[0 0 0.007]';
    t_sacroiliac=[0.00370292 -0.00222175 0.00740585]';
elseif contains(in_folder,'August2025') %rat24
    pos_hip_in_pelvis_frame=[0 0 0.00491681]';
    t_sacroiliac=[0.00409734 -0.0024584 0.00819468]';
elseif contains(in_folder,'rat25')
    pos_hip_in_pelvis_frame=[0 0 0.007]';
    t_sacroiliac=[0.00383521 -0.00230112 0.00767041]';
else
    keyboard;
end
R_spine_ground=R_pelvis*R_sacroiliac';

t_hip=pos_hip(1,:)'/1000-R_pelvis*pos_hip_in_pelvis_frame;
t_spine=t_hip-R_spine_ground*t_sacroiliac;

% spine_euler=rotm2eul(R_spine_ground,'ZXY')*180/pi;
spine_euler=rotm2eulZXY_fast(R_spine_ground)*180/pi;

end


function [hip_flexion, hip_introtation, pos_knee, Z_femur, nextGuess] = computeHipFlexion_Kneepos_FromDistances(R_spine_ground, R_pelvis_ground, hip_ankle_qs, pos_hip, pos_pelvis_bottom, pos_pelvis_top, pos_ankle, pos_ankle_1, pos_ankle_end, lfemur, ltibia,initialGuess);

% Hip-to-knee direction in the femur local frame.
    femurLocalAxis = [0; -1; 0];

hipAnkleDistance=norm(pos_ankle-pos_hip);
if hipAnkleDistance>lfemur+ltibia ||...
        hipAnkleDistance< abs(lfemur-ltibia)
    error('The measured ankle position is unreachable.');
end

hipAdduction=hip_ankle_qs(1);

residualFunction = @(q) reconstructionResiduals( ...
        q, pos_hip, pos_ankle, lfemur, ltibia, ...
        R_pelvis_ground, hipAdduction, femurLocalAxis);

options = optimoptions( ...
    'fsolve', ...
    'Display', 'off', ...
    'FunctionTolerance', 1e-12, ...
    'StepTolerance', 1e-12, ...
    'OptimalityTolerance', 1e-12);

if isempty(initialGuess)
    x0=[30*pi/180; 0];
else
    x0 = initialGuess;
end

[solution, residualVector, exitflag, solverOutput] = ...
    fsolve(residualFunction, x0, options);

if exitflag <= 0
    warning(['The hip reconstruction solver did not converge. ', ...
             'Exit flag: %d'], exitflag);
end

hipFlexionRad = solution(1);
hipInternalRotationRad = solution(2);
nextGuess=solution;

femurRotation = buildFemurRotation( ...
    R_pelvis_ground, hipFlexionRad, hipAdduction, hipInternalRotationRad);

pos_knee = pos_hip' ...
    + lfemur * femurRotation * femurLocalAxis;

residual = norm(pos_ankle' - pos_knee) - ltibia;

Z_femur=femurRotation(:,3);

hip_flexion=rad2deg(hipFlexionRad);
hip_introtation=rad2deg(hipInternalRotationRad);

end


function value = tibiaConstraint( ...
    hipFlexion, hipPosition, anklePosition, ...
    femurLength, tibiaLength, pelvisRotation, ...
    hipAdduction, hipInternalRotation, femurLocalAxis)

    femurRotation = buildFemurRotation( ...
        pelvisRotation, hipFlexion, hipAdduction, ...
        hipInternalRotation);

    kneePosition = hipPosition' ...
        + femurLength * femurRotation * femurLocalAxis;

    value = norm(anklePosition' - kneePosition)^2 ...
        - tibiaLength^2;
end

function R = buildFemurRotation( ...
    pelvisRotation, hipFlexion, ...
    hipAdduction, hipInternalRotation)
%BUILDFEMURROTATION Femur-local frame to global frame.
%
% Convention:
%   1. Flexion about the spine-fixed Z axis.
%   2. Adduction about the rotated local X axis.
%   3. Internal/external rotation about the femur-fixed local Y axis.
%
% Active rotations with column vectors:
%
%   R = R_spine * Rz(flexion) * Rx(adduction) * Ry(internalRotation)

    R = pelvisRotation ...
        * rotationZ(hipFlexion) ...
        * rotationX(hipAdduction) ...
        * rotationY(hipInternalRotation);
end

function R = rotationX(angle)
%ROTATIONX Active right-handed rotation about X.

    c = cos(angle);
    s = sin(angle);

    R = [1, 0,  0;
         0, c, -s;
         0, s,  c];
end


function R = rotationY(angle)
%ROTATIONY Active right-handed rotation about Y.

    c = cos(angle);
    s = sin(angle);

    R = [ c, 0, s;
          0, 1, 0;
         -s, 0, c];
end


function R = rotationZ(angle)
%ROTATIONZ Active right-handed rotation about Z.

    c = cos(angle);
    s = sin(angle);

    R = [c, -s, 0;
         s,  c, 0;
         0,  0, 1];
end

function ankle_flexion = computeAnkleFlexionAngle(hip_ankle_qs,pos_hip, pos_knee, pos_ankle, pos_ankle1, pos_ankle_end, P3, Z_femur)
    % Inputs:
    % Z_femur: global femur Z axis (shared with tibia and foot)

    ankleIntRad = hip_ankle_qs(2);

    %% Common Z axis
    Z_common = Z_femur/norm(Z_femur);

    %% 1. Tibia frame
    Y_tibia = (pos_knee - pos_ankle) / norm(pos_knee - pos_ankle);
    Z_tibia = Z_femur / norm(Z_femur);
    X_tibia = cross(Y_tibia, Z_tibia);
    X_tibia = X_tibia / norm(X_tibia);
    R_tibia = [X_tibia(:), Y_tibia(:), Z_tibia(:)];

    %% 2. Foot frame
   
    % Hip-knee-ankle plane
    femurVector = pos_knee - pos_hip;
    hipToAnkle  = pos_ankle - pos_hip;
    planeNormal = cross(femurVector, hipToAnkle);
    planeNormal = planeNormal/norm(planeNormal);
    distanceToPlane = dot(P3 - pos_ankle, planeNormal);

    footVector = pos_ankle_end - pos_ankle1;
    footVectorProjected = footVector ...
        - dot(footVector, planeNormal) * planeNormal;

    %Compute X_foot, Y_foot and Z_foot
    X_foot_tilde = footVectorProjected;
    X_foot_tilde = X_foot_tilde / norm(X_foot_tilde);
    
    ankle_flexion_aux=acos(dot(X_foot_tilde,Y_tibia));

    ankle_flexion = 2*pi/3 - ankle_flexion_aux; 

    %Construct the complete relative rotation
    % Adduction is explicitly zero.

    ankleAdductionRad = 0;

    R_relative = ...
        rotationZ(ankle_flexion) ...
        * rotationX(ankleAdductionRad) ...
        * rotationY(ankleIntRad);

    R_foot = R_tibia * R_relative;

    X_foot = R_foot(:,1);
    Y_foot = R_foot(:,2);
    Z_foot = R_foot(:,3);

    %% 3. Relative rotation and ZXY Euler angles
    % eulerZXY = rotm2eul(R_relative, 'ZXY');  % MATLAB intrinsic ZXY
    eulerZXY = rotm2eulZXY_fast (R_relative);
    ankle_flexion=eulerZXY(1);

end

function knee_angle = computeKneeAngle(P1, P4, P2)
    % Inputs:
    % P1 - Hip position [1x3]
    % P4 - Knee position [1x3]
    % P2 - Ankle position [1x3]

    % Step 1: Define Y axes
    Y_femur = (P1 - P4) / norm(P1 - P4);  % from knee to hip
    Y_tibia = (P4 - P2) / norm(P4 - P2);  % from ankle to knee

    % Step 2: Compute angle between Y axes
    cos_theta = dot(Y_femur, Y_tibia);  % dot product
    cos_theta = min(max(cos_theta, -1), 1);  % clamp for numerical safety
    knee_angle = acos(cos_theta);  % angle in radians
end

function residuals = reconstructionResiduals( ...
    q, pos_hip, pos_ankle, lfemur, ltibia, ...
    R_pelvis_ground, hipAdduction, femurLocalAxis)

    hipFlexion = q(1);
    hipInternalRotation = q(2);

    femurRotation = buildFemurRotation( ...
        R_pelvis_ground, ...
        hipFlexion, ...
        hipAdduction, ...
        hipInternalRotation);

    pos_knee = pos_hip' ...
        + lfemur * femurRotation * femurLocalAxis;

    Z_femur = femurRotation(:,3);
    Z_femur = Z_femur / norm(Z_femur);

    tibiaVector = pos_ankle' - pos_knee;

    % Constraint 1: correct tibia length.
    lengthResidual = norm(tibiaVector) - ltibia;

    % Constraint 2: tibia must be perpendicular to Z_femur.
    planeResidual = dot(Z_femur, tibiaVector);

    residuals = [ ...
        lengthResidual;
        planeResidual];
end

function solutions = solveHipAnglesAnalytical( ...
    R_pelvis_ground, ...
    hipAdduction, ...
    pos_hip, ...
    pos_ankle, ...
    lfemur, ...
    ltibia, ...
    previousSolution)
%SOLVEHIPANGLESANALYTICAL
% Analytically solves:
%
%   norm(pos_ankle-pos_knee) = ltibia
%   dot(Z_femur,pos_ankle-pos_knee) = 0
%
% Rotation convention:
%
%   R_femur = R_pelvis * Rz(hipFlexion) ...
%                       * Rx(hipAdduction) ...
%                       * Ry(hipInternalRotation)
%
% Femur local longitudinal axis:
%
%   femurLocalAxis = [0;-1;0]
%
% Inputs:
%   previousSolution: optional [hipFlexion; hipInternalRotation] used
%                     to choose the continuous solution branch.
%
% Output:
%   solutions is a structure containing:
%       hipFlexion
%       hipInternalRotation
%       posKnee
%       ZFemur
%       allCandidates

    if nargin < 8
        previousSolution = [];
    end

    femurLocalAxis = [0; -1; 0];

    posHip = pos_hip(:);
    posAnkle = pos_ankle(:);

    % Transform the hip-to-ankle vector to the pelvis frame.
    hipToAnklePelvis = ...
        R_pelvis_ground.' * (posAnkle-posHip);

    % Femur direction after fixed hip adduction, but before flexion.
    Rx = rotationX(hipAdduction);
    femurDirection0 = Rx*femurLocalAxis;

    px = hipToAnklePelvis(1);
    py = hipToAnklePelvis(2);
    pz = hipToAnklePelvis(3);

    wx = femurDirection0(1);
    wy = femurDirection0(2);
    wz = femurDirection0(3);

    % From:
    %
    % ||p - lfemur*Rz(phi)*w||^2 = ltibia^2
    %
    % obtain:
    %
    % A*cos(phi) + B*sin(phi) = D.

    A = px*wx + py*wy;
    B = -px*wy + py*wx;

    C = (dot(hipToAnklePelvis,hipToAnklePelvis) ...
        + lfemur^2-ltibia^2)/(2*lfemur);

    D = C-pz*wz;

    rho = hypot(A,B);

    tolerance = 1e-12;

    if rho < tolerance
        error(['Hip-flexion solution is geometrically degenerate: ', ...
               'the flexion-dependent projection has near-zero magnitude.']);
    end

    cosineArgument = D/rho;

    % Allow small floating-point violations.
    if cosineArgument > 1+tolerance || cosineArgument < -1-tolerance
        error(['No exact hip-flexion solution exists. ', ...
               'The requested segment lengths are incompatible ', ...
               'with the hip and ankle positions. D/rho = %.8f.'], ...
               cosineArgument);
    end

    cosineArgument = max(-1,min(1,cosineArgument));

    phase = atan2(B,A);
    angleOffset = acos(cosineArgument);

    hipFlexionCandidates = [
        wrapToPiLocal(phase + angleOffset)
        wrapToPiLocal(phase - angleOffset)
    ];

    allCandidates = repmat(struct( ...
        'hipFlexion',[], ...
        'hipInternalRotation',[], ...
        'posKnee',[], ...
        'ZFemur',[], ...
        'lengthResidual',[], ...
        'planeResidual',[]),4,1);

    candidateNumber = 0;

    for flexionIndex = 1:2
        hipFlexion = hipFlexionCandidates(flexionIndex);

        RbeforeInternalRotation = ...
            R_pelvis_ground ...
            * rotationZ(hipFlexion) ...
            * Rx;

        % Internal rotation does not change the knee location because
        % femurLocalAxis is parallel to the local Y axis.
        posKnee = posHip ...
            + lfemur*RbeforeInternalRotation*femurLocalAxis;

        tibiaVector = posAnkle-posKnee;

        % Express tibia vector in the femur frame before Ry(psi).
        tibiaPreInternal = ...
            RbeforeInternalRotation.'*tibiaVector;

        tx = tibiaPreInternal(1);
        tz = tibiaPreInternal(3);

        % For the standard active Ry:
        %
        % Ry(psi)*e_z = [sin(psi); 0; cos(psi)]
        %
        % Plane constraint:
        %
        % tx*sin(psi) + tz*cos(psi) = 0.
        %
        % Therefore:
        %
        % psi = atan2(-tz,tx) + k*pi.

        if hypot(tx,tz) < tolerance
            % The plane residual is zero for every internal rotation.
            internalRotationCandidates = [0; pi];
        else
            internalRotationBase = atan2(-tz,tx);

            internalRotationCandidates = [
                wrapToPiLocal(internalRotationBase)
                wrapToPiLocal(internalRotationBase + pi)
            ];
        end

        for rotationIndex = 1:2
            hipInternalRotation = ...
                internalRotationCandidates(rotationIndex);

            femurRotation = ...
                RbeforeInternalRotation ...
                * rotationY(hipInternalRotation);

            ZFemur = femurRotation(:,3);
            ZFemur = ZFemur/norm(ZFemur);

            candidateNumber = candidateNumber+1;

            allCandidates(candidateNumber).hipFlexion = ...
                hipFlexion;

            allCandidates(candidateNumber).hipInternalRotation = ...
                hipInternalRotation;

            allCandidates(candidateNumber).posKnee = posKnee;
            allCandidates(candidateNumber).ZFemur = ZFemur;

            allCandidates(candidateNumber).lengthResidual = ...
                norm(tibiaVector)-ltibia;

            allCandidates(candidateNumber).planeResidual = ...
                dot(ZFemur,tibiaVector);
        end
    end

    % Choose the candidate closest to the previous frame.
    if ~isempty(previousSolution)
        previousSolution = previousSolution(:);

        candidateDistance = zeros(4,1);

        for candidateIndex = 1:4
            flexionDifference = angularDifference( ...
                allCandidates(candidateIndex).hipFlexion, ...
                previousSolution(1));

            rotationDifference = angularDifference( ...
                allCandidates(candidateIndex).hipInternalRotation, ...
                previousSolution(2));

            candidateDistance(candidateIndex) = ...
                hypot(flexionDifference,rotationDifference);
        end

        [~,selectedIndex] = min(candidateDistance);
    else
        % Default branch selection:
        % choose the solution with the smallest total angle magnitude.
        candidateMagnitude = zeros(4,1);

        for candidateIndex = 1:4
            candidateMagnitude(candidateIndex) = hypot( ...
                allCandidates(candidateIndex).hipFlexion, ...
                allCandidates(candidateIndex).hipInternalRotation);
        end

        [~,selectedIndex] = min(candidateMagnitude);
    end

    selected = allCandidates(selectedIndex);

    solutions.hipFlexion = selected.hipFlexion;
    solutions.hipInternalRotation = ...
        selected.hipInternalRotation;
    solutions.posKnee = selected.posKnee;
    solutions.ZFemur = selected.ZFemur;
    solutions.lengthResidual = selected.lengthResidual;
    solutions.planeResidual = selected.planeResidual;
    solutions.allCandidates = allCandidates;
end

function angle = wrapToPiLocal(angle)
    angle = atan2(sin(angle),cos(angle));
end


function difference = angularDifference(angle1,angle2)
    difference = atan2( ...
        sin(angle1-angle2), ...
        cos(angle1-angle2));
end
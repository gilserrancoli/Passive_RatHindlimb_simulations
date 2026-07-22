rat25_baseline=load('sol_val_optJointPassive_Datarat25baselineFB_tolem5_Jminres777_JpenKDT1e-3_JminKDT1e-5_JminInertiaP1e-2_nokneemarker2D_pert2_wKDinteractionsnoanklehip_ParamID_normCostF.mat');

%% By pert (groups of 4)
Databypert_baseline.pert30=load('solution_rat25_nonorm_pert30.mat');
Databypert_baseline.pert20=load('solution_rat25_nonorm_pert20.mat');
Databypert_baseline.pert10=load('solution_rat25_nonorm_pert10.mat');
Databypert_baseline.pert0=load('solution_rat25_nonorm_pert0.mat');
Databypert_baseline.pertm10=load('solution_rat25_nonorm_pertm10.mat');

%Data by trial
Databytrial_baseline.pert30{1}=load('solution_rat25_nonorm_0_pert30.mat');
Databytrial_baseline.pert30{2}=load('solution_rat25_nonorm_9_pert30.mat');
Databytrial_baseline.pert30{3}=load('solution_rat25_nonorm_0_pert30_back.mat');
Databytrial_baseline.pert30{4}=load('solution_rat25_nonorm_9_pert30_back.mat');

Databytrial_baseline.pert20{1}=load('solution_rat25_nonorm_1_pert20.mat');
Databytrial_baseline.pert20{2}=load('solution_rat25_nonorm_8_pert20.mat');
Databytrial_baseline.pert20{3}=load('solution_rat25_nonorm_1_pert20_back.mat');
Databytrial_baseline.pert20{4}=load('solution_rat25_nonorm_8_pert20_back.mat');

Databytrial_baseline.pert10{1}=load('solution_rat25_nonorm_2_pert10.mat');
Databytrial_baseline.pert10{2}=load('solution_rat25_nonorm_7_pert10.mat');
Databytrial_baseline.pert10{3}=load('solution_rat25_nonorm_2_pert10_back.mat');
Databytrial_baseline.pert10{4}=load('solution_rat25_nonorm_7_pert10_back.mat');

Databytrial_baseline.pert0{1}=load('solution_rat25_nonorm_3_pert0.mat');
Databytrial_baseline.pert0{2}=load('solution_rat25_nonorm_6_pert0.mat');
Databytrial_baseline.pert0{3}=load('solution_rat25_nonorm_3_pertm0_back.mat');
Databytrial_baseline.pert0{4}=load('solution_rat25_nonorm_6_pertm0_back.mat');

Databytrial_baseline.pertm10{1}=load('solution_rat25_nonorm_4_pertm10.mat');
Databytrial_baseline.pertm10{2}=load('solution_rat25_nonorm_5_pertm10.mat');
Databytrial_baseline.pertm10{3}=load('solution_rat25_nonorm_4_pertm10_back.mat');
Databytrial_baseline.pertm10{4}=load('solution_rat25_nonorm_5_pertm10_back.mat');

%Plot data considering same inertia per all perturbations
PlotR2allpert(rat25_baseline);

%Plot data considering same inertia per perturbation location
c='k';
PlotR2bypert(Databypert_baseline.pertm10,'m10',c);
PlotR2bypert(Databypert_baseline.pert0,'0',c);
PlotR2bypert(Databypert_baseline.pert10,'10',c);
PlotR2bypert(Databypert_baseline.pert20,'20',c);
PlotR2bypert(Databypert_baseline.pert30,'30',c);

%Plot data considering same inertia per perturbation location
PlotR2bytrial(Databytrial_baseline.pertm10,'m10');
PlotR2bytrial(Databytrial_baseline.pert0,'0');
PlotR2bytrial(Databytrial_baseline.pert10,'10');
PlotR2bytrial(Databytrial_baseline.pert20,'20');
PlotR2bytrial(Databytrial_baseline.pert30,'30');




function PlotR2bypert(datapert,namepert,c)
ntrials=length(datapert.sol_val.out);

for i=1:ntrials
    %hip
    mdl=fitlm(datapert.sol_val.out{i}(:,8),datapert.sol_val.out_model{i}(:,8));
    R2hip(i)=mdl.Rsquared.Ordinary;
    
    %knee
    mdl=fitlm(datapert.sol_val.out{i}(:,11),datapert.sol_val.out_model{i}(:,11));
    R2knee(i)=mdl.Rsquared.Ordinary;
    
    %ankle
    mdl=fitlm(datapert.sol_val.out{i}(:,12),datapert.sol_val.out_model{i}(:,12));
    R2ankle(i)=mdl.Rsquared.Ordinary;
end

switch namepert
    case 'm10'
        x=1;
    case '0'
        x=2;
    case '10'
        x=3;
    case '20'
        x=4;
    case '30'
        x=5;
end

figure(1)
xj=x+0.15*(rand(ntrials,1)-0.5);
subplot(1,3,1);
plot(xj,R2hip,'o','Color',c,'MarkerFaceColor',c);
hold all;
set(gca,'XTick',[1:5],'XTickLabel',{'-10','0','10','20','30'});
ylim([0 1]);
xlim([0 6]);
title('hip moment');
ylabel('R2 []');
xlabel('pert');

subplot(1,3,2);
plot(xj,R2knee,'o','Color',c,'MarkerFaceColor',c);
set(gca,'XTick',[1:5],'XTickLabel',{'-10','0','10','20','30'});
ylim([0 1]);
xlim([0 6]);
title('knee moment');
ylabel('R2 []');
xlabel('pert');
hold all;

subplot(1,3,3);
plot(xj,R2ankle,'o','Color',c,'MarkerFaceColor',c);
set(gca,'XTick',[1:5],'XTickLabel',{'-10','0','10','20','30'});
ylim([0 1]);
xlim([0 6]);
set(gca,'XTick',[1:5],'XTickLabel',{'-10','0','10','20','30'});
hold all;
title('ankle moment');
ylabel('R2 []');
xlabel('pert');

end


function PlotR2bytrial(datapert,namepert)
ntrials=length(datapert);

for i=1:ntrials
    %hip
    mdl=fitlm(datapert{i}.sol_val.out{1}(:,8),datapert{i}.sol_val.out_model{1}(:,8));
    R2hip(i)=mdl.Rsquared.Ordinary;
    
    %knee
    mdl=fitlm(datapert{i}.sol_val.out{1}(:,11),datapert{i}.sol_val.out_model{1}(:,11));
    R2knee(i)=mdl.Rsquared.Ordinary;
    
    %ankle
    mdl=fitlm(datapert{i}.sol_val.out{1}(:,12),datapert{i}.sol_val.out_model{1}(:,12));
    R2ankle(i)=mdl.Rsquared.Ordinary;
end

switch namepert
    case 'm10'
        x=1;
    case '0'
        x=2;
    case '10'
        x=3;
    case '20'
        x=4;
    case '30'
        x=5;
end

blue=[0 0.45 0.74];
figure(1)
xj=x+0.15*(rand(ntrials,1)-0.5);
subplot(1,3,1);
plot(xj,R2hip,'o','Color',blue,'MarkerFaceColor',blue);
hold all;
set(gca,'XTick',[1:5],'XTickLabel',{'-10','0','10','20','30'});
ylim([0 1]);
xlim([0 6]);

subplot(1,3,2);
plot(xj,R2knee,'o','Color',blue,'MarkerFaceColor',blue);
set(gca,'XTick',[1:5],'XTickLabel',{'-10','0','10','20','30'});
ylim([0 1]);
xlim([0 6]);
hold all;

subplot(1,3,3);
plot(xj,R2ankle,'o','Color',blue,'MarkerFaceColor',blue);
set(gca,'XTick',[1:5],'XTickLabel',{'-10','0','10','20','30'});
ylim([0 1]);
xlim([0 6]);
set(gca,'XTick',[1:5],'XTickLabel',{'-10','0','10','20','30'});
hold all;



end

function PlotR2allpert(alldata)

datapert.sol_val.out(1:4)=alldata.sol_val.out([5,6,15,16]);
datapert.sol_val.out_model(1:4)=alldata.sol_val.out_model([5,6,15,16]);
PlotR2bypert(datapert,'m10','r');
datapert.sol_val.out(1:4)=alldata.sol_val.out([4,7,14,17]);
datapert.sol_val.out_model(1:4)=alldata.sol_val.out_model([4,7,14,17]);
PlotR2bypert(datapert,'0','r');
datapert.sol_val.out(1:4)=alldata.sol_val.out([3,8,13,18]);
datapert.sol_val.out_model(1:4)=alldata.sol_val.out_model([3,8,13,18]);
PlotR2bypert(datapert,'10','r');
datapert.sol_val.out(1:4)=alldata.sol_val.out([2,9,12,19]);
datapert.sol_val.out_model(1:4)=alldata.sol_val.out_model([2,9,12,19]);
PlotR2bypert(datapert,'20','r');
datapert.sol_val.out(1:4)=alldata.sol_val.out([1,10,11,20]);
datapert.sol_val.out_model(1:4)=alldata.sol_val.out_model([1,10,11,20]);
PlotR2bypert(datapert,'30','r');

end
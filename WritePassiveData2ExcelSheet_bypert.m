clear all;
clc;

weighting_acc2var=1;

if weighting_acc2var
    % R21_b=load('sol_val_optJointPassive_Datarat21baselineFB_tolem5_Jminres777_JpenKDT1e-3_JminKDT1e-5_JminInertiaP1e-2_nokneemarker2D_pert2_wKDinteractionsnoanklehip_ParamID_normCostF.mat');
    % R21_a=load('sol_val_optJointPassive_Datarat21achillescutFB_tolem5_Jminres777_JpenKDT1e-3_JminKDT1e-5_JminInertiaP1e-2_nokneemarker2D_pert2_wKDinteractionsnoanklehip_ParamID_normCostF.mat');
    % R22_b=load('sol_val_optJointPassive_Datarat22baselineFB_tolem5_Jminres777_JpenKDT1e-3_JminKDT1e-5_JminInertiaP1e-2_nokneemarker2D_pert2_wKDinteractionsnoanklehip_ParamID_normCostF.mat');
    % R22_a=load('sol_val_optJointPassive_Datarat22achillescutFB_tolem5_Jminres777_JpenKDT1e-3_JminKDT1e-5_JminInertiaP1e-2_nokneemarker2D_pert2_wKDinteractionsnoanklehip_ParamID_normCostF.mat');
    % R23_b=load('sol_val_optJointPassive_Datarat23baselineFB_tolem5_Jminres777_JpenKDT1e-3_JminKDT1e-5_JminInertiaP1e-2_nokneemarker2D_pert2_wKDinteractionsnoanklehip_ParamID_normCostF.mat');
    % R23_a=load('sol_val_optJointPassive_Datarat23achillescutFB_tolem5_Jminres777_JpenKDT1e-3_JminKDT1e-5_JminInertiaP1e-2_nokneemarker2D_pert2_wKDinteractionsnoanklehip_ParamID_normCostF.mat');
    % R24_b=load('sol_val_optJointPassive_Datarat24baselineFB_tolem5_Jminres777_JpenKDT1e-3_JminKDT1e-5_JminInertiaP1e-2_nokneemarker2D_pert2_wKDinteractionsnoanklehip_ParamID_normCostF.mat');
    % R24_a=load('sol_val_optJointPassive_Datarat24achillescutFB_tolem5_Jminres777_JpenKDT1e-3_JminKDT1e-5_JminInertiaP1e-2_nokneemarker2D_pert2_wKDinteractionsnoanklehip_ParamID_normCostF.mat');
    % R25_b=load('sol_val_optJointPassive_Datarat25baselineFB_tolem5_Jminres777_JpenKDT1e-3_JminKDT1e-5_JminInertiaP1e-2_nokneemarker2D_pert2_wKDinteractionsnoanklehip_ParamID_normCostF.mat');
    % R25_a=load('sol_val_optJointPassive_Datarat25achillescutFB_tolem5_Jminres777_JpenKDT1e-3_JminKDT1e-5_JminInertiaP1e-2_nokneemarker2D_pert2_wKDinteractionsnoanklehip_ParamID_normCostF.mat');


    R21_b=ReadData('21','norm','baseline'); %load('sol_val_rat21_norm.mat');
    R21_a=ReadData('21','norm','achillescut'); %load('sol_val_rat21_achillescut_norm.mat');
    R22_b=ReadData('22','norm','baseline');
    R22_a=ReadData('22','norm','achillescut');
    R23_b=ReadData('23','norm','baseline');
    R23_a=ReadData('23','norm','achillescut');
    R24_b=ReadData('24','norm','baseline');
    R24_a=ReadData('24','norm','achillescut');
    R25_b=ReadData('25','norm','baseline');
    R25_a=ReadData('25','norm','achillescut');

    xlsxname='all_JointPassiveParam_FB_weightvar_byPert.xlsx';

else
    R21_b=load('sol_val_optJointPassive_Datarat21baselineFB_tolem5_Jminres777_JpenKDT1e-3_JminKDT1e-5_JminInertiaP1e-2_nokneemarker2D_pert2_wKDinteractionsnoanklehip_ParamID.mat');
    R21_a=load('sol_val_optJointPassive_Datarat21achillescutFB_tolem5_Jminres777_JpenKDT1e-3_JminKDT1e-5_JminInertiaP1e-2_nokneemarker2D_pert2_wKDinteractionsnoanklehip_ParamID.mat');
    R22_b=load('sol_val_optJointPassive_Datarat22baselineFB_tolem5_Jminres777_JpenKDT1e-3_JminKDT1e-5_JminInertiaP1e-2_nokneemarker2D_pert2_wKDinteractionsnoanklehip_ParamID.mat');
    R22_a=load('sol_val_optJointPassive_Datarat22achillescutFB_tolem5_Jminres777_JpenKDT1e-3_JminKDT1e-5_JminInertiaP1e-2_nokneemarker2D_pert2_wKDinteractionsnoanklehip_ParamID.mat');
    R23_b=load('sol_val_optJointPassive_Datarat23baselineFB_tolem5_Jminres777_JpenKDT1e-3_JminKDT1e-5_JminInertiaP1e-2_nokneemarker2D_pert2_wKDinteractionsnoanklehip_ParamID.mat');
    R23_a=load('sol_val_optJointPassive_Datarat23achillescutFB_tolem5_Jminres777_JpenKDT1e-3_JminKDT1e-5_JminInertiaP1e-2_nokneemarker2D_pert2_wKDinteractionsnoanklehip_ParamID.mat');
    R24_b=load('sol_val_optJointPassive_Datarat24baselineFB_tolem5_Jminres777_JpenKDT1e-3_JminKDT1e-5_JminInertiaP1e-2_nokneemarker2D_pert2_wKDinteractionsnoanklehip_ParamID.mat');
    R24_a=load('sol_val_optJointPassive_Datarat24achillescutFB_tolem5_Jminres777_JpenKDT1e-3_JminKDT1e-5_JminInertiaP1e-2_nokneemarker2D_pert2_wKDinteractionsnoanklehip_ParamID.mat');
    R25_b=load('sol_val_optJointPassive_Datarat25baselineFB_tolem5_Jminres777_JpenKDT1e-3_JminKDT1e-5_JminInertiaP1e-2_nokneemarker2D_pert2_wKDinteractionsnoanklehip_ParamID.mat');
    R25_a=load('sol_val_optJointPassive_Datarat25achillescutFB_tolem5_Jminres777_JpenKDT1e-3_JminKDT1e-5_JminInertiaP1e-2_nokneemarker2D_pert2_wKDinteractionsnoanklehip_ParamID.mat');

    xlsxname='all_JointPassiveParam_FB_NOweightvar_byPert.xlsx';

end
table_passive=cell2table({'Khip -10';'Kknee -10';'Kankle -10';'Khip-knee -10';'Kknee-ankle -10';...    
    'Khip 0';'Kknee 0';'Kankle 0';'Khip-knee 0';'Kknee-ankle 0';...
    'Khip 10';'Kknee 10';'Kankle 10';'Khip-knee 10';'Kknee-ankle 10';...
    'Khip 20';'Kknee 20';'Kankle 20';'Khip-knee 20';'Kknee-ankle 20';...
    'Khip 30';'Kknee 30';'Kankle 30';'Khip-knee 30';'Kknee-ankle 30';...
    'Dhip -10';'Dknee -10';'Dankle -10';'Dhip-knee -10';'Dknee-ankle -10';...
    'Dhip 0';'Dknee 0';'Dankle 0';'Dhip-knee 0';'Dknee-ankle 0';...
    'Dhip 10';'Dknee 10';'Dankle 10';'Dhip-knee 10';'Dknee-ankle 10';...
    'Dhip 20';'Dknee 20';'Dankle 20';'Dhip-knee 20';'Dknee-ankle 20';...
    'Dhip 30';'Dknee 30';'Dankle 30';'Dhip-knee 30';'Dknee-ankle 30';...
    'M0hip -10';'M0knee -10';'M0ankle -10';...
    'M0hip 0';'M0knee 0';'M0ankle 0'; ...
    'M0hip 10';'M0knee 10';'M0ankle 10'; ...
    'M0hip 20';'M0knee 20';'M0ankle 20'; ...
    'M0hip 30';'M0knee 30';'M0ankle 30'});

table_passive=addColumn(R21_b,'rat21 baseline',table_passive);
table_passive=addColumn(R21_a,'rat21 achillescut',table_passive);
table_passive=addColumn(R22_b,'rat22 baseline',table_passive);
table_passive=addColumn(R22_a,'rat22 achillescut',table_passive);
table_passive=addColumn(R23_b,'rat23 baseline',table_passive);
table_passive=addColumn(R23_a,'rat23 achillescut',table_passive);
table_passive=addColumn(R24_b,'rat24 baseline',table_passive);
table_passive=addColumn(R24_a,'rat24 achillescut',table_passive);
table_passive=addColumn(R25_b,'rat25 baseline',table_passive);
table_passive=addColumn(R25_a,'rat25 achillescut',table_passive);

writetable(table_passive,xlsxname,'Sheet','Stiffness_damping');

names = {'m femur'; 'I_G femur'; 'd CoM femur'; ...
         'm tibia'; 'I_G tibia'; 'd CoM tibia'; ...
         'm foot'; 'I_G foot'; 'd CoM foot'};

suffixes = [30 20 10 0 -10];

table_inertia = cell(numel(names)*numel(suffixes),1);

k = 1;
for s = suffixes
    for i = 1:numel(names)
        table_inertia{k} = sprintf('%s %d', names{i}, s);
        k = k + 1;
    end
end

table_inertia = cell2table(table_inertia,'VariableNames',{'Parameter'});table_inertia=addColumn_inertia(R21_b,'rat21 baseline',table_inertia);
table_inertia=addColumn_inertia(R21_a,'rat21 achillescut',table_inertia);
table_inertia=addColumn_inertia(R22_b,'rat22 baseline',table_inertia);
table_inertia=addColumn_inertia(R22_a,'rat22 achillescut',table_inertia);
table_inertia=addColumn_inertia(R23_b,'rat23 baseline',table_inertia);
table_inertia=addColumn_inertia(R23_a,'rat23 achillescut',table_inertia);
table_inertia=addColumn_inertia(R24_b,'rat24 baseline',table_inertia);
table_inertia=addColumn_inertia(R24_a,'rat24 achillescut',table_inertia);
table_inertia=addColumn_inertia(R25_b,'rat25 baseline',table_inertia);
table_inertia=addColumn_inertia(R25_a,'rat25 achillescut',table_inertia);

writetable(table_inertia,xlsxname,'Sheet','inertia parameters');


%% Write R2
[R21_b.r2_aux, R21_b.R2_aux, R21_b.R2demean_aux]=ComputeR2values(R21_b);
[R21_a.r2_aux, R21_a.R2_aux, R21_a.R2demean_aux]=ComputeR2values(R21_a);
[R22_b.r2, R22_b.R2, R22_b.R2demean]=ComputeR2values(R22_b);
[R22_a.r2, R22_a.R2, R22_a.R2demean]=ComputeR2values(R22_a);
[R23_b.r2, R23_b.R2, R23_b.R2demean]=ComputeR2values(R23_b);
[R23_a.r2, R23_a.R2, R23_a.R2demean]=ComputeR2values(R23_a);
[R24_b.r2, R24_b.R2, R24_b.R2demean]=ComputeR2values(R24_b);
[R24_a.r2, R24_a.R2, R24_a.R2demean]=ComputeR2values(R24_a);
[R25_b.r2, R25_b.R2, R25_b.R2demean]=ComputeR2values(R25_b);
[R25_a.r2, R25_a.R2, R25_a.R2demean]=ComputeR2values(R25_a);

%modifications for R21
%in r^2
R21_b.r2(1:4,:)=R21_b.r2_aux(1:4,:);
R21_b.r2(5:6,:)=NaN;
R21_b.r2(7:14,:)=R21_b.r2_aux(5:12,:);
R21_b.r2(15:16,:)=NaN;
R21_b.r2(17:20,:)=R21_b.r2_aux(13:16,:);

R21_a.r2(1:4,:)=R21_a.r2_aux(1:4,:);
R21_a.r2(5:6,:)=NaN;
R21_a.r2(7:14,:)=R21_a.r2_aux(5:12,:);
R21_a.r2(15:16,:)=NaN;
R21_a.r2(17:20,:)=R21_a.r2_aux(13:16,:);

%in R^2
R21_b.R2(1:4,:)=R21_b.R2_aux(1:4,:);
R21_b.R2(5:6,:)=NaN;
R21_b.R2(7:14,:)=R21_b.R2_aux(5:12,:);
R21_b.R2(15:16,:)=NaN;
R21_b.R2(17:20,:)=R21_b.R2_aux(13:16,:);

R21_a.R2(1:4,:)=R21_a.R2_aux(1:4,:);
R21_a.R2(5:6,:)=NaN;
R21_a.R2(7:14,:)=R21_a.R2_aux(5:12,:);
R21_a.R2(15:16,:)=NaN;
R21_a.R2(17:20,:)=R21_a.R2_aux(13:16,:);

%in R2_dmean
R21_b.R2demean(1:4,:)=R21_b.R2demean_aux(1:4,:);
R21_b.R2demean(5:6,:)=NaN;
R21_b.R2demean(7:14,:)=R21_b.R2demean_aux(5:12,:);
R21_b.R2demean(15:16,:)=NaN;
R21_b.R2demean(17:20,:)=R21_b.R2demean_aux(13:16,:);

R21_a.R2demean(1:4,:)=R21_a.R2demean_aux(1:4,:);
R21_a.R2demean(5:6,:)=NaN;
R21_a.R2demean(7:14,:)=R21_a.R2demean_aux(5:12,:);
R21_a.R2demean(15:16,:)=NaN;
R21_a.R2demean(17:20,:)=R21_a.R2demean_aux(13:16,:);

if weighting_acc2var
    xlsxname='all_R2values_FB_weightvar_byPert.xlsx';
else
    xlsxname='all_R2values_FB_NOweightvar_byPert.xlsx';
end
%r^2
table_r2=cell2table({'F- r^2 hip 30';'F- r^2 knee 30';'F- r^2 ankle 30';'F- r^2 multi 30';...
    'F- r^2 hip 20';'F- r^2 knee 20'; 'F- r^2 ankle 20';'F- r^2 multi 20';...
    'F- r^2 hip 10';'F- r^2 knee 10'; 'F- r^2 ankle 10';'F- r^2 multi 10';...
    'F- r^2 hip 0';'F- r^2 knee 0'; 'F- r^2 ankle 0';'F- r^2 multi 0';...
    'F- r^2 hip -10';'F- r^2 knee -10'; 'F- r^2 ankle -10';'F- r^2 multi -10';...
    'F- r^2 hip -10';'F- r^2 knee -10'; 'F- r^2 ankle -10';'F- r^2 multi -10';...
    'F- r^2 hip 0';'F- r^2 knee 0'; 'F- r^2 ankle 0';'F- r^2 multi 0';...
    'F- r^2 hip 10';'F- r^2 knee 10'; 'F- r^2 ankle 10';'F- r^2 multi 10';...
    'F- r^2 hip 20';'F- r^2 knee 20'; 'F- r^2 ankle 20';'F- r^2 multi 20';...
    'F- r^2 hip 30';'F- r^2 knee 30';'F- r^2 ankle 30';'F- r^2 multi 30';...

    'B - r^2 hip 30';'B - r^2 knee 30';'B - r^2 ankle 30';'B- r^2 multi 30';...
    'B - r^2 hip 20';'B - r^2 knee 20'; 'B - r^2 ankle 20';'B- r^2 multi 20';...
    'B - r^2 hip 10';'B - r^2 knee 10'; 'B - r^2 ankle 10';'B- r^2 multi 10';...
    'B - r^2 hip 0';'B - r^2 knee 0'; 'B - r^2 ankle 0';'B- r^2 multi 0';...
    'B - r^2 hip -10';'B - r^2 knee -10'; 'B - r^2 ankle -10';'B- r^2 multi -10';...
    'B - r^2 hip -10';'B - r^2 knee -10'; 'B - r^2 ankle -10';'B- r^2 multi -10';...
    'B - r^2 hip 0';'B - r^2 knee 0'; 'B - r^2 ankle 0';'B- r^2 multi 0';...
    'B - r^2 hip 10';'B - r^2 knee 10'; 'B - r^2 ankle 10';'B- r^2 multi 10';...
    'B - r^2 hip 20';'B - r^2 knee 20'; 'B - r^2 ankle 20';'B- r^2 multi 20';...
    'B - r^2 hip 30';'B - r^2 knee 30';'B - r^2 ankle 30';'B- r^2 multi 30';});
table_r2.Properties.VariableNames{1}='trials';
aux=R21_b.r2';
table_r2.('R21_base')(1:80)=aux(:);
aux=R21_a.r2';
table_r2.('R21_accut')(1:80)=aux(:);
aux=R22_b.r2';
table_r2.('R22_base')(1:80)=aux(:);
aux=R22_a.r2';
table_r2.('R22_accut')(1:80)=aux(:);
aux=R23_b.r2';
table_r2.('R23_base')(1:80)=aux(:);
aux=R23_a.r2';
table_r2.('R23_accut')(1:80)=aux(:);
aux=R24_b.r2';
table_r2.('R24_base')(1:80)=aux(:);
aux=R24_a.r2';
table_r2.('R24_accut')(1:80)=aux(:);
aux=R25_b.r2';
table_r2.('R25_base')(1:80)=aux(:);
aux=R25_a.r2';
table_r2.('R25_accut')(1:80)=aux(:);
writetable(table_r2,xlsxname,'Sheet','table_r2');

% R^2
table_R2=cell2table({'F- R^2 hip 30';'F- R^2 knee 30';'F- R^2 ankle 30';'F- R^2 multi 30';...
    'F- R^2 hip 20';'F- R^2 knee 20'; 'F- R^2 ankle 20';'F- R^2 multi 20';...
    'F- R^2 hip 10';'F- R^2 knee 10'; 'F- R^2 ankle 10';'F- R^2 multi 10';...
    'F- R^2 hip 0';'F- R^2 knee 0'; 'F- R^2 ankle 0';'F- R^2 multi 0';...
    'F- R^2 hip -10';'F- R^2 knee -10'; 'F- R^2 ankle -10';'F- R^2 multi -10';...
    'F- R^2 hip -10';'F- R^2 knee -10'; 'F- R^2 ankle -10';'F- R^2 multi -10';...
    'F- R^2 hip 0';'F- R^2 knee 0'; 'F- R^2 ankle 0';'F- R^2 multi 0';...
    'F- R^2 hip 10';'F- R^2 knee 10'; 'F- R^2 ankle 10';'F- R^2 multi 10';...
    'F- R^2 hip 20';'F- R^2 knee 20'; 'F- R^2 ankle 20';'F- R^2 multi 20';...
    'F- R^2 hip 30';'F- R^2 knee 30';'F- R^2 ankle 30';'F- R^2 multi 30';...

    'B - R^2 hip 30';'B - R^2 knee 30';'B - R^2 ankle 30';'B- R^2 multi 30';...
    'B - R^2 hip 20';'B - R^2 knee 20'; 'B - R^2 ankle 20';'B- R^2 multi 20';...
    'B - R^2 hip 10';'B - R^2 knee 10'; 'B - R^2 ankle 10';'B- R^2 multi 10';...
    'B - R^2 hip 0';'B - R^2 knee 0'; 'B - R^2 ankle 0';'B- R^2 multi 0';...
    'B - R^2 hip -10';'B - R^2 knee -10'; 'B - R^2 ankle -10';'B- R^2 multi -10';...
    'B - R^2 hip -10';'B - R^2 knee -10'; 'B - R^2 ankle -10';'B- R^2 multi -10';...
    'B - R^2 hip 0';'B - R^2 knee 0'; 'B - R^2 ankle 0';'B- R^2 multi 0';...
    'B - R^2 hip 10';'B - R^2 knee 10'; 'B - R^2 ankle 10';'B- R^2 multi 10';...
    'B - R^2 hip 20';'B - R^2 knee 20'; 'B - R^2 ankle 20';'B- R^2 multi 20';...
    'B - R^2 hip 30';'B - R^2 knee 30';'B - R^2 ankle 30';'B- R^2 multi 30';});
table_R2.Properties.VariableNames{1}='trials';
aux=R21_b.R2';
table_R2.('R21_base')(1:80)=aux(:);
aux=R21_a.R2';
table_R2.('R21_accut')(1:80)=aux(:);
aux=R22_b.R2';
table_R2.('R22_base')(1:80)=aux(:);
aux=R22_a.R2';
table_R2.('R22_accut')(1:80)=aux(:);
aux=R23_b.R2';
table_R2.('R23_base')(1:80)=aux(:);
aux=R23_a.R2';
table_R2.('R23_accut')(1:80)=aux(:);
aux=R24_b.R2';
table_R2.('R24_base')(1:80)=aux(:);
aux=R24_a.R2';
table_R2.('R24_accut')(1:80)=aux(:);
aux=R25_b.R2';
table_R2.('R25_base')(1:80)=aux(:);
aux=R25_a.R2';
table_R2.('R25_accut')(1:80)=aux(:);
writetable(table_R2,xlsxname,'Sheet','table_R_2');

% R^2 demean and normalized
table_zR2=cell2table({'F- zR^2 hip 30';'F- zR^2 knee 30';'F- zR^2 ankle 30';'F- zR^2 multi 30';...
    'F- zR^2 hip 20';'F- zR^2 knee 20'; 'F- zR^2 ankle 20';'F- zR^2 multi 20';...
    'F- zR^2 hip 10';'F- zR^2 knee 10'; 'F- zR^2 ankle 10';'F- zR^2 multi 10';...
    'F- zR^2 hip 0';'F- zR^2 knee 0'; 'F- zR^2 ankle 0';'F- zR^2 multi 0';...
    'F- zR^2 hip -10';'F- zR^2 knee -10'; 'F- zR^2 ankle -10';'F- zR^2 multi -10';...
    'F- zR^2 hip -10';'F- zR^2 knee -10'; 'F- zR^2 ankle -10';'F- zR^2 multi -10';...
    'F- zR^2 hip 0';'F- zR^2 knee 0'; 'F- zR^2 ankle 0';'F- zR^2 multi 0';...
    'F- zR^2 hip 10';'F- zR^2 knee 10'; 'F- zR^2 ankle 10';'F- zR^2 multi 10';...
    'F- zR^2 hip 20';'F- zR^2 knee 20'; 'F- zR^2 ankle 20';'F- zR^2 multi 20';...
    'F- zR^2 hip 30';'F- zR^2 knee 30';'F- zR^2 ankle 30';'F- zR^2 multi 30';...

    'B - zR^2 hip 30';'B - zR^2 knee 30';'B - zR^2 ankle 30';'B- zR^2 multi 30';...
    'B - zR^2 hip 20';'B - zR^2 knee 20'; 'B - zR^2 ankle 20';'B- zR^2 multi 20';...
    'B - zR^2 hip 10';'B - zR^2 knee 10'; 'B - zR^2 ankle 10';'B- zR^2 multi 10';...
    'B - zR^2 hip 0';'B - zR^2 knee 0'; 'B - zR^2 ankle 0';'B- zR^2 multi 0';...
    'B - zR^2 hip -10';'B - zR^2 knee -10'; 'B - zR^2 ankle -10';'B- zR^2 multi -10';...
    'B - zR^2 hip -10';'B - zR^2 knee -10'; 'B - zR^2 ankle -10';'B- zR^2 multi -10';...
    'B - zR^2 hip 0';'B - zR^2 knee 0'; 'B - zR^2 ankle 0';'B- zR^2 multi 0';...
    'B - zR^2 hip 10';'B - zR^2 knee 10'; 'B - zR^2 ankle 10';'B- zR^2 multi 10';...
    'B - zR^2 hip 20';'B - zR^2 knee 20'; 'B - zR^2 ankle 20';'B- zR^2 multi 20';...
    'B - zR^2 hip 30';'B - zR^2 knee 30';'B - zR^2 ankle 30';'B- zR^2 multi 30';});
table_zR2.Properties.VariableNames{1}='trials';
aux=R21_b.R2demean';
table_zR2.('R21_base')(1:80)=aux(:);
aux=R21_a.R2demean';
table_zR2.('R21_accut')(1:80)=aux(:);
aux=R22_b.R2demean';
table_zR2.('R22_base')(1:80)=aux(:);
aux=R22_a.R2demean';
table_zR2.('R22_accut')(1:80)=aux(:);
aux=R23_b.R2demean';
table_zR2.('R23_base')(1:80)=aux(:);
aux=R23_a.R2demean';
table_zR2.('R23_accut')(1:80)=aux(:);
aux=R24_b.R2demean';
table_zR2.('R24_base')(1:80)=aux(:);
aux=R24_a.R2demean';
table_zR2.('R24_accut')(1:80)=aux(:);
aux=R25_b.R2demean';
table_zR2.('R25_base')(1:80)=aux(:);
aux=R25_a.R2demean';
table_zR2.('R25_accut')(1:80)=aux(:);
writetable(table_zR2,xlsxname,'Sheet','table_R2_fromdemandata');


%%For ratXX - joint moments
time=R24_b.sol_val.tgrid{1}(R24_b.sol_val.t2plot{1});
sheetnames={'F_Pert_30_1','F_Pert_20_1','F_Pert_10_1','F_Pert_0_1','F_Pert_-10_1','F_Pert_-10_2','F_Pert_0_2','F_Pert_10_2','F_Pert_20_2','F_Pert_30_2',...
            'B_Pert_30_1','B_Pert_20_1','B_Pert_10_1','B_Pert_0_1','B_Pert_-10_1','B_Pert_-10_2','B_Pert_0_2','B_Pert_10_2','B_Pert_20_2','B_Pert_30_2'};
for i=1:size(R24_b.sol_val.out,2)
    data2table=[time' R24_b.sol_val.out{i}(:,[8 11 12]) R24_b.sol_val.out_model{i}(:,[8 11 12])];
    table_JM=array2table(data2table);
    table_JM.Properties.VariableNames={'time','Hip exp [Nm]','Knee exp [Nm]','Ankle exp [Nm]','Hip model [Nm]','Knee model [Nm]','Ankle model [Nm]'};
    writetable(table_JM,xlsxname,'Sheet',['base_' sheetnames{i}]);
end
for i=1:size(R24_b.sol_val.out,2)
    data2table=[time' R24_a.sol_val.out{i}(:,[8 11 12]) R24_a.sol_val.out_model{i}(:,[8 11 12])];
    table_JM=array2table(data2table);
    table_JM.Properties.VariableNames={'time','Hip exp [Nm]','Knee exp [Nm]','Ankle exp [Nm]','Hip model [Nm]','Knee model [Nm]','Ankle model [Nm]'};
    writetable(table_JM,xlsxname,'Sheet',['accut_' sheetnames{i}]);
end


function  ratdata=ReadData(rat,norm,type)

if strcmp(type,'baseline')
    sufixtype='';
else
    sufixtype='achillescut_';
end

for i=1:5
    switch i
        case 1
            sufixpert='pert30';
            idx=1:5;
            idx_K0=1:3;
            if strcmp(rat,'21')
                idx_trial=[1 8 9 16];
            else
                idx_trial=[1 10 11 20];
            end
        case 2
            sufixpert='pert20';
            idx=6:10;
            idx_K0=4:6;
            if strcmp(rat,'21')
                idx_trial=[2 7 10 15];
            else
                idx_trial=[2 9 12 19];
            end
        case 3
            sufixpert='pert10';
            idx=11:15;
            idx_K0=7:9;
            if strcmp(rat,'21')
                idx_trial=[3 6 11 14];
            else
                idx_trial=[3 8 13 18];
            end
        case 4
            sufixpert='pert0';
            idx=16:20;
            idx_K0=10:12;
            if strcmp(rat,'21')
                idx_trial=[4 5 12 13];
            else
                idx_trial=[4 7 14 17];
            end
        case 5
            sufixpert='pertm10';
            idx=21:25;
            idx_K0=13:15;
            idx_trial=[5 6 15 16];
    end
    if strcmp(rat,'21') && strcmp(sufixpert,'pertm10')
    else
        data=load(['sol_val_rat' rat '_' sufixtype 'uniquepass_' norm '_' sufixpert '.mat']);

        
            
        ratdata.sol_val.Kstiff_unsc(idx,1)=data.sol_val.Kstiff_unsc;
        ratdata.sol_val.K0_unsc(idx_K0,1)=data.sol_val.K0_unsc;
        ratdata.sol_val.Ddamp_unsc(idx)=data.sol_val.Ddamp_unsc;
        ratdata.sol_val.inertiaParam_unsc(:,i)=data.sol_val.inertiaParam_unsc;
        
        ratdata.sol_val.out(idx_trial)=data.sol_val.out;
        ratdata.sol_val.out_model(idx_trial)=data.sol_val.out_model;

        ratdata.sol_val.tgrid(idx_trial)=data.sol_val.tgrid;
        ratdata.sol_val.t2plot(idx_trial)=data.sol_val.t2plot;
    end
    

end



end

function table_passive=addColumn(data,columnLabel,table_passive)
ncolumns_i=size(table_passive,2);
nrows=size(table_passive,1);
n=ncolumns_i+1;
table_passive.(columnLabel)=NaN(nrows,1);
if contains(columnLabel,'21')
    table_passive.(columnLabel)(1:20)=data.sol_val.Kstiff_unsc;
    table_passive.(columnLabel)(26:45)=data.sol_val.Ddamp_unsc;
    if isfield(data.sol_val,'theta0_unsc')
        table_passive.(columnLabel)(51:62)=data.sol_val.Kstiff_unsc([1:3 6:8 11:13 16:18]).*data.sol_val.theta0_unsc;
    else
        table_passive.(columnLabel)(51:62)=data.sol_val.K0_unsc;
    end
else
    table_passive.(columnLabel)(1:25)=data.sol_val.Kstiff_unsc;
    table_passive.(columnLabel)(26:50)=data.sol_val.Ddamp_unsc;
    if isfield(data.sol_val,'theta0_unsc')
        table_passive.(columnLabel)(51:65)=data.sol_val.Kstiff_unsc([1:3 6:8 11:13 16:18 20:22]).*data.sol_val.theta0_unsc;
    else
        table_passive.(columnLabel)(51:65)=data.sol_val.K0_unsc;
    end
end

end

function table_inertia=addColumn_inertia(data,columnLabel,table_inertia)
ncolumns_i=size(table_inertia,2);
nrows=size(table_inertia,1);
n=ncolumns_i+1;
table_inertia.(columnLabel)=NaN(nrows,1);
    for i=1:5
        if contains(columnLabel,'rat21')&&i==5
        else
            try
            table_inertia.(columnLabel)((i-1)*9+1:i*9)=data.sol_val.inertiaParam_unsc(:,i);
            catch
                keyboard;
            end
        end
    end

end

function [r2all, R2all, R2all_dmean]=ComputeR2values(data)

    for i=1:size(data.sol_val.out_model,2)
        I=[8,11,12]; %hip knee ankle
        for j=1:3 
            try
            r_aux=corrcoef(data.sol_val.out_model{i}(:,I(j)),data.sol_val.out{i}(:,I(j)));
            catch
                keyboard;
            end
            r2all(i,j)=r_aux(2)^2;
        end
        v1=data.sol_val.out_model{i}(:,[8 11 12]);
        v2=data.sol_val.out{i}(:,[8 11 12]);
        r_aux=corrcoef(v1(:),v2(:));
        % r_aux=corrcoef(data.sol_val.out_model{i}(:,[8 11 12]),data.sol_val.out{i}(:,[8 11 12]))
        r2all(i,4)=r_aux(2)^2;
    end
    

    for i=1:size(data.sol_val.out_model,2)
        I=[8,11,12]; %hip knee ankle
        for j=1:3
            y=data.sol_val.out{i}(:,I(j));
            yhat=data.sol_val.out_model{i}(:,I(j));
            R2all(i,j) = 1 - sum((y - yhat).^2) / sum((y - mean(y)).^2);
        end
        v1=data.sol_val.out_model{i}(:,[8 11 12]);
        v2=data.sol_val.out{i}(:,[8 11 12]);
        y=v2(:);
        yhat=v1(:);
        R2all(i,4)= 1 - sum((y - yhat).^2) / sum((y - mean(y)).^2);
    end
    
    %% Demean vectors and divide by standard deviation
    for i=1:size(data.sol_val.out_model,2)
        I=[8,11,12]; %hip knee ankle
        for j=1:3 
            % mean_model=mean(data.sol_val.out_model{i}(:,I(j)));
            % std_model=std(data.sol_val.out_model{i}(:,I(j)));
            mean_out=mean(data.sol_val.out{i}(:,I(j)));
            std_out=std(data.sol_val.out{i}(:,I(j)));
            r_aux=corrcoef((data.sol_val.out_model{i}(:,I(j))-mean_out)/std_out,(data.sol_val.out{i}(:,I(j))-mean_out)/std_out);
            r2all_dmean(i,j)=r_aux(2)^2;
        end
        % std_model=std(data.sol_val.out_model{i}(:,[8 11 12]));
        % mean_model=mean(data.sol_val.out_model{i}(:,[8 11 12]));
        std_out=std(data.sol_val.out{i}(:,[8 11 12]));
        mean_out=mean(data.sol_val.out{i}(:,[8 11 12]));
        v1=(data.sol_val.out_model{i}(:,[8 11 12])-mean_out)./std_out;
        v2=(data.sol_val.out{i}(:,[8 11 12])-mean_out)./std_out;
        r_aux=corrcoef(v1(:),v2(:));
        % r_aux=corrcoef(data.sol_val.out_model{i}(:,[8 11 12]),data.sol_val.out{i}(:,[8 11 12]))
        r2all_dmean(i,4)=r_aux(2)^2;
    end
    

    for i=1:size(data.sol_val.out_model,2)
        I=[8,11,12]; %hip knee ankle
        for j=1:3
            % y=data.sol_val.out{i}(:,I(j))-mean(data.sol_val.out{i}(:,I(j)));
            % yhat=data.sol_val.out_model{i}(:,I(j))-mean(data.sol_val.out_model{i}(:,I(j)));
            % mean_model=mean(data.sol_val.out_model{i}(:,I(j)));
            % std_model=std(data.sol_val.out_model{i}(:,I(j)));
            mean_out=mean(data.sol_val.out{i}(:,I(j)));
            std_out=std(data.sol_val.out{i}(:,I(j)));

            y=(data.sol_val.out{i}(:,I(j))-mean_out)./std_out
            yhat=(data.sol_val.out_model{i}(:,I(j))-mean_out)./std_out;
            R2all_dmean(i,j) = 1 - sum((y - yhat).^2) / sum((y - mean(y)).^2);
        end
        % std_model=std(data.sol_val.out_model{i}(:,[8 11 12]));
        % mean_model=mean(data.sol_val.out_model{i}(:,[8 11 12]));
        std_out=std(data.sol_val.out{i}(:,[8 11 12]));
        mean_out=mean(data.sol_val.out{i}(:,[8 11 12]));
        v1=(data.sol_val.out_model{i}(:,[8 11 12])-mean_out)./std_out;
        v2=(data.sol_val.out{i}(:,[8 11 12])-mean_out)./std_out;
        y=v2(:);
        yhat=v1(:);
        R2all_dmean(i,4)= 1 - sum((y - yhat).^2) / sum((y - mean(y)).^2);
    end

end
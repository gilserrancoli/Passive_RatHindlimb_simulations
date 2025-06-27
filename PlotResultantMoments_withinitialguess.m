S1=load('sol_val_optKinematics_optJointPassive_optInertia_DataMay2025IntactFB_tolem5_Jqtrack100_Jqtrackdot10_Jtrackqd2dot001_JpenKDT1e-3_JminKDT1e-5_JminInertiaP1e-2_rel2Dangles_N40.mat');
S2=load('sol_val_optKinematics_optJointPassive_optInertia_DataMay2025achillescutFB_tolem5_Jqtrack100_Jqtrackdot10_Jtrackqd2dot001_JpenKDT1e-3_JminKDT1e-5_JminInertiaP1e-2_rel2Dangles_N40.mat');

S1.sol_val.out_exp=CalculateMoments_initialguess(S1.sol_val,F);
S2.sol_val.out_exp=CalculateMoments_initialguess(S2.sol_val,F);

blue=[0 0.45 0.74];
red=[0.85 0.33 0.1];

for i=1:length(S1.sol_val.out_opt)
    figure(i)
    set(gcf,'Position',[360,97.7,749.7,216]);
    subplot(1,3,1)
    plot(S1.sol_val.tgrid_col{i},S1.sol_val.out_opt{i}(:,8),'Color',blue);
    hold all;
    plot(S1.sol_val.tgrid_col{i},S1.sol_val.out_exp{i}(:,8),'--','Color',blue);
    plot(S2.sol_val.tgrid_col{i},S2.sol_val.out_opt{i}(:,8),'Color',red);
    plot(S2.sol_val.tgrid_col{i},S2.sol_val.out_exp{i}(:,8),'--','Color',red);
    title('hip flexion');
    ylabel('M [Nm]')
    xlabel('time [s]')

    subplot(1,3,2)
    plot(S1.sol_val.tgrid_col{i},S1.sol_val.out_opt{i}(:,11),'Color',blue);
    hold all;
    plot(S1.sol_val.tgrid_col{i},S1.sol_val.out_exp{i}(:,11),'--','Color',blue);
    plot(S2.sol_val.tgrid_col{i},S2.sol_val.out_opt{i}(:,11),'Color',red);
    plot(S2.sol_val.tgrid_col{i},S2.sol_val.out_exp{i}(:,11),'--','Color',red);
    title('knee flexion')
    ylabel('M [Nm]')
    xlabel('time [s]')

    subplot(1,3,3)
    plot(S1.sol_val.tgrid_col{i},S1.sol_val.out_opt{i}(:,12),'Color',blue);
    hold all;
    plot(S1.sol_val.tgrid_col{i},S1.sol_val.out_exp{i}(:,12),'--','Color',blue);
    plot(S2.sol_val.tgrid_col{i},S2.sol_val.out_opt{i}(:,12),'Color',red);
    plot(S2.sol_val.tgrid_col{i},S2.sol_val.out_exp{i}(:,12),'--','Color',red);

    title('ankle flexion')
    ylabel('M [Nm]')
    xlabel('time [s]')
    legend({'intact','intact IG','achilles cut','achilles cut IG'},'box','off','Position',[0.843,0.1591,0.14,0.15])
end


function out_exp=CalculateMoments_initialguess(sol_val,F)
    import casadi.*
    scaling=sol_val.scaling;
    Options=sol_val.Options;
    forces_prescribed=sol_val.force_prescribed;

    for i=1:length(sol_val.QsQdot_prescribed)
        if ~isempty(sol_val.QsQdot_prescribed{i})
            all_QsQdot_exp{i}(:,1:14)=sol_val.QsQdot_prescribed{i}(sol_val.t2plot{i},:);
            all_QsQdot_exp{i}(:,15:28)=sol_val.guess.QsQdots{i}(sol_val.t2plot{i},:).*scaling.QsQdots; 
    
            all_Qd2dot_exp{i}(:,1:7)=sol_val.Qd2dot_prescribed{i}(sol_val.t2plot{i},:);
            all_Qd2dot_exp{i}(:,8:14)=sol_val.guess.Qd2dots{i}(sol_val.t2plot{i},:).*scaling.qd2dot;
    
            for k=1:size(all_Qd2dot_exp{i},1)
                if Options.optInertiaParam
                    out_exp{i}(k,:)=full(F([all_QsQdot_exp{i}(k,:)';all_Qd2dot_exp{i}(k,:)';forces_prescribed{i}(k,:)';sol_val.inertiaParam_unsc]));
                else
                    out_exp{i}(k,:)=full(F([all_QsQdot_exp{i}(k,:)';all_Qd2dot_exp{i}(k,:)';forces_prescribed{i}(k,:)']));
                end
            end
        end
    end

end

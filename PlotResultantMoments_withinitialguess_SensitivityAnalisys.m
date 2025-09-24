S1=load('sol_val_optJointPassive_optInertia_DataMay2025IntactFB_tolem5_Jqtrack100_Jqtrackdot10_JpenKDT1e-3_JminKDT1e-5_JminInertiaP1e-2_nokneemarker2D_N40.mat');
S2=load('sol_val_optJointPassive_optInertia_DataMay2025achillescutFB_tolem5_Jqtrack100_Jqtrackdot10_JpenKDT1e-3_JminKDT1e-5_JminInertiaP1e-2_nokneemarker2D_N40.mat');

[S1.sol_val.out_exp, S1.sol_val.out_exp_zeroqqdotforce, S1.sol_val.out_exp_onlyqdot, S1.sol_val.out_exp_onlyqd2dot, S1.sol_val.out_exp_onlyforces]=CalculateMoments_initialguess(S1.sol_val,F);
[S1.sol_val.out, S1.sol_val.out_zeroqqdotforce, S1.sol_val.out_onlyqdot, S1.sol_val.out_onlyqd2dot, S1.sol_val.out_onlyforces]=CalculateMoments(S1.sol_val,F);

% S2.sol_val.out_exp=CalculateMoments_initialguess(S2.sol_val,F);

blue=[0 0.45 0.74];
red=[0.85 0.33 0.1];
green=[0, 0.66,0.33];

for i=1:length(S1.sol_val.out_opt)
    figure(i)
    set(gcf,'Position',[360,97.7,564.3,542.6]);
    subplot(3,2,1)
    plot(S1.sol_val.tgrid_col{i},S1.sol_val.out_opt{i}(:,8),'LineWidth',2,'Color','k');
    hold all;
    plot(S1.sol_val.tgrid_col{i},S1.sol_val.out_zeroqqdotforce{i}(:,8),'Color',green);
    plot(S1.sol_val.tgrid_col{i},S1.sol_val.out_onlyforces{i}(:,8)-S1.sol_val.out_zeroqqdotforce{i}(:,8),'--','Color','r');
    plot(S1.sol_val.tgrid_col{i},S1.sol_val.out_onlyqdot{i}(:,8)-S1.sol_val.out_zeroqqdotforce{i}(:,8),'Color','b');
    plot(S1.sol_val.tgrid_col{i},S1.sol_val.out_onlyqd2dot{i}(:,8)-S1.sol_val.out_zeroqqdotforce{i}(:,8),'--','Color','b');
    title('hip flexion kin. optimized');

    subplot(3,2,2)
    plot(S1.sol_val.tgrid_col{i},S1.sol_val.out_exp{i}(:,8),'LineWidth',2,'Color','k');
    hold all;
    plot(S1.sol_val.tgrid_col{i},S1.sol_val.out_exp_zeroqqdotforce{i}(:,8),'Color',green);
    plot(S1.sol_val.tgrid_col{i},S1.sol_val.out_exp_onlyforces{i}(:,8)-S1.sol_val.out_exp_zeroqqdotforce{i}(:,8),'--','Color','r');
    plot(S1.sol_val.tgrid_col{i},S1.sol_val.out_exp_onlyqdot{i}(:,8)-S1.sol_val.out_exp_zeroqqdotforce{i}(:,8),'Color','b');
    plot(S1.sol_val.tgrid_col{i},S1.sol_val.out_exp_onlyqd2dot{i}(:,8)-S1.sol_val.out_exp_zeroqqdotforce{i}(:,8),'--','Color','b');
    % plot(S2.sol_val.tgrid_col{i},S2.sol_val.out_opt{i}(:,8),'Color',red);
    % plot(S2.sol_val.tgrid_col{i},S2.sol_val.out_exp{i}(:,8),'--','Color',red);
    title('hip flexion kin. initial guess');
    ylabel('M [Nm]')
    xlabel('time [s]')

    subplot(3,2,3)
    plot(S1.sol_val.tgrid_col{i},S1.sol_val.out_opt{i}(:,11),'LineWidth',2,'Color','k');
    hold all;
    plot(S1.sol_val.tgrid_col{i},S1.sol_val.out_zeroqqdotforce{i}(:,11),'Color',green);
    plot(S1.sol_val.tgrid_col{i},S1.sol_val.out_onlyforces{i}(:,11)-S1.sol_val.out_zeroqqdotforce{i}(:,11),'--','Color','r');
    plot(S1.sol_val.tgrid_col{i},S1.sol_val.out_onlyqdot{i}(:,11)-S1.sol_val.out_zeroqqdotforce{i}(:,11),'Color','b');
    plot(S1.sol_val.tgrid_col{i},S1.sol_val.out_onlyqd2dot{i}(:,11)-S1.sol_val.out_zeroqqdotforce{i}(:,11),'--','Color','b');
    title('knee flexion kin. optimized');

    subplot(3,2,4)
    plot(S1.sol_val.tgrid_col{i},S1.sol_val.out_exp{i}(:,11),'LineWidth',2,'Color','k');
    hold all;
    plot(S1.sol_val.tgrid_col{i},S1.sol_val.out_exp_zeroqqdotforce{i}(:,11),'Color',green);
    plot(S1.sol_val.tgrid_col{i},S1.sol_val.out_exp_onlyforces{i}(:,11)-S1.sol_val.out_exp_zeroqqdotforce{i}(:,11),'--','Color','r');
    plot(S1.sol_val.tgrid_col{i},S1.sol_val.out_exp_onlyqdot{i}(:,11)-S1.sol_val.out_exp_zeroqqdotforce{i}(:,11),'Color','b');
    plot(S1.sol_val.tgrid_col{i},S1.sol_val.out_exp_onlyqd2dot{i}(:,11)-S1.sol_val.out_exp_zeroqqdotforce{i}(:,11),'--','Color','b');
    % plot(S2.sol_val.tgrid_col{i},S2.sol_val.out_opt{i}(:,11),'Color',red);
    % plot(S2.sol_val.tgrid_col{i},S2.sol_val.out_exp{i}(:,11),'--','Color',red);
    title('knee flexion kin. initial guess')
    ylabel('M [Nm]')
    xlabel('time [s]')

    subplot(3,2,5)
    plot(S1.sol_val.tgrid_col{i},S1.sol_val.out_opt{i}(:,12),'LineWidth',2,'Color','k');
    hold all;
    plot(S1.sol_val.tgrid_col{i},S1.sol_val.out_zeroqqdotforce{i}(:,12),'Color',green);
    plot(S1.sol_val.tgrid_col{i},S1.sol_val.out_onlyforces{i}(:,12)-S1.sol_val.out_zeroqqdotforce{i}(:,12),'--','Color','r');
    plot(S1.sol_val.tgrid_col{i},S1.sol_val.out_onlyqdot{i}(:,12)-S1.sol_val.out_zeroqqdotforce{i}(:,12),'Color','b');
    plot(S1.sol_val.tgrid_col{i},S1.sol_val.out_onlyqd2dot{i}(:,12)-S1.sol_val.out_zeroqqdotforce{i}(:,12),'--','Color','b');
    title('ankle flexion kin. optimized')

    subplot(3,2,6)
    plot(S1.sol_val.tgrid_col{i},S1.sol_val.out_exp{i}(:,12),'LineWidth',2,'Color','k');
    hold all;
    plot(S1.sol_val.tgrid_col{i},S1.sol_val.out_exp_zeroqqdotforce{i}(:,12),'Color',green);
    plot(S1.sol_val.tgrid_col{i},S1.sol_val.out_exp_onlyforces{i}(:,12)-S1.sol_val.out_exp_zeroqqdotforce{i}(:,12),'--','Color','r');
    plot(S1.sol_val.tgrid_col{i},S1.sol_val.out_exp_onlyqdot{i}(:,12)-S1.sol_val.out_exp_zeroqqdotforce{i}(:,12),'Color','b');
    plot(S1.sol_val.tgrid_col{i},S1.sol_val.out_exp_onlyqd2dot{i}(:,12)-S1.sol_val.out_exp_zeroqqdotforce{i}(:,12),'--','Color','b');
    % plot(S2.sol_val.tgrid_col{i},S2.sol_val.out_opt{i}(:,12),'Color',red);
    % plot(S2.sol_val.tgrid_col{i},S2.sol_val.out_exp{i}(:,12),'--','Color',red);
    title('ankle flexion kin. initial guess')
    ylabel('M [Nm]')
    xlabel('time [s]')
legend({'$total$', '$G(q)$', '$\tau_{GRF}$', '$C(q,\dot{q})$', '$M\ddot{q}$'}, ...
       'Interpreter', 'latex', 'Box', 'off', 'Position', [0.843, 0.1591, 0.14, 0.15]);
end


function [out_exp, out_exp_zeroqqdotforce, out_exp_onlyqdot, out_exp_onlyqd2dot, out_exp_onlyforces]=CalculateMoments_initialguess(sol_val,F)
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
                    qsqdot=all_QsQdot_exp{i}(k,:)';
                    qsqdot(2:2:end)=0;
                    out_exp_zeroqqdotforce{i}(k,:)=full(F([qsqdot;zeros(14,1);zeros(9,1);sol_val.inertiaParam_unsc]));
                    out_exp_onlyqdot{i}(k,:)=full(F([all_QsQdot_exp{i}(k,:)';zeros(14,1);zeros(9,1);sol_val.inertiaParam_unsc]));
                    out_exp_onlyqd2dot{i}(k,:)=full(F([qsqdot;all_Qd2dot_exp{i}(k,:)';zeros(9,1);sol_val.inertiaParam_unsc]));
                    out_exp_onlyforces{i}(k,:)=full(F([qsqdot;zeros(14,1);forces_prescribed{i}(k,:)';sol_val.inertiaParam_unsc]));
                else
                    out_exp{i}(k,:)=full(F([all_QsQdot_exp{i}(k,:)';all_Qd2dot_exp{i}(k,:)';forces_prescribed{i}(k,:)']));
                    out_exp_onlyqdot=[];
                    out_exp_onlyqd2dot=[];
                end
            end
        end
    end

end
function [out, out_zeroqqdotforce, out_onlyqdot, out_onlyqd2dot, out_onlyforces]=CalculateMoments(sol_val,F)
import casadi.*
    scaling=sol_val.scaling;
    Options=sol_val.Options;
    forces_prescribed=sol_val.force_prescribed;

    for i=1:length(sol_val.QsQdot_prescribed)
        if ~isempty(sol_val.QsQdot_prescribed{i})
            all_QsQdot{i}(:,1:14)=sol_val.QsQdot_prescribed{i}(sol_val.t2plot{i},:);
            all_QsQdot{i}(:,15:28)=sol_val.QsQdots_col_unsc{i}'; 
    
            all_Qd2dot{i}(:,1:7)=sol_val.Qd2dot_prescribed{i}(sol_val.t2plot{i},:);
            all_Qd2dot{i}(:,8:14)=sol_val.Qd2dot_col_unsc{i}';
    
            for k=1:size(all_Qd2dot{i},1)
                if Options.optInertiaParam
                    out{i}(k,:)=full(F([all_QsQdot{i}(k,:)';all_Qd2dot{i}(k,:)';forces_prescribed{i}(k,:)';sol_val.inertiaParam_unsc]));
                    qsqdot=all_QsQdot{i}(k,:)';
                    qsqdot(2:2:end)=0;
                    out_zeroqqdotforce{i}(k,:)=full(F([qsqdot;zeros(14,1);zeros(9,1);sol_val.inertiaParam_unsc]));
                    out_onlyqdot{i}(k,:)=full(F([all_QsQdot{i}(k,:)';zeros(14,1);zeros(9,1);sol_val.inertiaParam_unsc]));
                    out_onlyqd2dot{i}(k,:)=full(F([qsqdot;all_Qd2dot{i}(k,:)';zeros(9,1);sol_val.inertiaParam_unsc]));
                    out_onlyforces{i}(k,:)=full(F([qsqdot;zeros(14,1);forces_prescribed{i}(k,:)';sol_val.inertiaParam_unsc]));
                else
                    out{i}(k,:)=full(F([all_QsQdotp{i}(k,:)';all_Qd2dot{i}(k,:)';forces_prescribed{i}(k,:)']));
                    out_onlyqdot=[];
                    out_onlyqd2dot=[];
                    out_onlyforces=[];
                end
            end
        end
    end

end
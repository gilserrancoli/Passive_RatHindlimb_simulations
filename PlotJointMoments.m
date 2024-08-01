function PlotJointMoments(sol_val,i,Options)

for j=8:size(sol_val.out_opt{i},2)
    subplot(4,2,j-7)
    plot(sol_val.tgrid_col{i}, sol_val.out_opt{i}(:,j),'LineWidth',2); % total moment
    hold all;
    if Options.optimizeMuscleProp
        switch j %muscle moment
            case 8 
                plot(sol_val.tgrid_col{i},sum(sol_val.T_hip_flx_opt,2));
            case 9
                plot(sol_val.tgrid_col{i},sum(sol_val.T_hip_add_opt,2));
            case 10
                plot(sol_val.tgrid_col{i},sum(sol_val.T_hip_int_opt,2));
            case 11
                plot(sol_val.tgrid_col{i},sum(sol_val.T_knee_flx_opt,2));
            case 12
                plot(sol_val.tgrid_col{i},sum(sol_val.T_ankle_flx_opt,2));
            case 13
                plot(sol_val.tgrid_col{i},sum(sol_val.T_ankle_add_opt,2));
            case 14
                plot(sol_val.tgrid_col{i},sum(sol_val.T_ankle_int_opt,2));
        end
    end
    if Options.optimizePassiveJointEl
        switch j %Joint passive moment
            case 8
                plot(sol_val.tgrid_col{i},sol_val.PassiveM_hip_flx_opt);
            case 9
                plot(sol_val.tgrid_col{i},sol_val.PassiveM_hip_add_opt);
            case 10
                plot(sol_val.tgrid_col{i},sol_val.PassiveM_hip_int_opt);
            case 11
                plot(sol_val.tgrid_col{i},sol_val.PassiveM_knee_flx_opt);
            case 12
                plot(sol_val.tgrid_col{i},sol_val.PassiveM_ankle_flx_opt);
            case 13
                plot(sol_val.tgrid_col{i},sol_val.PassiveM_ankle_add_opt);
            case 14
                plot(sol_val.tgrid_col{i},sol_val.PassiveM_ankle_int_opt);
        end
    end
    plot(sol_val.tgrid_col{i},sol_val.res_col_unsc{i}(j-7,:)');
    title(sol_val.name_dofs{j-7});
end
if Options.optimizeMuscleProp&~Options.optimizePassiveJointEl
    legend('total moment','sum muscle moments','residual');
elseif  ~Options.optimizeMuscleProp&Options.optimizePassiveJointEl
    legend('total moment','sum joint passive moments','residual');
elseif Options.optimizeMuscleProp&Options.optimizePassiveJointEl
    legend('total moment','sum muscle moments','sum joint passive moments','residual');
end

    
    
function PlotJointMoments_SagittalPlane(sol_val,i,Options)

ii=1;
for j=[8 11 12]
    subplot(3,1,ii)
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
                plot(sol_val.tgrid_col{i},sol_val.PassiveM_hip_flx_opt{i});
            case 9
                plot(sol_val.tgrid_col{i},sol_val.PassiveM_hip_add_opt{i});
            case 10
                plot(sol_val.tgrid_col{i},sol_val.PassiveM_hip_int_opt{i});
            case 11
                plot(sol_val.tgrid_col{i},sol_val.PassiveM_knee_flx_opt{i});
            case 12
                plot(sol_val.tgrid_col{i},sol_val.PassiveM_ankle_flx_opt{i});
            case 13
                plot(sol_val.tgrid_col{i},sol_val.PassiveM_ankle_add_opt{i});
            case 14
                plot(sol_val.tgrid_col{i},sol_val.PassiveM_ankle_int_opt{i});
        end
    end
    if isfield(sol_val,'res_col_unsc')
        plot(sol_val.tgrid_col{i},sol_val.res_col_unsc{i}(j-7,:)');
    end
    title(sol_val.name_dofs{j-7});
    ylabel('Moment [Nm]');
    ii=ii+1;
end
if Options.optimizeMuscleProp&~Options.optimizePassiveJointEl&isfield(sol_val,'res_col_unsc')
    legend('total moment','sum muscle moments','residual');
elseif  ~Options.optimizeMuscleProp&Options.optimizePassiveJointEl&isfield(sol_val,'res_col_unsc')
    legend('total moment','sum joint passive moments','residual');
elseif Options.optimizeMuscleProp&Options.optimizePassiveJointEl&isfield(sol_val,'res_col_unsc')
    legend('total moment','sum muscle moments','sum joint passive moments','residual');
elseif Options.optimizeMuscleProp&~Options.optimizePassiveJointEl&~isfield(sol_val,'res_col_unsc')
    legend('total moment','sum muscle moments');
elseif  ~Options.optimizeMuscleProp&Options.optimizePassiveJointEl&~isfield(sol_val,'res_col_unsc')
    legend('total moment','sum joint passive moments');
elseif Options.optimizeMuscleProp&Options.optimizePassiveJointEl&~isfield(sol_val,'res_col_unsc')
    legend('total moment','sum muscle moments','sum joint passive moments');
end

    
    
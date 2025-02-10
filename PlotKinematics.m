function PlotKinematics(sol_val)


%plot kinematics
name_dofs={'hip flx', 'hip add', 'hip int', 'knee flex', 'ankle flex', 'ankle add', 'ankle int'};

for i=1:length(sol_val.QsQdots_unsc)
    if ~isempty(sol_val.QsQdots_unsc{i})
        figure(i*3-2)
        t=sol_val.tgrid{i};
        for j=1:7
            subplot(3,3,j);
            plot(t(sol_val.t2plot{i}), sol_val.guess.QsQdots{i}(sol_val.t2plot{i},(j-1)*2+1)*sol_val.scaling.q*180/pi,'LineWidth',2);
            hold all;
            % plot(t(sol_val.t2plot{i}), sol_val.QsQdots_col_unsc{i}((j-1)*2+1,:)*180/pi,'LineWidth',2);
            plot(t(1:4:end-1), sol_val.QsQdots_unsc{i}((j-1)*2+1,1:end-1)*180/pi,'LineWidth',2);
            title(name_dofs{j});
            ylabel('[º]')
            xlabel('time [s]')
        end
        legend({'experimental','model'})

        figure(i*3-1)
        t=sol_val.tgrid{i};
        for j=1:7
            subplot(3,3,j);
            plot(t(sol_val.t2plot{i}), sol_val.guess.QsQdots{i}(sol_val.t2plot{i},(j-1)*2+2)*sol_val.scaling.qdot,'LineWidth',2);
            hold all;
            % plot(t(sol_val.t2plot{i}), sol_val.QsQdots_col_unsc{i}((j-1)*2+2,:),'LineWidth',2);
            plot(t(1:4:end-1), sol_val.QsQdots_unsc{i}((j-1)*2+2,1:end-1),'LineWidth',2);

            title(name_dofs{j});
            ylabel('[rad/s]')
            xlabel('time [s]')
        end
        legend({'experimental','model'})

        figure(i*3)
        t=sol_val.tgrid{i};
        for j=1:7
            subplot(3,3,j);
            plot(t(sol_val.t2plot{i}), sol_val.guess.Qd2dots{i}(sol_val.t2plot{i},j)*sol_val.scaling.qd2dot,'LineWidth',2);
            hold all;
            plot(t(sol_val.t2plot{i}), sol_val.Qd2dot_col{i}(j,:)*sol_val.scaling.qd2dot,'LineWidth',2);
            title(name_dofs{j});
            ylabel('[rad/s^2]')
            xlabel('time [s]')
        end
        legend({'experimental','model'})

    end
end
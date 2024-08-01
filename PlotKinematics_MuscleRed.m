function PlotKinematics(sol_val)


%plot kinematics

for i=1:length(sol_val.QsQdots_col_unsc)
    if ~isempty(sol_val.QsQdots_col_unsc{i})
        f1=figure
        t=sol_val.tgrid{i};
        for j=1:7
            subplot(4,2,j);
            plot(sol_val.tgrid_col{i}, sol_val.QsQdots_col_unsc{i}(:,j*2-1)*180/pi);
            if j==6|j==7
                xlabel('time [s]')
            end
            title(sol_val.name_dofs{j})
            ylabel('angle [º]')
        end
        f2=figure
        for j=1:7
            subplot(4,2,j);
            plot(sol_val.tgrid_col{i}, sol_val.QsQdots_col_unsc{i}(:,j*2)*180/pi);
            if j==6|j==7
                xlabel('time [s]')
            end
            title(sol_val.name_dofs{j})
            ylabel('angle vel [º/s]')
        end
        
        
    end
end

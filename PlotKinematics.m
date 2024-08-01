function PlotKinematics(sol_val)


%plot kinematics

for i=1:length(sol_val.QsQdots_unsc)
    if ~isempty(sol_val.QsQdots_unsc{i})
        figure
        t=sol_val.tgrid{i};
        for j=1:7
            subplot(3,3,j);
            plot(t(sol_val.t2plot), sol_val.QsQdots_col_unsc{i}((j-1)*2+1,:)*180/pi);
            hold all;
            plot(t(sol_val.t2plot), sol_val.guess.QsQdots{i}(sol_val.t2plot,(j-1)*2+1)*180/pi);
        end
    end
end
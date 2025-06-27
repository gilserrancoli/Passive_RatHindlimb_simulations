function PlotKinematics_error(sol_val,guess,t_col_grid,Options)

for i=1:length(sol_val.QsQdots)
    errors(i,:)=ComputeErrors(sol_val.QsQdots_col_unsc{i},guess.QsQdots{i}(t_col_grid,:)');
end
bar(mean(errors));
hold all
errorbar(mean(errors),std(errors),'LineStyle','none','LineWidth',2);
labels={'hipl flex','hip add','hip int','knee flex','ankle flx','ankle add','ankle int'};
set(gca,'XTickLabels',labels)
ylabel('RMSE angles [º]');
end

function errors=ComputeErrors(QsQdots_col_unsc,expQsQdots)


errors=rms([QsQdots_col_unsc(1:2:end,:)-expQsQdots(1:2:end,:)]'*180/pi);
end
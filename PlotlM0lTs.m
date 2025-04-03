%Plot lM0
scaling=sol_val.scaling;

h=bar([sol_val.bounds.lM0.lower*scaling.lM0 sol_val.MTparam(2,:)', sol_val.lM0_unsc sol_val.bounds.lM0.upper*scaling.lM0]);
set(gca,'XTick',1:length(sol_val.muscle_names),'XTickLabels',sol_val.muscle_names);
ylabel('l^M_0 [m]')
set(h([1 4]),'FaceColor','k')
set(h([2]),'FaceColor','b')
set(h([3]),'FaceColor','r')
legend({'lb','ref','model','ub'})


figure
h=bar([sol_val.bounds.lTs.lower*scaling.lTs sol_val.lTs_unsc sol_val.bounds.lTs.upper*scaling.lTs]);
set(gca,'XTick',1:length(sol_val.muscle_names),'XTickLabels',sol_val.muscle_names);
ylabel('l^T_s [m]')
set(h([1 3]),'FaceColor','k')
set(h([2]),'FaceColor','r')
legend({'lb','model','ub'})


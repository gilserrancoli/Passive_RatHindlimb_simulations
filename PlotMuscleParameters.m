function PlotMuscleParameters(sol_val,bounds,scaling,muscle_names)

subplot(2,1,1);
h=bar([bounds.lM0.lower*scaling.lM0 sol_val.lM0_unsc bounds.lM0.upper*scaling.lM0]);
set([h(1) h(3)],'FaceColor',[0.5 0.5 0.5]);
set(gca,'XTick',1:length(muscle_names),'XTickLabels',muscle_names);
title('lM0')
ylabel('l^M_0 [m]');

subplot(2,1,2);
h=bar([bounds.lTs.lower*scaling.lTs sol_val.lTs_unsc bounds.lTs.upper*scaling.lTs]);
set([h(1) h(3)],'FaceColor',[0.5 0.5 0.5]);
set(gca,'XTick',1:length(muscle_names),'XTickLabels',muscle_names);
title('lTs')
ylabel('l^T_s [m]');
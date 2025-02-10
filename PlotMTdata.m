function PlotMTdata(sol_val,bounds,scaling,muscle_names,i)

    meanlMT=mean(sol_val.lMT{i});
    stdlMT=std(sol_val.lMT{i});
   
%     mat2comp=[bounds.lM0.lower*scaling.lM0+bounds.lTs.lower*scaling.lTs ...
%         sol_val.lM0_unsc+sol_val.lTs_unsc...
%         bounds.lM0.upper*scaling.lM0+bounds.lTs.upper*scaling.lTs meanlMT'];
%     subplot(2,1,1);
%     h=bar(mat2comp(1:19,:));
%     set(gca,'XTick',1:19,'XTickLabels',muscle_names(1:19));
%     subplot(2,1,2);
%     h2=bar(mat2comp(20:end,:));
%     set(gca,'XTick',1:19,'XTickLabels',muscle_names(20:end));
    
    mat2comp2=[bounds.lM0.lower*scaling.lM0 sol_val.lM0_unsc bounds.lM0.upper*scaling.lM0 ...
        bounds.lTs.lower*scaling.lTs sol_val.lTs_unsc bounds.lTs.upper*scaling.lTs meanlMT'];
    subplot(2,1,1);
    h=bar(mat2comp2(1:19,:));
    set(gca,'XTick',1:19,'XTickLabels',muscle_names(1:19));
    legend({'lM0_{LB}','lM0','lM0_{UB}','lTs_{LB}','lTs','lTs_{UB}','mean lMT'});
    ylabel('[m]');
    subplot(2,1,2);
    h2=bar(mat2comp2(20:end,:));
    set(gca,'XTick',1:19,'XTickLabels',muscle_names(20:end));
    ylabel('[m]');
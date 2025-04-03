%% Plot passive properties

%% Different solution for each trial
%Case 5 trials
x=[15 0 -15 0 15];

maxylimK(1:3)=1e-3;
maxylimD(1:3)=1e-3;
maxylimt(1:3)=0;
minylimt(1:3)=0;
for i=1:length(sol_val.nametrials)
    subplot(3,3,1);
    plot(x(i),sol_val.Kstiff{i}(1),'o','MarkerFaceColor',[0.00,0.45,0.74],'MarkerEdgeColor','none');
    title('Stiffness hip');
    set(gca,'Xtick',[-15 0 15]);
    maxylimK(1)=max(maxylimK(1),sol_val.Kstiff{i}(1));
    maxylimK(1)=ceil(maxylimK(1)*100)/100;
    ylim([0 maxylimK(1)]);
    ylabel('K [Nm/rad]');
    hold all;

    subplot(3,3,4);
    plot(x(i),sol_val.Kstiff{i}(2),'o','MarkerFaceColor',[0.00,0.45,0.74],'MarkerEdgeColor','none');
    title('Stiffness knee');
    set(gca,'Xtick',[-15 0 15]);
    maxylimK(2)=max(maxylimK(2),sol_val.Kstiff{i}(2));
    maxylimK(2)=ceil(maxylimK(2)*100)/100;
    ylim([0 maxylimK(2)]);
    ylabel('K [Nm/rad]');
    hold all;

    subplot(3,3,7);
    plot(x(i),sol_val.Kstiff{i}(3),'o','MarkerFaceColor',[0.00,0.45,0.74],'MarkerEdgeColor','none');
    title('Stiffness ankle');
    set(gca,'Xtick',[-15 0 15]);
    maxylimK(3)=max(maxylimK(3),sol_val.Kstiff{i}(3));
    maxylimK(3)=ceil(maxylimK(3)*100)/100;
    ylim([0 maxylimK(3)]);
    ylabel('K [Nm/rad]');
    hold all;


    %theta 0
    subplot(3,3,2);
    plot(x(i),sol_val.theta0{i}(1)*180/pi,'o','MarkerFaceColor',[0.00,0.45,0.74],'MarkerEdgeColor','none');
    title('theta_0');
    set(gca,'Xtick',[-15 0 15]);
    maxylimt(1)=max(maxylimt(1),sol_val.theta0{i}(1)*180/pi);
    maxylimt(1)=ceil(maxylimt(1)/10)*10;
    minylimt(1)=min(minylimt(1),sol_val.theta0{i}(1)*180/pi);
    minylimt(1)=floor(minylimt(1)/10)*10;
    ylim([minylimt(1) maxylimt(1)]);
    ylabel('theta_0 [º]');
    hold all;

    subplot(3,3,5);
    plot(x(i),sol_val.theta0{i}(2)*180/pi,'o','MarkerFaceColor',[0.00,0.45,0.74],'MarkerEdgeColor','none');
    hold all;
    title('theta_0');
    set(gca,'Xtick',[-15 0 15]);
    maxylimt(2)=max(maxylimt(2),sol_val.theta0{i}(2)*180/pi);
    maxylimt(2)=ceil(maxylimt(2)/10)*10;
    minylimt(2)=min(minylimt(2),sol_val.theta0{i}(2)*180/pi);
    minylimt(2)=floor(minylimt(2)/10)*10;
    ylim([minylimt(2) maxylimt(2)]);
    ylabel('theta_0 [º]');
    

    subplot(3,3,8);
    plot(x(i),sol_val.theta0{i}(3)*180/pi,'o','MarkerFaceColor',[0.00,0.45,0.74],'MarkerEdgeColor','none');
    title('theta_0');
    set(gca,'Xtick',[-15 0 15]);
    maxylimt(3)=max(maxylimt(3),sol_val.theta0{i}(3)*180/pi);
    maxylimt(3)=ceil(maxylimt(3)/10)*10;
    minylimt(3)=min(minylimt(3),sol_val.theta0{i}(3)*180/pi);
    minylimt(3)=floor(minylimt(3)/10)*10;
    ylim([minylimt(3) maxylimt(3)]);
    ylabel('theta_0 [º]');
    hold all;

    %Damping
    subplot(3,3,3);
    plot(x(i),sol_val.Ddamp{i}(1),'o','MarkerFaceColor',[0.00,0.45,0.74],'MarkerEdgeColor','none');
    title('Damping hip');
    set(gca,'Xtick',[-15 0 15]);
    maxylimD(1)=max(maxylimD(1),sol_val.Ddamp{i}(1));
    maxylimD(1)=ceil(maxylimD(1)*100)/100;
    ylim([0 maxylimD(1)]);
    ylabel('D [Nms/rad]');
    hold all;

    subplot(3,3,6);
    plot(x(i),sol_val.Ddamp{i}(2),'o','MarkerFaceColor',[0.00,0.45,0.74],'MarkerEdgeColor','none');
    title('Damping knee');
    set(gca,'Xtick',[-15 0 15]);
    maxylimD(2)=max(maxylimD(2),sol_val.Ddamp{i}(2));
    maxylimD(2)=ceil(maxylimD(2)*100)/100;
    ylim([0 maxylimD(2)]);
    ylabel('D [Nms/rad]');
    hold all;

    subplot(3,3,9);
    plot(x(i),sol_val.Ddamp{i}(3),'o','MarkerFaceColor',[0.00,0.45,0.74],'MarkerEdgeColor','none');
    title('Damping ankle');
    set(gca,'Xtick',[-15 0 15]);
    maxylimD(3)=max(maxylimD(3),sol_val.Ddamp{i}(3));
    maxylimD(3)=ceil(maxylimD(3)*100)/100;
    ylim([0 maxylimD(3)]);
    ylabel('D [Nms/rad]');
    hold all;

    

end

%Case 3 trials
x=[10 0 -10];

maxylimK(1:3)=1e-3;
maxylimD(1:3)=1e-3;
maxylimt(1:3)=0;
minylimt(1:3)=0;
for i=1:length(sol_val.nametrials)
    subplot(3,3,1);
    plot(x(i),sol_val.Kstiff{i}(1),'o','MarkerFaceColor',[0.00,0.45,0.74],'MarkerEdgeColor','none');
    title('Stiffness hip');
    set(gca,'Xtick',[-10 0 10]);
    xlim([-15 15]);
    maxylimK(1)=max(maxylimK(1),sol_val.Kstiff{i}(1));
    maxylimK(1)=ceil(maxylimK(1)*100)/100;
    ylim([0 maxylimK(1)]);
    ylabel('K [Nm/rad]');
    hold all;

    subplot(3,3,4);
    plot(x(i),sol_val.Kstiff{i}(2),'o','MarkerFaceColor',[0.00,0.45,0.74],'MarkerEdgeColor','none');
    title('Stiffness knee');
    set(gca,'Xtick',[-10 0 10]);
    xlim([-15 15]);
    maxylimK(2)=max(maxylimK(2),sol_val.Kstiff{i}(2));
    maxylimK(2)=ceil(maxylimK(2)*100)/100;
    ylim([0 maxylimK(2)]);
    ylabel('K [Nm/rad]');
    hold all;

    subplot(3,3,7);
    plot(x(i),sol_val.Kstiff{i}(3),'o','MarkerFaceColor',[0.00,0.45,0.74],'MarkerEdgeColor','none');
    title('Stiffness ankle');
    set(gca,'Xtick',[-10 0 10]);
    xlim([-15 15]);
    maxylimK(3)=max(maxylimK(3),sol_val.Kstiff{i}(3));
    maxylimK(3)=ceil(maxylimK(3)*100)/100;
    ylim([0 maxylimK(3)]);
    ylabel('K [Nm/rad]');
    hold all;


    %theta 0
    subplot(3,3,2);
    plot(x(i),sol_val.theta0{i}(1)*180/pi,'o','MarkerFaceColor',[0.00,0.45,0.74],'MarkerEdgeColor','none');
    title('theta_0');
    set(gca,'Xtick',[-10 0 10]);
    xlim([-15 15]);
    maxylimt(1)=max(maxylimt(1),sol_val.theta0{i}(1)*180/pi);
    maxylimt(1)=ceil(maxylimt(1)/10)*10;
    minylimt(1)=min(minylimt(1),sol_val.theta0{i}(1)*180/pi);
    minylimt(1)=floor(minylimt(1)/10)*10;
    ylim([minylimt(1) maxylimt(1)]);
    ylabel('theta_0 [º]');
    hold all;

    subplot(3,3,5);
    plot(x(i),sol_val.theta0{i}(2)*180/pi,'o','MarkerFaceColor',[0.00,0.45,0.74],'MarkerEdgeColor','none');
    hold all;
    title('theta_0');
    set(gca,'Xtick',[-10 0 10]);
    xlim([-15 15]);
    maxylimt(2)=max(maxylimt(2),sol_val.theta0{i}(2)*180/pi);
    maxylimt(2)=ceil(maxylimt(2)/10)*10;
    minylimt(2)=min(minylimt(2),sol_val.theta0{i}(2)*180/pi);
    minylimt(2)=floor(minylimt(2)/10)*10;
    ylim([minylimt(2) maxylimt(2)]);
    ylabel('theta_0 [º]');
    

    subplot(3,3,8);
    plot(x(i),sol_val.theta0{i}(3)*180/pi,'o','MarkerFaceColor',[0.00,0.45,0.74],'MarkerEdgeColor','none');
    title('theta_0');
    set(gca,'Xtick',[-10 0 10]);
    xlim([-15 15]);
    maxylimt(3)=max(maxylimt(3),sol_val.theta0{i}(3)*180/pi);
    maxylimt(3)=ceil(maxylimt(3)/10)*10;
    minylimt(3)=min(minylimt(3),sol_val.theta0{i}(3)*180/pi);
    minylimt(3)=floor(minylimt(3)/10)*10;
    ylim([minylimt(3) maxylimt(3)]);
    ylabel('theta_0 [º]');
    hold all;

    %Damping
    subplot(3,3,3);
    plot(x(i),sol_val.Ddamp{i}(1),'o','MarkerFaceColor',[0.00,0.45,0.74],'MarkerEdgeColor','none');
    title('Damping hip');
    set(gca,'Xtick',[-10 0 10]);
    xlim([-15 15]);
    maxylimD(1)=max(maxylimD(1),sol_val.Ddamp{i}(1));
    maxylimD(1)=ceil(maxylimD(1)*100)/100;
    ylim([0 maxylimD(1)]);
    ylabel('D [Nms/rad]');
    hold all;

    subplot(3,3,6);
    plot(x(i),sol_val.Ddamp{i}(2),'o','MarkerFaceColor',[0.00,0.45,0.74],'MarkerEdgeColor','none');
    title('Damping knee');
    set(gca,'Xtick',[-10 0 10]);
    xlim([-15 15]);
    maxylimD(2)=max(maxylimD(2),sol_val.Ddamp{i}(2));
    maxylimD(2)=ceil(maxylimD(2)*100)/100;
    ylim([0 maxylimD(2)]);
    ylabel('D [Nms/rad]');
    hold all;

    subplot(3,3,9);
    plot(x(i),sol_val.Ddamp{i}(3),'o','MarkerFaceColor',[0.00,0.45,0.74],'MarkerEdgeColor','none');
    title('Damping ankle');
    set(gca,'Xtick',[-10 0 10]);
    xlim([-15 15]);
    maxylimD(3)=max(maxylimD(3),sol_val.Ddamp{i}(3));
    maxylimD(3)=ceil(maxylimD(3)*100)/100;
    ylim([0 maxylimD(3)]);
    ylabel('D [Nms/rad]');
    hold all;

    

end








%% One solution for all trials
subplot(3,1,1)
bar(sol_val.Kstiff);
set(gca,'XTick',1:3,'XTickLabels',{'hip','knee','ankle'});
ylabel('K [Nm/rad]');
subplot(3,1,2)
bar(sol_val.theta0);
set(gca,'XTick',1:3,'XTickLabels',{'hip','knee','ankle'});
ylabel('\theta_0 [rad]');
subplot(3,1,3)
bar(sol_val.Ddamp);
set(gca,'XTick',1:3,'XTickLabels',{'hip','knee','ankle'});
ylabel('D [Nm/(rad/s)]');
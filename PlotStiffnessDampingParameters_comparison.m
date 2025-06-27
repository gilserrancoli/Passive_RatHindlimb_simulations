S1=load('sol_val_optJointPassive_optInertia_DataMay2025IntactFB_tolem5_JpenKDT1e-3_JminKDT1e-5_JminInertiaP1e-2_rel2Dangles_N40.mat');
S2=load('sol_val_optJointPassive_optInertia_DataMay2025achilelscutFB_tolem5_JpenKDT1e-3_JminKDT1e-5_JminInertiaP1e-2_rel2Dangles_N40.mat');

%% case 7+7 trials
% x=[30 15 0 -15 0 15 30];
x=obtainPerturbIDs(S1.sol_val.nametrials);
maxylimK(1:3)=0.005;
maxylimD(1:3)=0.03;
maxylimt(1:3)=0;
maxylimt(2)=-40;
minylimt(1:3)=0;
minylimt(2)=-60;
maxylimI(1:3)=1e-5;
if Options.individualpassiveprop==1
    ntrials4passprop=size(S1.sol_val.nametrials,1);
    list4passprop=[];
elseif contains(main_folder,'May2025')
    ntrials4passprop=4;
    list4passprop=[4 3 2 1 1 2 3 4 4 3 2 1 1 2 3 4]; %perturbations at 30 20 10 0 0 10...
    if size(S1.sol_val.Kstiff_unsc,2)==1
        S1.sol_val.Kstiff_unsc=reshape(S1.sol_val.Kstiff_unsc,3,4);
        S2.sol_val.Kstiff_unsc=reshape(S2.sol_val.Kstiff_unsc,3,4);
        S1.sol_val.Ddamp_unsc=reshape(S1.sol_val.Ddamp_unsc,3,4);
        S2.sol_val.Ddamp_unsc=reshape(S2.sol_val.Ddamp_unsc,3,4);
        S1.sol_val.theta0_unsc=reshape(S1.sol_val.theta0_unsc,3,4);
        S2.sol_val.theta0_unsc=reshape(S2.sol_val.theta0_unsc,3,4);
        if Options.optInertiapassiveParam
            S1.sol_val.InertiapassiveParam_unsc=reshape(S1.sol_val.InertiapassiveParam_unsc,3,4);
            S2.sol_val.InertiapassiveParam_unsc=reshape(S2.sol_val.InertiapassiveParam_unsc,3,4);
        end
    end
else
    %define how many different passive prop are considered
    keyboard;
end

for i=1:size(S1.sol_val.Kstiff_unsc,2)
    if i<=7
        X=x(i);
        FA=1;
    else
        X=x(i-7);
        FA=0.5;
    end
    MFC1=[0.00,0.45,0.74];
    MFC2=[1,0.45,0.24];

    hout1(i)=PlotTrials(X,x,S1.sol_val,MFC1,i,maxylimK,maxylimD,maxylimt,minylimt,FA,ntrials4passprop,S1.sol_val.Options,list4passprop);
    hout2(i)=PlotTrials(X,x,S2.sol_val,MFC2,i,maxylimK,maxylimD,maxylimt,minylimt,FA,ntrials4passprop,S2.sol_val.Options,list4passprop);
end
if Options.individualpassiveprop
    legend([hout1(1), hout1(8),hout2(1),hout2(8)],{'Intact For','Intact Bac','AchillesCut For','AchillesCut Bac'},'Orientation','horizontal','box','off')
else
    legend([hout1(1), hout2(1)],{'Intact','Achilles Cut'},'Orientation','horizontal','box','off')
end

if isfield(S1.sol_val.Options,'optInertiapassiveParam')
    if S1.sol_val.Options.optInertiapassiveParam
        figure
        for i=1:size(S1.sol_val.InertiapassiveParam_unsc,2)
            if i<=7
                X=x(i);
                FA=1;
            else
                X=x(i-7);
                FA=0.5;
            end
        MFC1=[0.00,0.45,0.74];
        MFC2=[1,0.45,0.24];
    
        hout1(i)=PlotTrials_Inertia(X,x,S1.sol_val,MFC1,i,maxylimI,FA,ntrials4passprop,S1.sol_val.Options,list4passprop);
        hout2(i)=PlotTrials_Inertia(X,x,S2.sol_val,MFC2,i,maxylimI,FA,ntrials4passprop,S2.sol_val.Options,list4passprop);
    
        end
    end
end



function hout=PlotTrials(X,x,sol_val,MFC,i,maxylimK,maxylimD,maxylimt,minylimt,FA,ntrials4passprop,Options,list4passprop)
xticks=unique(sort(x));    

    subplot(3,3,1);
    if Options.individualpassiveprop
        Itrial4passprop=i;
        scatter(X,sol_val.Kstiff_unsc{i}(1),'o','MarkerFaceColor',MFC,'MarkerEdgeColor','none','MarkerFaceAlpha',FA);
    else
        Itrial4passprop=list4passprop(i);
        scatter(X,sol_val.Kstiff_unsc(1,Itrial4passprop),'o','MarkerFaceColor',MFC,'MarkerEdgeColor','none','MarkerFaceAlpha',FA);
    end
    title('Stiffness hip');
    set(gca,'Xtick',xticks);
    if i>1
        y=ylim;
        if Options.individualpassiveprop
            maxylimK(1)=max([y(2),maxylimK(1),sol_val.Kstiff_unsc{i}(1)]);
        else
            maxylimK(1)=max([y(2),maxylimK(1),sol_val.Kstiff_unsc(1,Itrial4passprop)]);
        end
    else
        if Options.individualpassiveprop
            maxylimK(1)=max(maxylimK(1),sol_val.Kstiff_unsc{i}(1));
        else
            maxylimK(1)=max(maxylimK(1),sol_val.Kstiff_unsc(1,Itrial4passprop));
        end
    end
    maxylimK(1)=ceil(maxylimK(1)*100)/100;
    ylim([0 maxylimK(1)]);
    ylabel('K [Nm/rad]');
    hold all;
    xlim([min(x) max(x)]);

    subplot(3,3,4);
    if Options.individualpassiveprop
        Itrial4passprop=i;
        scatter(X,sol_val.Kstiff_unsc{i}(2),'o','MarkerFaceColor',MFC,'MarkerEdgeColor','none','MarkerFaceAlpha',FA);
    else
        Itrial4passprop=list4passprop(i);
        scatter(X,sol_val.Kstiff_unsc(2,Itrial4passprop),'o','MarkerFaceColor',MFC,'MarkerEdgeColor','none','MarkerFaceAlpha',FA);
    end
    title('Stiffness knee');
    if i>1
        y=ylim;
        if Options.individualpassiveprop
            maxylimK(2)=max([y(2),maxylimK(2),sol_val.Kstiff_unsc{i}(2)]);
        else
            maxylimK(2)=max([y(2),maxylimK(2),sol_val.Kstiff_unsc(2,Itrial4passprop)]);
        end
    else
        if Options.individualpassiveprop
            maxylimK(2)=max(maxylimK(2),sol_val.Kstiff_unsc{i}(2));
        else
            maxylimK(2)=max(maxylimK(2),sol_val.Kstiff_unsc(2,Itrial4passprop));
        end
    end
    set(gca,'Xtick',xticks);
    maxylimK(2)=ceil(maxylimK(2)*100)/100;
    ylim([0 maxylimK(2)]);
    ylabel('K [Nm/rad]');
    hold all;
    xlim([min(x) max(x)]);

    subplot(3,3,7);
    if Options.individualpassiveprop
        Itrial4passprop=i;
        scatter(X,sol_val.Kstiff_unsc{i}(3),'o','MarkerFaceColor',MFC,'MarkerEdgeColor','none','MarkerFaceAlpha',FA);
    else
        Itrial4passprop=list4passprop(i);
        scatter(X,sol_val.Kstiff_unsc(3,Itrial4passprop),'o','MarkerFaceColor',MFC,'MarkerEdgeColor','none','MarkerFaceAlpha',FA);
    end
    title('Stiffness ankle');
    if i>1
        y=ylim;
        if Options.individualpassiveprop
            maxylimK(3)=max([y(2),maxylimK(3),sol_val.Kstiff_unsc{i}(3)]);
        else
            maxylimK(3)=max([y(2),maxylimK(3),sol_val.Kstiff_unsc(3,Itrial4passprop)]);
        end
    else
        if Options.individualpassiveprop
            maxylimK(3)=max(maxylimK(3),sol_val.Kstiff_unsc{i}(3));
        else
            maxylimK(3)=max(maxylimK(3),sol_val.Kstiff_unsc(3,Itrial4passprop));            
        end
    end
    set(gca,'Xtick',xticks);
    maxylimK(3)=ceil(maxylimK(3)*100)/100;
    ylim([0 maxylimK(3)]);
    ylabel('K [Nm/rad]');
    hold all;
    xlim([min(x) max(x)]);

    %theta 0
    subplot(3,3,2);
    if i>1
        y=ylim;
    end
    if Options.individualpassiveprop
        Itrial4passprop=i;
        scatter(X,sol_val.theta0_unsc{i}(1)*180/pi,'o','MarkerFaceColor',MFC,'MarkerEdgeColor','none','MarkerFaceAlpha',FA);
        if i>1
            maxylimt(1)=max([y(2),maxylimt(1),sol_val.theta0_unsc{i}(1)*180/pi]);
            minylimt(1)=min([y(1),minylimt(1),sol_val.theta0_unsc{i}(1)*180/pi]);
        else
            maxylimt(1)=max(maxylimt(1),sol_val.theta0_unsc{i}(1)*180/pi);
            minylimt(1)=min(minylimt(1),sol_val.theta0_unsc{i}(1)*180/pi);
        end
    else
        Itrial4passprop=list4passprop(i);
        scatter(X,sol_val.theta0_unsc(1,Itrial4passprop)*180/pi,'o','MarkerFaceColor',MFC,'MarkerEdgeColor','none','MarkerFaceAlpha',FA);
        if i>1
            maxylimt(1)=max([y(2),maxylimt(1),sol_val.theta0_unsc(1,Itrial4passprop)*180/pi]);
            minylimt(1)=min([y(1),minylimt(1),sol_val.theta0_unsc(1,Itrial4passprop)*180/pi]);
        else
            maxylimt(1)=max(maxylimt(1),sol_val.theta0_unsc(1,Itrial4passprop)*180/pi);
            minylimt(1)=min(minylimt(1),sol_val.theta0_unsc(1,Itrial4passprop)*180/pi);
        end
    end
    title('theta_0');
    set(gca,'Xtick',xticks);
    maxylimt(1)=ceil(maxylimt(1)/10)*10;
    minylimt(1)=floor(minylimt(1)/10)*10;
    ylim([minylimt(1) maxylimt(1)]);
    ylabel('theta_0 [º]');
    hold all;
    xlim([min(x) max(x)]);

    subplot(3,3,5);
    if Options.individualpassiveprop
        scatter(X,sol_val.theta0_unsc{i}(2)*180/pi,'o','MarkerFaceColor',MFC,'MarkerEdgeColor','none','MarkerFaceAlpha',FA);
        maxylimt(2)=max(maxylimt(2),sol_val.theta0_unsc{i}(2)*180/pi);
        minylimt(2)=min(minylimt(2),sol_val.theta0_unsc{i}(2)*180/pi);
    else
        Itrial4passprop=list4passprop(i);
        scatter(X,sol_val.theta0_unsc(2,Itrial4passprop)*180/pi,'o','MarkerFaceColor',MFC,'MarkerEdgeColor','none','MarkerFaceAlpha',FA);
        maxylimt(2)=max(maxylimt(2),sol_val.theta0_unsc(2,Itrial4passprop)*180/pi);
        minylimt(2)=min(minylimt(2),sol_val.theta0_unsc(2,Itrial4passprop)*180/pi);
    end
    hold all;
    title('theta_0');
    set(gca,'Xtick',xticks);
    maxylimt(2)=ceil(maxylimt(2)/10)*10;
    minylimt(2)=floor(minylimt(2)/10)*10;
    ylim([minylimt(2) maxylimt(2)]);
    ylabel('theta_0 [º]');
    xlim([min(x) max(x)]);

    subplot(3,3,8);
    if i>1
        y=ylim;
    end
    if Options.individualpassiveprop
        scatter(X,sol_val.theta0_unsc{i}(3)*180/pi,'o','MarkerFaceColor',MFC,'MarkerEdgeColor','none','MarkerFaceAlpha',FA);
        if i>1
            maxylimt(3)=max([y(2),maxylimt(3),sol_val.theta0_unsc{i}(3)*180/pi]);
            minylimt(3)=min([y(1),minylimt(3),sol_val.theta0_unsc{i}(3)*180/pi]);
        else
            maxylimt(3)=max(maxylimt(3),sol_val.theta0_unsc{i}(3)*180/pi);
            minylimt(3)=min(minylimt(3),sol_val.theta0_unsc{i}(3)*180/pi);
        end
    else
        Itrial4passprop=list4passprop(i);
        scatter(X,sol_val.theta0_unsc(3,Itrial4passprop)*180/pi,'o','MarkerFaceColor',MFC,'MarkerEdgeColor','none','MarkerFaceAlpha',FA);
        if i>1
            maxylimt(3)=max([y(2),maxylimt(3),sol_val.theta0_unsc(3,Itrial4passprop)*180/pi]);
            minylimt(3)=min([y(1),minylimt(3),sol_val.theta0_unsc(3,Itrial4passprop)*180/pi]);
        else
            maxylimt(3)=max(maxylimt(3),sol_val.theta0_unsc(3,Itrial4passprop)*180/pi);
            minylimt(3)=min(minylimt(3),sol_val.theta0_unsc(3,Itrial4passprop)*180/pi);
        end
    end
    title('theta_0');
    set(gca,'Xtick',xticks);
    maxylimt(3)=ceil(maxylimt(3)/10)*10;
    minylimt(3)=floor(minylimt(3)/10)*10;
    ylim([minylimt(3) maxylimt(3)]);
    ylabel('theta_0 [º]');
    hold all;
    xlim([min(x) max(x)]);

    %Damping
    subplot(3,3,3);
    if Options.individualpassiveprop
        scatter(X,sol_val.Ddamp_unsc{i}(1),'o','MarkerFaceColor',MFC,'MarkerEdgeColor','none','MarkerFaceAlpha',FA);
    else
        Itrial4passprop=list4passprop(i);
        scatter(X,sol_val.Ddamp_unsc(1,Itrial4passprop),'o','MarkerFaceColor',MFC,'MarkerEdgeColor','none','MarkerFaceAlpha',FA);
    end
    title('Damping hip');
    if i>1
        y=ylim;
        if Options.individualpassiveprop
            maxylimD(1)=max([y(2),maxylimD(1),sol_val.Ddamp_unsc{i}(1)]);
        else
            maxylimD(1)=max([y(2),maxylimD(1),sol_val.Ddamp_unsc(1,Itrial4passprop)]);
        end
    else
        if Options.individualpassiveprop
            maxylimD(1)=max(maxylimD(1),sol_val.Ddamp_unsc{i}(1));
        else
            maxylimD(1)=max(maxylimD(1),sol_val.Ddamp_unsc(1,Itrial4passprop));
        end
    end
    set(gca,'Xtick',xticks);
    maxylimD(1)=ceil(maxylimD(1)*100)/100;
    ylim([0 maxylimD(1)]);
    ylabel('D [Nms/rad]');
    hold all;
    xlim([min(x) max(x)]);

    subplot(3,3,6);
    if Options.individualpassiveprop
        scatter(X,sol_val.Ddamp_unsc{i}(2),'o','MarkerFaceColor',MFC,'MarkerEdgeColor','none','MarkerFaceAlpha',FA);
    else
        Itrial4passprop=list4passprop(i);
        scatter(X,sol_val.Ddamp_unsc(2,Itrial4passprop),'o','MarkerFaceColor',MFC,'MarkerEdgeColor','none','MarkerFaceAlpha',FA);
    end
    title('Damping knee');
    if i>1
        y=ylim;
        if Options.individualpassiveprop
            maxylimD(2)=max([y(2),maxylimD(2),sol_val.Ddamp_unsc{i}(2)]);
        else
            maxylimD(2)=max([y(2),maxylimD(2),sol_val.Ddamp_unsc(2,Itrial4passprop)]);
        end
    else
        if Options.individualpassiveprop
            maxylimD(2)=max(maxylimD(2),sol_val.Ddamp_unsc{i}(2));
        else
            maxylimD(2)=max(maxylimD(2),sol_val.Ddamp_unsc(2,Itrial4passprop));
        end
    end
    set(gca,'Xtick',xticks);
    maxylimD(2)=ceil(maxylimD(2)*100)/100;
    ylim([0 maxylimD(2)]);
    ylabel('D [Nms/rad]');
    hold all;
    xlim([min(x) max(x)]);

    subplot(3,3,9);
    if Options.individualpassiveprop
        hout=scatter(X,sol_val.Ddamp_unsc{i}(3),'o','MarkerFaceColor',MFC,'MarkerEdgeColor','none','MarkerFaceAlpha',FA);
    else
        Itrial4passprop=list4passprop(i);
        hout=scatter(X,sol_val.Ddamp_unsc(3,Itrial4passprop),'o','MarkerFaceColor',MFC,'MarkerEdgeColor','none','MarkerFaceAlpha',FA);
    end
    title('Damping ankle');
    if i>1
        y=ylim;
        if Options.individualpassiveprop
            maxylimD(3)=max([y(2),maxylimD(3),sol_val.Ddamp_unsc{i}(3)]);
        else
            maxylimD(3)=max([y(2),maxylimD(3),sol_val.Ddamp_unsc(3,Itrial4passprop)]);
        end
    else
        if Options.individualpassiveprop
            maxylimD(3)=max(maxylimD(3),sol_val.Ddamp_unsc{i}(3));
        else
            maxylimD(3)=max(maxylimD(3),sol_val.Ddamp_unsc(3,Itrial4passprop));
        end
    end
    set(gca,'Xtick',xticks);
    maxylimD(3)=ceil(maxylimD(3)*100)/100;
    ylim([0 maxylimD(3)]);
    ylabel('D [Nms/rad]');
    hold all;
    xlim([min(x) max(x)]);
end

function hout=PlotTrials_Inertia(X,x,sol_val,MFC,i,maxylimI,FA,ntrials4passprop,Options,list4passprop)
    xticks=unique(sort(x));    

    subplot(3,1,1);
    if Options.individualpassiveprop
        Itrial4passprop=i;
        scatter(X,sol_val.InertiapassiveParam_unsc{i}(1),'o','MarkerFaceColor',MFC,'MarkerEdgeColor','none','MarkerFaceAlpha',FA);
    else
        Itrial4passprop=list4passprop(i);
        scatter(X,sol_val.InertiapassiveParam_unsc(1,Itrial4passprop),'o','MarkerFaceColor',MFC,'MarkerEdgeColor','none','MarkerFaceAlpha',FA);
    end
    title('"Inertia" hip');
    set(gca,'Xtick',xticks);
    if i>1
        y=ylim;
        if Options.individualpassiveprop
            maxylimI(1)=max([y(2),maxylimI(1),sol_val.InertiapassiveParam_unsc{i}(1)]);
        else
            maxylimI(1)=max([y(2),maxylimI(1),sol_val.InertiapassiveParam_unsc(1,Itrial4passprop)]);
        end
    else
        if Options.individualpassiveprop
            maxylimI(1)=max(maxylimI(1),sol_val.InertiapassiveParam_unsc{i}(1));
        else
            maxylimI(1)=max(maxylimI(1),sol_val.InertiapassiveParam_unsc(1,Itrial4passprop));
        end
    end
    maxylimI(1)=ceil(maxylimI(1)*100000)/100000;
    ylim([0 maxylimI(1)]);
    ylabel('I [Nm/rad^s]');
    hold all;
    xlim([min(x) max(x)]);


    subplot(3,1,2);
    if Options.individualpassiveprop
        Itrial4passprop=i;
        scatter(X,sol_val.InertiapassiveParam_unsc{i}(2),'o','MarkerFaceColor',MFC,'MarkerEdgeColor','none','MarkerFaceAlpha',FA);
    else
        Itrial4passprop=list4passprop(i);
        scatter(X,sol_val.InertiapassiveParam_unsc(2,Itrial4passprop),'o','MarkerFaceColor',MFC,'MarkerEdgeColor','none','MarkerFaceAlpha',FA);
    end
    title('"Inertia" knee');
    set(gca,'Xtick',xticks);
    if i>1
        y=ylim;
        if Options.individualpassiveprop
            maxylimI(2)=max([y(2),maxylimI(2),sol_val.InertiapassiveParam_unsc{i}(2)]);
        else
            maxylimI(2)=max([y(2),maxylimI(2),sol_val.InertiapassiveParam_unsc(2,Itrial4passprop)]);
        end
    else
        if Options.individualpassiveprop
            maxylimI(2)=max(maxylimI(2),sol_val.InertiapassiveParam_unsc{i}(2));
        else
            maxylimI(2)=max(maxylimI(2),sol_val.InertiapassiveParam_unsc(2,Itrial4passprop));
        end
    end
    maxylimI(2)=ceil(maxylimI(2)*100000)/100000;
    ylim([0 maxylimI(2)]);
    ylabel('I [Nm/rad^s]');
    hold all;
    xlim([min(x) max(x)]);

    subplot(3,1,3);
    if Options.individualpassiveprop
        Itrial4passprop=i;
        hout=scatter(X,sol_val.InertiapassiveParam_unsc{i}(3),'o','MarkerFaceColor',MFC,'MarkerEdgeColor','none','MarkerFaceAlpha',FA);
    else
        Itrial4passprop=list4passprop(i);
        hout=scatter(X,sol_val.InertiapassiveParam_unsc(3,Itrial4passprop),'o','MarkerFaceColor',MFC,'MarkerEdgeColor','none','MarkerFaceAlpha',FA);
    end
    title('"Inertia" ankle');
    set(gca,'Xtick',xticks);
    if i>1
        y=ylim;
        if Options.individualpassiveprop
            maxylimI(3)=max([y(2),maxylimI(3),sol_val.InertiapassiveParam_unsc{i}(3)]);
        else
            maxylimI(3)=max([y(2),maxylimI(3),sol_val.InertiapassiveParam_unsc(3,Itrial4passprop)]);
        end
    else
        if Options.individualpassiveprop
            maxylimI(3)=max(maxylimI(3),sol_val.InertiapassiveParam_unsc{i}(3));
        else
            maxylimI(3)=max(maxylimI(3),sol_val.InertiapassiveParam_unsc(3,Itrial4passprop));
        end
    end
    maxylimI(3)=ceil(maxylimI(3)*100000)/100000;
    ylim([0 maxylimI(3)]);
    ylabel('I [Nm/rad^s]');
    hold all;
    xlim([min(x) max(x)]);


end

function x=obtainPerturbIDs(nametrials)

    for i=1:length(nametrials)
        S=strsplit(nametrials{i},'_');
        pert=S{3};
        pert=strrep(pert,'pos','');
        if strcmp(pert(1),'m')
            x(i)=-str2num(pert(2:end));
        else
            x(i)=str2num(pert(1:end));
        end
    end


end
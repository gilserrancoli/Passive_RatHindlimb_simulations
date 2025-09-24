%% Lengths, velocities and accelerations from Muscle Analysis
for i=1:12
    if ~contains(folders(i).name,'.')
        cd(folders(i).name)
folder='C:\Gil\Collaborations\MatthewTresch\Zhong_ParameterEstimation\Optimization\MuscleAnalysisPointKin_formuscleacc_August2025\1_forward10mm_at20';
lengths=importdata([folder '\RatHindlimbRight-scaled_MuscleAnalysis_Length.sto']);

subplot(3,4,1)
plot(lengths.data(:,1),lengths.data(:,30));
ylabel('length [m]')
hold all;
title('LG');
subplot(3,4,2)
plot(lengths.data(:,1),lengths.data(:,29));
hold all;
title('MG');
subplot(3,4,3)
plot(lengths.data(:,1),lengths.data(:,31));
hold all;
title('Pla');
subplot(3,4,4)
plot(lengths.data(:,1),lengths.data(:,32));
hold all;
title('Sol');

[B,A]=butter(3,1000/5000);
lengths.data(:,2:end)=filtfilt(B,A,lengths.data(:,2:end));

subplot(3,4,1)
plot(lengths.data(:,1),lengths.data(:,30));
ylabel('length [m]')
hold all;
title('LG');
subplot(3,4,2)
plot(lengths.data(:,1),lengths.data(:,29));
title('MG');
subplot(3,4,3)
plot(lengths.data(:,1),lengths.data(:,31));
title('Pla');
subplot(3,4,4)
plot(lengths.data(:,1),lengths.data(:,32));
title('Sol');


splinesM(1)=spline(lengths.data(:,1),lengths.data(:,30));
splinesM(2)=spline(lengths.data(:,1),lengths.data(:,29));
splinesM(3)=spline(lengths.data(:,1),lengths.data(:,31));
splinesM(4)=spline(lengths.data(:,1),lengths.data(:,32));
for i=1:4
    vsplinesM(i)=fnder(splinesM(i),1);
    asplinesM(i)=fnder(splinesM(i),2);
end

t=lengths.data(:,1);

subplot(3,4,5)
plot(lengths.data(:,1),ppval(vsplinesM(1),t));
ylabel('vel [m/s]')
subplot(3,4,6)
plot(lengths.data(:,1),ppval(vsplinesM(2),t));
subplot(3,4,7)
plot(lengths.data(:,1),ppval(vsplinesM(3),t));
subplot(3,4,8)
plot(lengths.data(:,1),ppval(vsplinesM(4),t));

subplot(3,4,9)
plot(lengths.data(:,1),ppval(asplinesM(1),t));
ylabel('acc [m/s^2]');
xlabel('time [s]');
subplot(3,4,10)
plot(lengths.data(:,1),ppval(asplinesM(2),t));
xlabel('time [s]');
subplot(3,4,11)
plot(lengths.data(:,1),ppval(asplinesM(3),t));
xlabel('time [s]');
subplot(3,4,12)
plot(lengths.data(:,1),ppval(asplinesM(4),t));
xlabel('time [s]');

%% Position, velocity and acceleration of the center point of the muscle
pos_data(1)=importdata([folder '\RatHindlimbRight-scaled_PointKinematics_LGorigin_pos.sto']);
pos_data(2)=importdata([folder '\RatHindlimbRight-scaled_PointKinematics_LGinsertion_pos.sto']);
pos_data(3)=importdata([folder '\RatHindlimbRight-scaled_PointKinematics_MGorigin_pos.sto']);
pos_data(4)=importdata([folder '\RatHindlimbRight-scaled_PointKinematics_MGinsertion_pos.sto']);
pos_data(5)=importdata([folder '\RatHindlimbRight-scaled_PointKinematics_Plaorigin_pos.sto']);
pos_data(6)=importdata([folder '\RatHindlimbRight-scaled_PointKinematics_Plainsertion_pos.sto']);
pos_data(7)=importdata([folder '\RatHindlimbRight-scaled_PointKinematics_Solorigin_pos.sto']);
pos_data(8)=importdata([folder '\RatHindlimbRight-scaled_PointKinematics_Solinsertion_pos.sto']);

t=pos_data(1).data(:,1);
posvec_LG=pos_data(2).data(:,2:4)-pos_data(1).data(:,2:4);
posvec_MG=pos_data(4).data(:,2:4)-pos_data(3).data(:,2:4);
posvec_Pla=pos_data(6).data(:,2:4)-pos_data(5).data(:,2:4);
posvec_Sol=pos_data(8).data(:,2:4)-pos_data(7).data(:,2:4);

for i=1:3
    possplines(i,1)=spline(t,posvec_LG(:,i));
    possplines(i,2)=spline(t,posvec_MG(:,i));
    possplines(i,3)=spline(t,posvec_Pla(:,i));
    possplines(i,4)=spline(t,posvec_Sol(:,i));

    for j=1:4
        vsplines(i,j)=fnder(possplines(i,j),1);
        asplines(i,j)=fnder(possplines(i,j),2);
    end
end

%plot components of acceleration and resultant value
for j=1:4
    subplot(1,4,j)
    for i=1:3
        switch i
            case 1
                c='r';
            case 2
                c='g';
            case 3
                c='b';
        end
        plot(t,ppval(asplines(i,j),t),'Color',c);
        hold all;
    end
    plot(t,sqrt(sum(ppval(asplines(1,j),t).^2+...
        ppval(asplines(2,j),t).^2+...
        ppval(asplines(3,j),t).^2,2)),'k','LineWidth',2);
    if j==1
        title('LG');
    elseif j==2
        title('MG');
    elseif j==3
        title('Pla');
    elseif j==4
        title('Sol');
    end
    xlabel('time [s]');
    ylabel('acc ([m/s^2]');
end
legend({'a_x','a_y','a_z','|a|'})


%% Plot moment arms
ma_ankle=importdata([folder '\RatHindlimbRight-scaled_MuscleAnalysis_MomentArm_ankle_flx.sto']);
plot(ma_ankle.data(:,1),ma_ankle.data(:,[30 29 31 32]));
xlabel('time');
ylabel('ma [m]');
legend({'LG', 'MG','Pla','Sol'});

%% Plot "muscle moment due to mass acceleration"
mass=2e-3;
MLG=-mass*ppval(asplinesM(1),t).*ma_ankle.data(:,30);
MMG=-mass*ppval(asplinesM(2),t).*ma_ankle.data(:,29);
MPla=-mass*ppval(asplinesM(3),t).*ma_ankle.data(:,31);
MSol=-mass*ppval(asplinesM(4),t).*ma_ankle.data(:,32);

plot(t,MLG)
hold all
plot(t,MMG)
plot(t,MPla)
plot(t,MSol)
plot(t,MLG+MMG+MPla+MSol,'k','LineWidth',2);
title('total -m*acc*ma');
xlabel('time');
ylabel('moment [Nm]')
legend({'LG','MG','Pla','Sol','Total'})


figure(64)
    plot(t,MLG+MMG+MPla+MSol);
    hold all;

    else
    end
    
end
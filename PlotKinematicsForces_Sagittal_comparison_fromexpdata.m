folder_baseline='C:\Gil\Collaborations\MatthewTresch\Zhong_ParameterEstimation\Optimization\DataAugust2025\baseline_forward_2mm\';
folder_achillescut='C:\Gil\Collaborations\MatthewTresch\Zhong_ParameterEstimation\Optimization\DataAugust2025\anklecut_forward_2mm\';
cd([folder_baseline 'kinematics']);
kin_baseline_files=dir(['*.mot']);
cd([folder_achillescut 'kinematics']);
kin_achillescut_files=dir(['*.mot']);

for i=1:length(kin_baseline_files)
    data_baseline=importdata([kin_baseline_files(i).folder '\' kin_baseline_files(i).name]);
    data_achillescut=importdata([kin_achillescut_files(i).folder '\' kin_achillescut_files(i).name]);

    figure(i);
    set(gcf,'Position',[360,97.7,745,237.3]);
    subplot(1,3,1)
    plot(data_baseline.data(:,1),data_baseline.data(:,9));
    hold all;
    plot(data_achillescut.data(:,1),data_achillescut.data(:,9));
    title('hip flexion');
    ylabel('[º]');
    xlabel('time [s]');

    subplot(1,3,2)
    plot(data_baseline.data(:,1),data_baseline.data(:,12));
    hold all;
    plot(data_achillescut.data(:,1),data_achillescut.data(:,12));
    title('knee flexion');
    ylabel('[º]');
    xlabel('time [s]');

    subplot(1,3,3)
    plot(data_baseline.data(:,1),data_baseline.data(:,13));
    hold all;
    plot(data_achillescut.data(:,1),data_achillescut.data(:,13));
    title('ankle flexion');
    ylabel('[º]');
    xlabel('time [s]');
end

for i=1:length(S1.sol_val.forces_prescribed)
    figure(i);
    set(gcf,'Position',[360,97.7,745,237.3]);
    subplot(1,3,1)
    plot(S1.sol_val.tgrid{i}(S1.sol_val.t2plot{i}),S1.sol_val.forces_prescribed{i}(:,1));
    hold all;
    plot(S2.sol_val.tgrid{i}(S2.sol_val.t2plot{i}),S2.sol_val.forces_prescribed{i}(:,1));
    title('F_x');
    ylabel('[N]');
    xlabel('time [s]');

    subplot(1,3,2)
    plot(S1.sol_val.tgrid{i}(S1.sol_val.t2plot{i}),S1.sol_val.forces_prescribed{i}(:,2));
    hold all;
    plot(S2.sol_val.tgrid{i}(S2.sol_val.t2plot{i}),S2.sol_val.forces_prescribed{i}(:,2));
    title('F_y');
    ylabel('[N]');
    xlabel('time [s]');

    subplot(1,3,3)
    plot(S1.sol_val.tgrid{i}(S1.sol_val.t2plot{i}),S1.sol_val.forces_prescribed{i}(:,3));
    hold all;
    plot(S2.sol_val.tgrid{i}(S2.sol_val.t2plot{i}),S2.sol_val.forces_prescribed{i}(:,3));
    title('F_z');
    ylabel('[N]');
    xlabel('time [s]');
end
filepath=fileparts(matlab.desktop.editor.getActiveFilename)
cd(filepath);

%% Inverse Kinematics
clear all;
close all;
root_folder=pwd;
import org.opensim.modeling.*
generic_IKTool=InverseKinematicsTool('Setup_IK_generic.xml');

cd('kinematics');
path_Measurements=pwd;
folders=dir;
folders=folders(~ismember({folders.name},{'.','..'}));

for i=2:length(folders)
    cd([folders(i).folder '\' folders(i).name]);
    movS=dir('*.trc');
    
    for j=1:length(movS)
        data=importdata(movS(j).name,'\t',5);
        IKTool=generic_IKTool;
        IKTool.set_time_range(0,data.data(1,2));
        IKTool.set_time_range(1,data.data(end,2));
        IKTool.set_marker_file(movS(j).name);
        IKTool.setOutputMotionFileName(strrep(movS(j).name,'trc','mot'));
        IKTool.set_model_file(['..\..\rat_hindlimb_v44_kneepin_scaled.osim']);
        IKTool.setResultsDir(pwd);
        IKTool.run();
        fprintf(['Finished to process ' movS(j).name '\n']);
        
    end
end
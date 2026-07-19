

%% Inverse Kinematics
clear all;
close all;

filepath=fileparts(matlab.desktop.editor.getActiveFilename);
cd(filepath);

root_folder=pwd;
import org.opensim.modeling.*
generic_IKTool=InverseKinematicsTool('test_Setup_IK.xml');

in_folder=[filepath '\DataMay2025\achillescut_ForwardOnly\kinematics'];
out_folder=[filepath '\DataMay2025\achillescut_ForwardOnly_rel2Dangles\kinematics']; 
                                  
path_Measurements=pwd;

cd(in_folder);
movS=dir('*.trc');
    
for j=1:length(movS)
    cd(in_folder);
    data=importdata(movS(j).name,'\t',5);
    cd(out_folder);
    IKTool=generic_IKTool;
    IKTool.set_time_range(0,data.data(1,2));
    IKTool.set_time_range(1,data.data(end,2));
    IKTool.set_marker_file([in_folder '\' movS(j).name]);
    IKTool.setOutputMotionFileName([out_folder '\' strrep(movS(j).name,'trc','mot')]);
    IKTool.set_model_file([filepath '\rat_hindlimb_v44_kneepin_scaled_May2025_using2Dangles.osim']);
    IKTool.set_coordinate_file(strrep(movS(j).name,'.trc','_2Dangles.mot'));
    IKTool.setResultsDir(out_folder);
    IKTool.run();
    IKTool.print('test.xml');
    fprintf(['Finished to process ' movS(j).name '\n']);
    
end

% This function returns the indices of the muscles actuating the different 
% joints (for use with the moment arms) and the indices of the muscles in 
% the vector containing all muscle names.
%
% Author: Antoine Falisse, modified by Gil Serrancolí
% Date: 31/07/2022
% 
function [Indmusi,mai] = ...
    MomentArmIndices_3D(muscleNames,muscle_spanning_joint_INFO)

NMuscle = length(muscleNames);
for i = 1:length(muscleNames)
    % Muscles are ordered as left first and right second
    Indmusi.(muscleNames{i}(1:end)) = ...
        find(strcmp(muscleNames,muscleNames{i}));
    % We select the muscles that are actuating the joints
    % In muscle_spanning_joint_INFO: the columns are organized as follows:
    % 1) hip flex, 2) hip add, 3) hip rot, 4) knee flex, 5) ankle flex, 
    % 6) sub, 7) trunk ext, 8) trunk ben, 9) trunk rot
    for j = 1:size(muscle_spanning_joint_INFO,2)
        if (muscle_spanning_joint_INFO(i,j) == 1)
           mai(j).mus(1,i) = Indmusi.(muscleNames{i}(1:end));
        end        
    end    
end

for j = 1:size(muscle_spanning_joint_INFO,2)
    mai(j).mus(mai(j).mus == 0) = [];
end
end

%% Inverse Kinematics
clear all;
close all;

filepath=fileparts(matlab.desktop.editor.getActiveFilename);
cd(filepath);
prefix='anklecut_forward_2mm';
in_folder=[filepath '\DataAugust2025\' prefix '\kinematics'];
in_folder_nokneemarker=[filepath '\DataAugust2025\' prefix '\kinematics'];
in_folder_perturbation=[filepath '\DataAugust2025\' prefix '\perturbation'];

lfemur=38; %mm
ltibia=42; %mm

path_Measurements=pwd;

cd(in_folder);
movS=dir('*.trc');

for j=1:length(movS)
    trc_data=readtable([in_folder '\' movS(j).name],'FileType','text');
    cd(in_folder_nokneemarker);
    in_mot_data=importdata([strrep(movS(j).name,'.trc','.mot')]);
    q.labels=in_mot_data.colheaders;
    q.data(:,1)=in_mot_data.data(:,1);
    q.data(:,2:7)=repmat(in_mot_data.data(1,2:7),size(q.data,1),1); %sacrum_data
    q.data(:,8)=in_mot_data.data(:,8);
    pelvis_euler=q.data(:,2:4);

    %index position of markers in the matrix
    I_pelvis_top=find(contains(trc_data.Properties.VariableNames,'pelvis_top'));
    I_pelvis_bottom=find(contains(trc_data.Properties.VariableNames,'pelvis_bottom'));
    I_hip=find(contains(trc_data.Properties.VariableNames,'hip'));
    % I_knee=find(contains(trc_data.Properties.VariableNames,'knee'));
    I_ankle=find(contains(trc_data.Properties.VariableNames,'ankle'));
    I_mtp=find(contains(trc_data.Properties.VariableNames,'mtp'));

    nf=size(trc_data,1)-1;
    %position hip and ankle markers
    pos_pelvis_top=table2array(trc_data(2:end,I_pelvis_top:(I_pelvis_top+2)));
    pos_pelvis_bottom=table2array( trc_data(2:end,I_pelvis_bottom:(I_pelvis_bottom+2)));
    pos_hip=repmat(table2array(trc_data(2,I_hip:(I_hip+2))),nf,1);
    pos_ankle_0=table2array(trc_data(2:end,I_ankle:(I_ankle+2)));
    pos_mtp=table2array(trc_data(2:end,I_mtp:(I_mtp+2)));

    %create new ankle position
    perturbation=importdata([in_folder_perturbation '\' strrep(strrep(movS(j).name,'kinematics','motor'),'.trc','.csv')]);
    perturbation_name='filtered_position';
    Iperturbation=find(contains(perturbation.colheaders,perturbation_name));
    perturbation_in=perturbation.data(:,Iperturbation);
    p_centered=pos_ankle_0-mean(pos_ankle_0,1);
    [~,~,V]=svd(p_centered,'econ');
    v_ankle=V(:,1);
    if v_ankle(1)<0
        pos_ankle=pos_ankle_0(1,:)-(perturbation_in-perturbation_in(1))*v_ankle';
    else
        pos_ankle=pos_ankle_0(1,:)+(perturbation_in-perturbation_in(1))*v_ankle';
    end


    for i=1:nf
        [eulerZXY(i,:), pos_knee(i,:), Z_femur(i,:)] = computeHipEulerAnglesFromDistances(pos_hip(i,:), pos_pelvis_bottom(i,:), pos_pelvis_top(i,:), pos_ankle(i,:), pos_ankle(1,:), pos_ankle(end,:), lfemur, ltibia, pelvis_euler(i,:));

        knee_angle_deg(i,:) = -computeKneeAngle(pos_hip(i,:), pos_knee(i,:), pos_ankle(i,:));
        
        eulerZXY_ankle(i,:) = computeAnkleEulerAngles(pos_hip(i,:), pos_pelvis_bottom(i,:), pos_pelvis_top(i,:), pos_knee(i,:), pos_ankle(i,:), pos_ankle(1,:), pos_ankle(end,:), pos_mtp(i,:),Z_femur(i,:));
    end

    q.data(:,9:11)=eulerZXY;
    q.data(:,12)=knee_angle_deg;
    q.data(:,13:15)=eulerZXY_ankle;

    write_motionFile(q,[in_folder_nokneemarker '\' strrep(movS(j).name,'.trc','_2D.mot')])
    
    fprintf(['processed ' movS(j).name '\n']);
end


function pos_knee = ComputeKneePos(pos_hip, pos_ankle, pos_mtp, lfemur, ltibia)
% Finds points P4 at distance l1 from P1 and l2 from P2 in the plane defined by P1, P2, P3

    % Create orthonormal basis in the plane
    v1 = pos_ankle - pos_hip;
    v1 = v1 / norm(v1);  % x-axis

    normal = cross(pos_ankle - pos_hip, pos_mtp - pos_hip);
    normal = normal / norm(normal);

    v2 = cross(normal, v1);  % y-axis in the plane

    % Represent P4 in this 2D basis: P4 = P1 + x*v1 + y*v2
    % So define coordinates (x, y) in the plane
    
    % Setup the equations:
    % (x*v1 + y*v2) should be at distance l1 from origin (P1)
    % (x*v1 + y*v2 - (P2-P1)) should be at distance l2 from vector (P2-P1)

    % Let Q = x*v1 + y*v2 (in-plane vector from P1 to P4)
    % Then:
    % ||Q|| = l1
    % ||Q - (P2 - P1)|| = l2

    % Convert to a 2D problem in coordinates
    P2_plane = [dot(pos_ankle - pos_hip, v1,2), dot(pos_ankle - pos_hip, v2,2)];

    % Solve system:
    % (1) x^2 + y^2 = l1^2
    % (2) (x - x2)^2 + (y - y2)^2 = l2^2
    syms x y real
    eq1 = x^2 + y^2 == lfemur^2;
    eq2 = (x - P2_plane(1))^2 + (y - P2_plane(2))^2 == ltibia^2;
    sol = solve([eq1, eq2], [x, y], 'Real', true);

    % Back to 3D:
    P4_1 = pos_hip + double(sol.x(1)) * v1 + double(sol.y(1)) * v2;
    P4_2 = pos_hip + double(sol.x(2)) * v1 + double(sol.y(2)) * v2;
    
    if (P4_1(1)>pos_ankle(1))&&(P4_2(1)<=pos_ankle(1))
        pos_knee=P4_1;
    elseif (P4_1(1)<pos_ankle(1))&&(P4_2(1)>=pos_ankle(1))
        pos_knee=P4_2;
    else
        keyboard;
    end
end


function eul_deg_hip = computeHipEulerZXY(P1, P1a, P1b, P4, P2)
    % P1: Hip joint center
    % P1a, P1b: define pelvis sagittal plane
    % P4: Knee joint
    % P2: Ankle joint (defines femur movement plane)

   % Computes ZXY Euler angles of femur relative to pelvis

    %% === Pelvis Frame ===
    % X_pelvis: from hip to landmark (e.g., right ASIS)
    X_pelvis = P1b - P1;
    X_pelvis = X_pelvis / norm(X_pelvis);

    % Z_pelvis: normal to sagittal plane (P1, P1a, P1b)
    v1 = P1a - P1;
    v2 = P1b - P1;
    Z_pelvis = cross(v1, v2);
    Z_pelvis = Z_pelvis / norm(Z_pelvis);

    % Y_pelvis: completes right-handed system
    Y_pelvis = cross(Z_pelvis, X_pelvis);
    Y_pelvis = Y_pelvis / norm(Y_pelvis);

    % Rotation matrix of pelvis (columns are its axes)
    R_pelvis = [X_pelvis(:), Y_pelvis(:), Z_pelvis(:)];

    %% === Femur Frame ===
    % Y_femur: from knee to hip (longitudinal axis)
    Y_femur = P1 - P4;
    Y_femur = Y_femur / norm(Y_femur);

    % Z_femur: normal to hip-knee-ankle plane
    Z_femur = cross(P4 - P1, P2 - P1);
    Z_femur = Z_femur / norm(Z_femur);

    % Ensure Z_femur points in the same general direction as Z_pelvis
    if dot(Z_femur, Z_pelvis) < 0
        Z_femur = -Z_femur;
    end

    % X_femur: orthogonal axis
    X_femur = cross(Y_femur, Z_femur);
    X_femur = X_femur / norm(X_femur);

    % Re-orthogonalize (optional for numerical stability)
    Z_femur = cross(X_femur, Y_femur);
    Z_femur = Z_femur / norm(Z_femur);

    % Rotation matrix of femur
    R_femur = [X_femur(:), Y_femur(:), Z_femur(:)];

    %% === Relative Rotation ===
    R_relative = R_pelvis' * R_femur;

    %% === ZXY Euler angles ===
    eul_rad = rotm2eul(R_relative, 'ZXY');  % [Z, X, Y]
    eul_deg_hip = rad2deg(eul_rad);
end


function [eulerZXY, selectedP4, Z_femur] = computeHipEulerAnglesFromDistances(P1, P1a, P1b, P2, P2a, P2b, l1, l2, pelvis_euler)
    % Step 1: Define the plane from P2a, P2b, P1
    plane_origin = P2a;
    v1 = P2b - P2a;
    v2 = P1 - P2a;
    plane_normal = cross(v1, v2);
    plane_normal = plane_normal / norm(plane_normal);

    % Step 2: Orthonormal basis for the plane
    ex = v1 / norm(v1);
    ez = plane_normal;
    ey = cross(ez, ex);
    
    % Transformation matrix to plane basis
    R_plane = [ex(:), ey(:), ez(:)];
    T_plane = @(p) R_plane' * (p(:) - plane_origin(:));
    T_world = @(q) (R_plane * q(:)) + plane_origin(:);
    
    % Step 3: Project points into the plane frame
    P1p = T_plane(P1);  % Hip in plane frame
    P2p = T_plane(P2);  % Ankle in plane frame

    % Step 4: Intersect 2 circles in 2D (on plane)
    [P4a_2d, P4b_2d] = intersectCircles2D(P1p(1:2), l1, P2p(1:2), l2);

    if isempty(P4a_2d)
        error('No valid P4: spheres do not intersect the plane with consistent distances.');
    end

    % Step 5: Lift back to 3D
    P4a = T_world([P4a_2d 0]);
    P4b = T_world([P4b_2d 0]);
    
    if P4a(1) > P2(1)
        selectedP4 = P4a;
    elseif P4b(1) > P2(1)
        selectedP4 = P4b;
    else
        warning('Neither P4 candidate has x > P2.x. Defaulting to P4a.');
        selectedP4 = P4a;
    end

    % Choose one candidate — e.g., the one with Z_femur aligned positively with Z_pelvis
    [eulerZXY, ok, Z_femur] = tryEuler(P1, P1a, P1b, P2, P2a, P2b, selectedP4, pelvis_euler);

end

function [thetaZXY, valid, Z_femur] = tryEuler(P1, P1a, P1b, P2, P2a, P2b, P4, pelvis_euler)
    
    P4=P4(:)';
    valid = true;

    %wrong way, not accurate, better to take the orientation from IK
    % X_pelvis = (P1b - P1) / norm(P1b - P1);
    % Z_pelvis = cross(P1b - P1a, P1 - P1a); Z_pelvis = Z_pelvis / norm(Z_pelvis);
    % Y_pelvis = cross(Z_pelvis, X_pelvis); Y_pelvis = Y_pelvis / norm(Y_pelvis);
    % R_pelvis = [X_pelvis(:), Y_pelvis(:), Z_pelvis(:)];
    R_pelvis=eul2rotm(pelvis_euler*pi/180,'ZXY');
    Z_pelvis=R_pelvis(:,3);

    Y_femur = (P1 - P4) / norm(P1 - P4);
    Z_femur = cross(P2a - P1, P2b - P1);
    Z_femur = Z_femur / norm(Z_femur);
    % Ensure consistent direction with pelvis Z axis
    if dot(Z_femur, Z_pelvis) < 0
        Z_femur = -Z_femur;
    end
    % Compute orthogonal X axis
    X_femur = cross(Y_femur, Z_femur);
    X_femur = X_femur / norm(X_femur);
    % Re-orthogonalize Z in case of numerical error
    Z_femur = cross(X_femur, Y_femur);
    Z_femur = Z_femur / norm(Z_femur);

    R_femur = [X_femur(:), Y_femur(:), Z_femur(:)];

    R = R_pelvis' * R_femur;
    thetaZXY=rotm2eul(R,'ZXY')*180/pi;
end

function [p1, p2] = intersectCircles2D(c1, r1, c2, r2)
    c1 = c1(:)';  % Convert to row
    c2 = c2(:)';
    % Compute intersection of two circles in 2D
    d = norm(c2 - c1);
    if d > r1 + r2 || d < abs(r1 - r2)
        p1 = []; p2 = []; return;
    end
    a = (r1^2 - r2^2 + d^2) / (2*d);
    h = sqrt(r1^2 - a^2);
    p_mid = c1 + a*(c2 - c1)/d;
    
    dir = (c2 - c1) / d;  % unit vector from c1 to c2
    perp = [-dir(2), dir(1)];

    p1 = p_mid + h*perp;
    p2 = p_mid - h*perp;
end

function knee_angle_deg = computeKneeAngle(P1, P4, P2)
    % Inputs:
    % P1 - Hip position [1x3]
    % P4 - Knee position [1x3]
    % P2 - Ankle position [1x3]

    % Step 1: Define Y axes
    Y_femur = (P1 - P4) / norm(P1 - P4);  % from knee to hip
    Y_tibia = (P4 - P2) / norm(P4 - P2);  % from ankle to knee

    % Step 2: Compute angle between Y axes
    cos_theta = dot(Y_femur, Y_tibia);  % dot product
    cos_theta = min(max(cos_theta, -1), 1);  % clamp for numerical safety
    angle_rad = acos(cos_theta);  % angle in radians

    % Step 3: Convert to degrees
    knee_angle_deg = rad2deg(angle_rad);
end


function eulerZXY = computeAnkleEulerAngles(P1, P1a, P1b, P4, P2, P2a, P2b, P3, Z_femur)
    % Inputs:
    % - P1a, P1b: define plane with P2 and P4
    % - P4: toe (world)
    % - P2: ankle (world)
    % - P3_local: MTP in foot local frame
    % Z_femur: global femur Z axis (shared with tibia and foot)

    %% 1. Tibia frame
    Y_tibia = (P4 - P2) / norm(P4 - P2);
    Z_tibia = Z_femur / norm(Z_femur);
    X_tibia = cross(Y_tibia, Z_tibia);
    X_tibia = X_tibia / norm(X_tibia);
    R_tibia = [X_tibia(:), Y_tibia(:), Z_tibia(:)];

    %% 2. Foot frame
    % Z axis of foot is parallel to Z_femur
    Z_foot = Z_femur / norm(Z_femur);

    % Plane defined by P2a, P2b, P1
    v1 = P2b - P2a;
    v2 = P1 - P2a;
    plane_normal = cross(v1, v2);  % orthogonal to plane
    plane_normal = plane_normal / norm(plane_normal);
    % Vector from point on plane to P3
    vec = P3 - P1a;
    % Distance from P3 to plane along the normal
    d = dot(vec, plane_normal);
    % Project P3 onto plane by subtracting normal component
    P3_proj = P3 - d * plane_normal;

    %Compute X_foot, Y_foot and Z_foot
    X_foot_tilde = P3_proj - P2;
    X_foot_tilde = X_foot_tilde / norm(X_foot_tilde);
    theta = deg2rad(30);  % rotation angle
    Rz = axang2rotm([Z_foot(:)' theta]);  % rotation matrix around Z_foot
    X_foot = Rz * X_foot_tilde(:);  % rotated vector
     Y_foot = cross(Z_foot, X_foot);
    Y_foot = Y_foot / norm(Y_foot);
    X_foot = cross(Y_foot, Z_foot);
    X_foot = X_foot / norm(X_foot);

    R_foot = [X_foot(:), Y_foot(:), Z_foot(:)];

    %% 3. Relative rotation and ZXY Euler angles
    R_relative = R_tibia' * R_foot;
    eulerZXY = rad2deg(rotm2eul(R_relative, 'ZXY'));  % MATLAB intrinsic ZXY
end

function R = R_from_axes(Z_axis, P1a, P1b, P2, P4)
    % Construct X in the plane defined by P1a, P1b, P2, P4
    v1 = P1b - P1a;
    v2 = P4 - P1a;
    plane_normal = cross(v1, v2);
    plane_normal = plane_normal / norm(plane_normal);

    Z = Z_axis / norm(Z_axis);
    X = cross(plane_normal, Z); X = X / norm(X);
    Y = cross(Z, X);

    R = [X(:), Y(:), Z(:)];
end
function eul = rotm2eulZXY_fast(R)
%ROTM2EULZXY_FAST Convert rotation matrices to intrinsic ZXY Euler angles.
%
%   eul = rotm2eulZXY_fast(R)
%
% Input
%   R   3x3 or 3x3xN rotation matrix array
%
% Output
%   eul Nx3 array containing [zAngle, xAngle, yAngle], in radians.
%
% This implements the principal ZXY solution, with:
%   x in [-pi/2, pi/2]
%   z and y in [-pi, pi]
%
% R must use MATLAB's premultiplication convention.

    arguments
        R (3,3,:) double
    end

    % Protect asin from small floating-point violations such as 1+eps.
    sx = reshape(R(3,2,:), [], 1);
    sx = max(-1, min(1, sx));

    x = asin(sx);
    cx = sqrt(max(0, 1 - sx.^2));

    r12 = reshape(R(1,2,:), [], 1);
    r22 = reshape(R(2,2,:), [], 1);
    r31 = reshape(R(3,1,:), [], 1);
    r33 = reshape(R(3,3,:), [], 1);

    z = atan2(-r12, r22);
    y = atan2(-r31, r33);

    % Gimbal-lock handling: cos(x) approximately zero.
    singular = cx < 1e-10;

    if any(singular)
        r11 = reshape(R(1,1,:), [], 1);
        r21 = reshape(R(2,1,:), [], 1);

        % At the singularity, z and y are not independently identifiable.
        % Choose y = 0 and absorb the combined rotation into z.
        z(singular) = atan2(r21(singular), r11(singular));
        y(singular) = 0;
    end

    eul = [z, x, y];
end
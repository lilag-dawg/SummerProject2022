clc
clear
close all

imu_data = readtable('magneto_accel_measurements.csv');
ref_quat_data = readtable('quaternionRef.csv');

ref_quat = [ref_quat_data{:,"w"}, ref_quat_data{:,"x"}, ref_quat_data{:,"y"}, ref_quat_data{:,"z"}];

magnetometer = [imu_data{:,"magnetX"}, imu_data{:,"magnetY"}, imu_data{:,"magnetZ"}, imu_data{:,"Time"}];
accelerometer = [imu_data{:,"accelX"}, imu_data{:,"accelY"}, imu_data{:,"accelZ"}, imu_data{:,"Time"}];

% Loop
n = height(imu_data);

% grahpic stuff
eulerEstimate = zeros(n,3);
eulerRef = zeros(n,3);

for i = 1:n
    dcm = ecompass(magnetometer(i,1:3),accelerometer(i,1:3))
    q = mat2quat(dcm);

    % graphic stuff
    [alpha,beta,gamma] = quat2euler(q);
    % XYZ
    eulerEstimate(i,:) = 180/pi*[alpha, beta, gamma];
    
    [alpha,beta,gamma] = quat2euler(ref_quat(i,:));
    eulerRef(i,:) = 180/pi*[alpha, beta, gamma];
end

tiledlayout(3,1);

ax1 = nexttile;
plot(ref_quat_data{:,"Time"},eulerRef(:,1),imu_data{:,"Time"},eulerEstimate(:,1))
title(ax1,'Angle \phi')
legend('Real','Estimated')

ax2 = nexttile;
plot(ref_quat_data{:,"Time"},eulerRef(:,2),imu_data{:,"Time"},eulerEstimate(:,2))
title(ax2,'Angle \theta')
legend('Real','Estimated')

ax3 = nexttile;
plot(ref_quat_data{:,"Time"},eulerRef(:,3),imu_data{:,"Time"},eulerEstimate(:,3))
title(ax3,'Angle \psi')
legend('Real','Estimated')


function dcm = ecompass(mag,acc)
    mag_unit = mag./norm(mag);
    acc_unit = acc./norm(acc);

    D = -acc_unit;
    N = mag_unit;
    E = cross(D,N);
    
    dcm = [N;E;D];
    dcm = dcm./norm(dcm);
end

function [alpha,beta,gamma] = quat2euler(q)
    R = quat2mat(q);

    % convert using XYZ convention
    
    alpha = atan2(-R(2,3),R(3,3));

    if (1-R(1,3)^2 < 0)
         beta = atan2(R(1,3),0);
    else
         beta = atan2(R(1,3),sqrt(1-R(1,3)^2));
    end
    gamma = atan2(-R(1,2),R(1,1));
end

function q = mat2quat(R)
    q = zeros(4,1);

    tr = R(1,1) + R(2,2) + R(3,3);
    
    if (tr > 0)
        S = sqrt(tr+1)*2;
        q(1) = S/4;
        q(2) = (R(3,2)-R(2,3))/S;
        q(3) = (R(1,3)-R(3,1))/S;
        q(4) = (R(2,1)-R(1,2))/S;
    elseif(R(1,1)>R(2,2) && R(1,1)>R(3,3))
        S = sqrt(1+R(1,1)-R(2,2)-R(3,3))*2;
        q(1) = (R(3,2)-R(2,3))/S;
        q(2) = S/4;
        q(3) = (R(1,2)+R(2,1))/S;
        q(4) = (R(1,3)+R(3,1))/S; 
    elseif(R(2,2)>R(3,3))
        S = sqrt(1+R(2,2)-R(1,1)-R(3,3))*2;
        q(1) = (R(1,3)-R(3,1))/S;
        q(2) = (R(1,2)+R(2,1))/S;
        q(3) = S/4;
        q(4) = (R(2,3)+R(3,2))/S;         
    else
        S = sqrt(1+R(3,3)-R(1,1)-R(2,2))*2;
        q(1) = (R(2,1)-R(1,2))/S;
        q(2) = (R(1,3)+R(3,1))/S;
        q(3) = (R(2,3)+R(3,2))/S; 
        q(4) = S/4;
    end
end


function R = quat2mat(q)
    R = zeros(3,3);

    R(1,1) = q(1)^2+q(2)^2-q(3)^2-q(4)^2;
    R(1,2) = 2*(q(2)*q(3)-q(1)*q(4));
    R(1,3) = 2*(q(1)*q(3)+q(2)*q(4));
    R(2,1) = 2*(q(2)*q(3)+q(1)*q(4));
    R(2,2) = q(1)^2-q(2)^2+q(3)^2-q(4)^2;
    R(2,3) = 2*(q(3)*q(4)-q(1)*q(2));
    R(3,1) = 2*(q(2)*q(4)-q(1)*q(3));
    R(3,2) = 2*(q(1)*q(2)+q(3)*q(4));
    R(3,3) = q(1)^2-q(2)^2-q(3)^2+q(4)^2;
end
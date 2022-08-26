clc
clear
close all

imu_data = readtable('gyroData.csv');
ref_quat_data = readtable('quaternionRef.csv');

gyroscope = [imu_data{:,"gyroX"}, imu_data{:,"gyroY"}, imu_data{:,"gyroZ"}, imu_data{:,"Time"}];
ref_quat = [ref_quat_data{:,"w"}, ref_quat_data{:,"x"}, ref_quat_data{:,"y"}, ref_quat_data{:,"z"}];

n = height(imu_data);
q = [1 0 0 0]';

% grahpic stuff
eulerEstimate = zeros(n,3);
eulerRef = zeros(n,3);


for i = 1:n
    if i-1 ~= 0
        h = gyroscope(i,4) - gyroscope(i-1,4);
        q = poissonIntegration(gyroscope(i,1:3),q,h);
    else
        h = 0;
        q = poissonIntegration(gyroscope(i,1:3),q,h);
    end

    % graphic stuff
    [phi,teta,psi] = quat2euler(q);
    % XYZ
    eulerEstimate(i,:) = 180/pi*[phi,teta,psi];
    
    [phi,teta,psi] = quat2euler(ref_quat(i,:));
    eulerRef(i,:) = 180/pi*[phi,teta,psi];
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




function q_k = poissonIntegration(w,q_k_1,h)

    w = -w;
    Q = [0 -w(1) -w(2) -w(3);
         w(1) 0 w(3) -w(2);
         w(2) -w(3) 0 w(1);
         w(3) w(2) -w(1) 0];
    
    if norm(w) ~=0
        q_k = (cos(h/2*norm(w))*eye(4) + sin(h/2*norm(w))/norm(w)*Q)*q_k_1;
    else
        q_k = q_k_1;
    end

    q_k = q_k./norm(q_k);
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
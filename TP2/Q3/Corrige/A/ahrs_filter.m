clc
clear
close all

imu_data = readtable('imu_data_a.csv');
ref_quat_data = readtable('quaternionRef.csv');

ref_quat = [ref_quat_data{:,"w"}, ref_quat_data{:,"x"}, ref_quat_data{:,"y"}, ref_quat_data{:,"z"}];

magnetometer = [imu_data{:,"magnetX"}, imu_data{:,"magnetY"}, imu_data{:,"magnetZ"}, imu_data{:,"Time"}];
accelerometer = [imu_data{:,"accelX"}, imu_data{:,"accelY"}, imu_data{:,"accelZ"}, imu_data{:,"Time"}];
gyroscope = [imu_data{:,"gyroX"}, imu_data{:,"gyroY"}, imu_data{:,"gyroZ"}, imu_data{:,"Time"}];

% EKF initialisation

    % constante
    Ha = eye(3);
    Hw = eye(3);
    Fw = zeros(3,3);
    Fa = zeros(3,3);
    g = [0 0 9.81]';
    m = [1 0 0]';

    % noise Parameter
    gyroVar = 1.45e-5; 
    gyroAddVar = 7.615e-4;
    NEDVar = 1e-5;
    accelVar = 1e-5;
    accelAddVar = 9e-5;
    magnetoVar = 5.3e-5;
    
    Q = [NEDVar*eye(3) zeros(3,3) zeros(3,3) zeros(3,3);
         zeros(3,3) gyroVar*eye(3) zeros(3,3) zeros(3,3);
         zeros(3,3) zeros(3,3) gyroAddVar*eye(3) zeros(3,3);
         zeros(3,3) zeros(3,3) zeros(3,3) accelAddVar*eye(3)];

    x_k = zeros(9,1);
    sigma_k = eye(9);

% Loop
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
        dcm = ecompass(magnetometer(i,1:3),accelerometer(i,1:3));
        q = mat2quat(dcm);
        q = poissonIntegration(gyroscope(i,1:3),q,h);
    end

    % EKF

    R_bn = quat2mat(q);

    R_nb = R_bn';

    A = [zeros(3,3) -R_bn*Hw zeros(3,3);
         zeros(3,3) Fw zeros(3,3);
         zeros(3,3) zeros(3,3) Fa];
    B = [-eye(3) -R_bn zeros(3,3) zeros(3,3);
         zeros(3,3) zeros(3,3) eye(3) zeros(3,3);
         zeros(3,3) zeros(3,3) zeros(3,3) eye(3)];

    if ~isnan(accelerometer(i,:))
        % predict
        x_k_k_1 = A*x_k;
        sig_k_k_1 = A*sigma_k*A' + B*Q*B';

        % update for accelerometer measurement
        skewG = [0 g(3) 0;
                -g(3) 0 0;
                 0 0 0];
        rotHa = R_bn*Ha;

        C = [skewG zeros(3,3) rotHa];
        D = R_bn;

        Pyy = C*sig_k_k_1*C' + D*(accelVar*eye(3))*D';
        Pxy = sig_k_k_1*C';
        Pyx = C*sig_k_k_1;
        sigma_k = sig_k_k_1 - (Pxy/Pyy)*Pyx;

        y_k = g - R_bn*(Ha*x_k(7:9)-accelerometer(i,1:3)');
        
        x_k = Pxy/Pyy * (y_k - C*x_k_k_1);
        %todo change see chapter 10 p.6
    end

    if ~isnan(magnetometer(i,:))
        % predict
        x_k_k_1 = A*x_k;
        sig_k_k_1 = A*sigma_k*A' + B*Q*B';

        % update for magnetometer measurement
        skewM = [0 0 0;
                 0 0 m(1);
                 0 -m(1) 0];
    
        C = [skewM zeros(3,3) zeros(3,3)];
        D = -R_bn;
    
        Pyy = C*sig_k_k_1*C' + D*(magnetoVar*eye(3))*D';
        Pxy = sig_k_k_1*C';
        Pyx = C*sig_k_k_1;
        sigma_k = sig_k_k_1 - (Pxy/Pyy)*Pyx;
        
        y_k = m - R_bn*magnetometer(i,1:3)';

        x_k = Pxy/Pyy * (y_k - C*x_k_k_1);
    end

    if any(isnan(accelerometer(i,:))) && any(isnan(magnetometer(i,:)))
        x_k = A*x_k;
    end

    M = rodrigue(x_k(1:3));

    R = R_nb*M;
    q = mat2quat(R');

    x_k(1:3) = 0;
    
    % graphic stuff
    [phi,teta,psi] = quat2euler(q);
    % XYZ
    eulerEstimate(i,:) = 180/pi*[phi,teta,psi];
end

% Graphic stuff
for i = 1:n
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




function dcm = ecompass(mag,acc)
    mag_unit = mag./norm(mag);
    acc_unit = acc./norm(acc);

    D = -acc_unit;
    N = mag_unit;
    E = cross(D,N);
    E = E./norm(E);
    
    dcm = [N;E;D];
    dcm = dcm./norm(dcm);
end

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

function M = rodrigue(p)
    skewP = [0 -p(3) p(2);
             p(3) 0 -p(1);
            -p(2) p(1) 0];
    
    normP = norm(p);
    M = eye(3) - skewP/normP*sin(normP) +  skewP^2/normP^2*(1-cos(normP));
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

% function q = mat2quat(R)
%     q = zeros(4,1);
% 
%     q(1) = sqrt(1 + R(1,1) + R(2,2) + R(3,3));
%     q(2) = sign(R(3,2)-R(2,3))*sqrt(1 + R(1,1) - R(2,2) - R(3,3));
%     q(3) = sign(R(1,3)-R(3,1))*sqrt(1 - R(1,1) + R(2,2) - R(3,3));
%     q(4) = sign(R(2,1)-R(1,2))*sqrt(1 - R(1,1) - R(2,2) + R(3,3));
% 
%     q = 1/2*q;
% end
% 
% function y = sign(x)
%     if(x >= 0)
%         y = 1;
%     else
%         y = -1;
%     end
% end

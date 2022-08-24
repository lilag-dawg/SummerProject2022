clc
clear
close all

odometryPosition = readtable("estimated_odometry_position.csv");
odometryData=readtable('odometry_data.csv');
positionFixingData = readtable('laser_data.csv');
truePosition = readtable('true_pose.csv');

% constantes
wheel_radius = 0.0975;
robot_width = 0.1655;
sampleTime = odometryData{2,3} - odometryData{1,3};

% bruit des donnees
varianceOdometry = 0.01;
variancePositionFixing = 0.05;


% intialisation du filtre
Q = eye(2)*varianceOdometry;
R = eye(2)*variancePositionFixing;
C = [1 0 0; 0 1 0];

x_k = [0,0,0]';
sigma_k = eye(3);

n = height(odometryData);

n_hat = 0;
e_hat = 0;
psi_hat = 0;
time = 0;

% EKF

for k = 1:n
    u_hat = wheel_radius/2 * (odometryData{k,"u_L"}+odometryData{k,"u_R"});
    
    n_hat = odometryPosition{k,'X'};
    e_hat = odometryPosition{k,'Y'};
    psi_hat = odometryPosition{k,'Psi'};

    A = [0 0 -u_hat*sin(psi_hat);
         0 0 u_hat*cos(psi_hat);
         0 0 0];

    B = [1/2*cos(psi_hat) 1/2*cos(psi_hat);
         1/2*sin(psi_hat) 1/2*sin(psi_hat);
         -1/(2*robot_width) 1/(2*robot_width)]; % [wr, wl]


    % predict
    x_k_k_1 = A*x_k;
    sig_k_k_1 = A*sigma_k*A' + B*Q*B';

    % update

    Pyy = C*sig_k_k_1*C' + R;
    Pxy = sig_k_k_1*C';
    Pyx = C*sig_k_k_1;
    sig_k = sig_k_k_1 - (Pxy/Pyy)*Pyx;

    if ~isnan(positionFixingData{k,'X'})

        y_k = [positionFixingData{k,'X'} - n_hat;
               positionFixingData{k,'Y'} - e_hat];
    
        x_k = x_k_k_1 + Pxy/Pyy * (y_k - C*x_k_k_1);
    end

    odometryPosition{k,"X"} = n_hat + x_k(1);
    odometryPosition{k,"Y"} = e_hat + x_k(2);
    odometryPosition{k,"Psi"} = psi_hat + x_k(3); % x_k(3) should always be 0.
end

% erreurs moyenne
ex = zeros(1,n);
ey = zeros(1,n);
e_psi = zeros(1,n);

for i = 1:n
    ex(i) = odometryPosition{i,'X'} - truePosition{i,'X'};
    ey(i) = odometryPosition{i,'Y'}  - truePosition{i,'Y'};
    e_psi(i) = odometryPosition{i,'Psi'}  - truePosition{i,'Psi'};
    
end

ex_moy = sum(ex)/n
ey_moy = sum(ey)/n
e_psi_moy = sum(e_psi)/n

figure
plot(odometryPosition{:,"X"},odometryPosition{:,"Y"},truePosition{:,'X'},truePosition{:,'Y'})
legend('calculé','réelle')
title('Position du robot')
xlabel('x (m)')
ylabel('y (m)')


figure
plot(positionFixingData{:,'X'},positionFixingData{:,'Y'},truePosition{:,'X'},truePosition{:,'Y'})
legend('laser','réelle')


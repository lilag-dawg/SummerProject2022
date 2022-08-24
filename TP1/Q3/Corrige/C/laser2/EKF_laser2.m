clc
clear
close all

odometryPosition = readtable("estimated_odometry_position.csv");
odometryData=readtable('odometry_data.csv');
lookupTable  = readtable('lookupTable.csv');
positionFixingData = readtable('laser_data.csv');
truePosition = readtable('true_pose.csv');


% constantes
wheel_radius = 0.0975;
robot_width = 0.1655;
sampletime = 0.05;

% bruit des donnees
varianceOdometry = 0.01;
variancePositionFixing = 0.05;

% intialisation du filtre
Q = eye(2)*varianceOdometry;
R = variancePositionFixing;

x_k = [0,0,0]';
sigma_k = eye(3);

n = height(odometryData);

n_hat = 0;
e_hat = 0;
psi_hat = 0;
time = 0;

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
    lia = ismembertol(positionFixingData{:,"Time"},time,1e-5);
    rows = find(lia == 1);
    m = length(rows);
    
    for i = 1:m
        [C,D] = getMeasurementMatrices(psi_hat,positionFixingData,rows(i));
    
        Pyy = C*sig_k_k_1*C' + D*R*D';
        Pxy = sig_k_k_1*C';
        Pyx = C*sig_k_k_1;
        sigma_k = sig_k_k_1 - (Pxy/Pyy)*Pyx;
    
        y_k = getMeasurement(n_hat,e_hat,psi_hat,positionFixingData,lookupTable,rows(i));
    
        x_k = x_k_k_1 + Pxy/Pyy * (y_k - C*x_k_k_1);

        n_hat = n_hat + x_k(1);
        e_hat = e_hat + x_k(2);
        psi_hat = psi_hat + x_k(3);
    end

    odometryPosition{k,"X"} = n_hat;
    odometryPosition{k,"Y"} = e_hat;
    odometryPosition{k,"Psi"} = psi_hat;

    time  = time + sampletime;
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



function [C,D] = getMeasurementMatrices(psi_hat,data,index)

    d = data{index,"Distance"};
    teta = data{index,"Teta"};

    xi = d*cos(teta);
    yi = d*sin(teta);

    C = [1 0 -xi*sin(psi_hat)-yi*cos(psi_hat);
         0 1 xi*cos(psi_hat)-yi*sin(psi_hat)];

    D = [cos(teta)*cos(psi_hat)-sin(teta)*sin(psi_hat);
         cos(teta)*sin(psi_hat)+sin(teta)*cos(psi_hat)];
end

function y = getMeasurement(n_hat,e_hat,psi_hat,data,tab,index)

    d = data{index,"Distance"};
    teta = data{index,"Teta"};

    R = [cos(psi_hat) -sin(psi_hat);
         sin(psi_hat) cos(psi_hat)];

    lia = ismember(tab{:,"Key"},data{index,"Key"});
    row = find(lia == 1);
    
    P = [tab{row,"X"} tab{row,"Y"}]';
    P_hat = R*[d*cos(teta);d*sin(teta)] + [n_hat;e_hat];

    y = P-P_hat;
end


clc
clear
close all

odometryData=readtable('odometry_data.csv');
truePosition = readtable('true_pose.csv');

% constantes
wheel_radius = 0.0975;
robot_width = 0.1655;
sampleTime = odometryData{2,3} - odometryData{1,3};

n = height(odometryData);

position = zeros(2,n);
position(:,1) = [-1.6312 0]';

psi = zeros(1,n);
psi(1) = 0;

for k = 1:n-1
    u_hat = 1/2 * (odometryData{k,"u_R"}+odometryData{k,"u_L"});
    w_hat = 1/(2*robot_width)*(odometryData{k,"u_L"}-odometryData{k,"u_R"}); % uL-uR
    
    position(1,k+1) = position(1,k) + u_hat*cos(psi(k))*sampleTime;
    position(2,k+1) = position(2,k) + u_hat*sin(psi(k))*sampleTime;

    psi(k+1) = psi(k) + w_hat*sampleTime;
end

M(:,1:2) = position';
M(:,3) = psi';
M(:,4) = odometryData{:,3};

T = array2table(M);
T.Properties.VariableNames(1:4) = {'X','Y','Psi','Time'};
writetable(T,'estimated_odometry_position.csv');

figure
plot(T{:,'X'},T{:,'Y'},truePosition{:,'X'},truePosition{:,'Y'})
legend('laser','réelle')

% erreurs moyenne
ex = zeros(1,n);
ey = zeros(1,n);
e_psi = zeros(1,n);

for i = 1:n
    ex(i) = T{i,'X'} - truePosition{i,'X'};
    ey(i) = T{i,'Y'}  - truePosition{i,'Y'};
    e_psi(i) = T{i,'Psi'}  - truePosition{i,'Psi'};
    
end

ex_moy = sum(ex)/n
ey_moy = sum(ey)/n
e_psi_moy = sum(e_psi)/n




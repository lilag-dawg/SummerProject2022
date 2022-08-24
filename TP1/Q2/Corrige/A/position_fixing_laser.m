clc
clear
close all


T=readtable('laser_data.csv');

position = [];

time = 0;
x = [-1,-1]';
p = [];
s = [];
i = 1;
while i < height(T)
    for k = i:height(T)
        if T{k,4} <= time
            p = [p;T{k,3}];
            s = [s;[T{k,1},T{k,2}]];
        else
            time = T{k,4};
            break
        end
    end

    if length(p) >=2
        J = findJ(x,s,p);
        g = findG(x,s,p);
        h = -J\g;
        
        while norm(h) >= 1e-8
            J = findJ(x,s,p);
            g = findG(x,s,p);
            h = -J\g;
            x = x+h;
        end
    end
    position = [position;x'];

    p = [];
    s = [];
    i = k;
end

figure
plot(position(:,1),position(:,2))

title('Position du robot')
xlabel('x (m)')
ylabel('y (m)')


function J = findJ(x,s,p)
    J = [];
    for i = 1:length(p)
        J(i,:) = ((x-s(i,:)')/norm(x-s(i,:)'))';
    end
end

function g = findG(x,s,p)
    g = [];
    for i = 1:length(p)
        g(i,:) = norm(x-s(i,:)')-p(i);
    end
end
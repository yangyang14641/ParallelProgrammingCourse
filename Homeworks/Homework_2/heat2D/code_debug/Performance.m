%-----------------------------------------
%%%% This program is designed for parallel code performance plot...
%%%% YANG YANG
%%%% Peking Univerisity
%%%% Jan 5 2017
%------------------------------------------
%% data
numThread = 6;   % number of computing threads
N = [120^2, 600^2, 1200^2, 2400^2];  % Problem scale
speedUp=[0.716638, 2.132516, 1.773081, 1.696943];
parallelCost = [12.2399506860, 110.7361543450, ...
    549.4440846960, 2490.0252091020];
serialCost = [8.7716123170, 236.1466142150, ...
    974.2090435070, 4225.4302848400];


%% data plot
% speed up plot
figure('Color',[1 1 1]);
plot(N,speedUp,'ro-','linewidth',2,'markersize',20);
xlabel('Probelm Scale');
ylabel('Speed Up');
title('Speed Up figure');

% time cost plot
figure('Color',[1 1 1]);
plot(N,parallelCost,'ro-','linewidth',2,'markersize',20);
hold on;
plot(N,serialCost,'bx-','linewidth',2,'markersize',20);
hold off;
xlabel('Probelm Scale');
ylabel('Time Cost (s)');
title('Time Cost figure');
legend('Parallel','Serial');

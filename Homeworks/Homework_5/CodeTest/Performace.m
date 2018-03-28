%-------------------------------------------------------------------------
%%%% Progrm is designed for plot time cost in code performance test...
%%%% YANG YANG
%%%% Peking University
%%%% Jan 5 2017
%------------------------------------------
%% data
% problem scale
N = [120^2, 600^2, 1200^2];           
% time cost
serial_cost = [14.198089, 357.691399, 1428.974339]; 
pthread_cost = [12.167301, 128.956751, 466.923103];
mpi_cost = [1.9713819026947021E+01, 4.8343813109397888E+02, 1.9419226799011230E+03];

pthread_speedup = serial_cost ./ pthread_cost;
mpi_speedup = serial_cost ./ mpi_cost;

%% plot data
% speed up plot
figure('Color',[1 1 1]);
plot(N,pthread_speedup,'ro-','linewidth',2,'markersize',10);
hold on;
plot(N,mpi_speedup,'bx-','linewidth',2,'markersize',10);
hold off;
xlabel('Problem scale');
ylabel('Speed Up');
legend('Pthread','MPI');
title('Speed Up Figure');

% time cost plot
figure('Color',[1 1 1]);
plot(N,serial_cost,'m*-','linewidth',2,'markersize',10);
hold on;
plot(N,pthread_cost,'ro-','linewidth',2,'markersize',10);
plot(N,mpi_cost,'b.-','linewidth',2,'markersize',10);
hold off;
xlabel('Problem Scale');
ylabel('Time Cost (sec)');
title('Time Cost Figure');
legend('Serial','Pthread','MPI');

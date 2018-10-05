% copyright Luca Daniel, MIT 2018
% modified by Zijing 2018/9/22
clear;
close all;

%select input evaluation functions
eval_u = 'eval_u_step';
Nx=15;
Ny=15;
% 2D Heat Process
tmsize = 4; %Dimension of tumor square
[p,x_start,t_start,t_stop,max_dt_FE] = getParam_2DHeating(1,1,Nx,Ny,tmsize);
eval_f = 'eval_f_LinearSystem';

p.Nx=Nx;
p.Ny=Ny;
% test FE function
timestep=0.02;
t_stop=100;

figure(1)
visualize=1;
[X] = ForwardEuler2D(eval_f,x_start,eval_u,p,t_start,t_stop,timestep,visualize);
%% Jacobian Calculation
Jf_analytic = eval_Jf_LinearSystem(x_start,1,p);
dx=0.001;
N=Nx*Ny;
Jf_Dif=zeros(N);
for i=1:N
    x_new=x_start;
    x_new(i)=x_new(i)+dx;
    temp=(eval_f_LinearSystem(x_new,1,p)-eval_f_LinearSystem(x_start,1,p))./dx;
    Jf_Dif(:,i)=temp;
end
error=Jf_analytic-Jf_Dif;
disp(mean(error(:)));



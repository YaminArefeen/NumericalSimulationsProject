%In this script, Yamin explores the difference between Forward
%Euler and Trapezoidal integration for our project. 

%% 0) Initializions
% clear;
close all;
addpath funcs;
Nx=128;   % larger matrix size shows higher resolution
Ny=128;
load('Image_mask_128_128.mat')

figure('units','normalized','outerposition',[0 0 1 1]); 
subplot(1,4,1); imshow(T1_image,[]); title('MR Image');
subplot(1,4,2); imshow(Tumor_mask,[]); title('Tumor Mask');
subplot(1,4,3); imshow(brainmask,[]); title('Brain Mask');
subplot(1,4,4); imshow(mask,[]); title('Different Masks');
%%
location_tumor = [98,80];
Input_power=15;
[p] = getParam_2DHeating_Brain_Jacobian_v2(1,1,Nx,Ny,mask,location_tumor,Input_power);

sigma = 1e-7; %scaling for our nonlinear term.
F_jac = @(T,t) evalSysandJac(T,p,sigma);
%forward function which evaluates our system and Jaobian
F_nojac = @(T,t) eval_f_nonlinearSystem(T,p,sigma);
%forward function which evaluates just our system.

%Computing reference steady state solution
tol = 1e-8;
maxiters = 50;
x_ref = NewtonSolve(F_jac,zeros(size(p.A,1),1),p,sigma,tol,maxiters);
%% 1) Determining value of delta t for which FE goes unstable
close all;
x_start = zeros(size(p.A,1),1); 
timestep=0.001; %FE time step
t_start = 0; %FE start time
t_stop= 240; %FE stop time
max_temp=10; %Maximum temperature we will allow our system to reach
visualize=1; %Flag indicating wheter we wish to visualize
step_interval_show=1000;   
%% 2) Implementing TRAP solve and showing when it is faster
dt_trap = 1;
t_stop = 300;
tolerances_newt = [1e-10,1e-10,1e0];
visualize=1;
max_temp=10;
tic
% trap_x = trap(x_start,dt_trap,t_stop,F_jac,tolerances_newt,x_ref);
[trap_x,t_all] = trap_withVisualization(x_start,dt_trap,t_stop,F_jac,tolerances_newt,x_ref,T1_image,brainmask,max_temp,visualize);
trap_time = toc;
[ Period_acc ] = Heating_period_cal( trap_x,t_all,7,mask,180, 120);
 
accuracy_trap = norm(x_ref - trap_x(:,end))/norm(x_ref);

fprintf('Accuracy Trap, %f, achieved in %f seconds\n'...
, accuracy_trap,trap_time)

figure
plot(trap_x(:,end))
hold on
plot(x_ref)

%% 3) Implementing Dynamic TRAP
dt_trap_dyn = 1;
t_stop = 300;
dt_max = 5; %10=> max dt such that I still keep desired accuracy
%at a given time
visualize=1;
max_temp=10;
tolerances_newt = [1e-10,1e-10,1e0];

tic
[dyntrap_x,t_all] = dyn_trap_withVisualization(x_start,dt_trap_dyn,t_stop,F_jac,...
    tolerances_newt,x_ref,dt_max,T1_image,brainmask,max_temp,visualize);
dyn_trap_time = toc;
t_all=t_all(t_all~=0);
%%
% load('all.mat')
[ Period_acc ] = Heating_period_cal( dyntrap_x,t_all,4,mask,180, 120);
accuracy_trap = norm(x_ref - dyntrap_x(:,end))/norm(x_ref);
%%
fprintf('Accuracy Trap, %f, achieved in %f seconds\n'...
, accuracy_trap,dyn_trap_time)
figure
plot(trap_x(:,end))
hold on
plot(x_ref)
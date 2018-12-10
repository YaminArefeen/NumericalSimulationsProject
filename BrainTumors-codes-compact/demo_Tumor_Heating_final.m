%In this script, Yamin explores the difference between Forward
%Euler and Trapezoidal integration for our project. 

%% 0) Initializions and Generate Parameters
clear;
close all;
warning off;
% addpath funcs;
Nx=128;   % larger matrix size shows higher resolution
Ny=128;
load('Image_mask_128_128.mat')
set(0,'DefaultFigureWindowStyle','docked');
figure('units','normalized','outerposition',[0 0 1 1]); 
subplot(1,4,1); imshow(T1_image,[]); title('MR Image','FontSize',16);
subplot(1,4,2); imshow(Tumor_mask,[]); title('Tumor Mask','FontSize',16);
subplot(1,4,3); imshow(brainmask,[]); title('Brain Mask','FontSize',16);
subplot(1,4,4); imshow(mask,[]); title('Different Masks','FontSize',16);

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
x_start = zeros(size(p.A,1),1); 

%% 2)  Dynamic TRAP with Full model
close all;
dt_trap_dyn = 1;
t_stop = 2;
dt_max = 1; %10=> max dt such that I still keep desired accuracy
%at a given time
visualize=1;
max_temp=10;
tolerances_newt = [1e-10,1e-10,1e0];

tic
[dyntrap_x,~] = dyn_trap_withVisualization(x_start,dt_trap_dyn,t_stop,F_jac,...
    tolerances_newt,dt_max,T1_image,brainmask,max_temp,visualize);
time_ori=toc;
disp(['Simulation time of Full-system for 2 iteration: ',num2str(time_ori),' s']);
%% 3) Non-linear Model Reduction Proper Orthogonal Decomposition 
load('Basis_POD_q20.mat','S','U');
q=20;
dt_trap_dyn = 1;
dt_max = 1; %10=> max dt such that I still keep desired accuracy
max_temp=10;

t_stop = 300;
tolerances_newt = [1e-10,1e-10,1e0];
%Calculate V for the reduced system
p.c=ones(size(p.B,1),1);
p.c=diag(p.c);

mask_brain = p.mask(p.mask>0); 
%Portions of the brain over which we simulate, will be 
%N x 1.  Essentially we are picking out portion of mask
%which just corresponds to decision voxels and not air.

vessel_indices = find(mask_brain>2);
%Indices corresponding to vessel voxels in the temperature
%vector and the Jacobian Matrix

vessel_mask = zeros(size(p.A));
vessel_mask(vessel_indices) = 1; 
%mask all places where we have a vessel 

% We're going to linearize around 0, as our heat fluctuation around that
% point will tend to be trivial.
J = p.A;
for ii = 1:length(vessel_indices)
    vidx = vessel_indices(ii);
    J(vidx,vidx) = J(vidx,vidx) + sigma .* 2;
end

init_solve = J\p.B;
V = U(:,1:q);


m = p;
m.A = V'*J*V;
m.B = V'*p.B;
m.c = V'*p.c;

fprintf('Starting TRapezoidal Integration\n')
%We're going to use the linear model everywhere but at the blood vessels
Fr_jac = @(T,t) evalRedSysandJac(T,m,sigma);
%forward function which evaluates our system and Jaobian

xr_start = zeros(size(m.A,1),1); 

close all;
visualize=1;
interval_steps_showResults = 3;
tic;
[trap_xr,t_all] = dyn_trap_withVisualization_MOR(xr_start,dt_trap_dyn,t_stop,Fr_jac,...
    tolerances_newt,dt_max,T1_image,brainmask,max_temp,m.c,interval_steps_showResults,visualize);
trapr_time = toc;
disp(['Simulation time of Reduced-system for 300 iteration: ',num2str(trapr_time),' s']);

trap_xr_y = trap_xr;
trap_xr_y = m.c'*trap_xr_y;

%% 4) Compare reduced model vs. full model
close all;
load('Im_compare_ModelReduction.mat')
figure('units','normalized','outerposition',[0 0 1 1]); 
subplot(1,3,1);
imshow(Temperature_ref,[0,15]); colormap('hot');
title('Original model','FontSize',10);
subplot(1,3,2);
imshow(Temperature_im,[0,15]); colormap('hot');
title('Non-linear MOR with POD (q=20)','FontSize',10);
subplot(1,3,3);
imshow(1000*abs(Temperature_ref-Temperature_im),[0,15]); colormap('hot');
title('1000X difference','FontSize',10);
%% 5) Evaluation map
close all;
figure('units','normalized','outerposition',[0 0 1 1]); 
imshow(T1_image,[]); colormap('gray');
title('MR T1-weighted image','FontSize',14);

pause(1);

figure('units','normalized','outerposition',[0 0 1 1]); 
subplot(1,2,1);
imshow(Temperature_im,[0,10]); colormap('hot'); colorbar;
title('Temperature at 300s','FontSize',14);

[ Period_acc ] = Heating_period_cal( trap_xr_y,t_all,4,mask,180, 120);

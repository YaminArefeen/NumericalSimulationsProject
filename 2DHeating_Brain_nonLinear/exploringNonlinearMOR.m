%This script demonstrates the ability of a reduced system to speed up the solver for
%acceptable accuracy

%% 0) Initializions
% clear;
clear all;
close all;
addpath funcs;
load('MRI_tumor_image.mat');
Nx=64;   % larger matrix size shows higher resolution
Ny=64;

T1_image = real(ifft2c(crop(fft2c(double(T1_image)),[Nx,Ny])));
Tumor_mask = real(ifft2c(crop(fft2c(double(Tumor_mask)),[Nx,Ny])));
brainmask = real(ifft2c(crop(fft2c(double(brainmask)),[Nx,Ny])));

Tumor_mask=Tumor_mask>5*mean(Tumor_mask(:));
brainmask=brainmask>1*mean(brainmask(:));
figure; 
subplot(1,3,1); imshow(T1_image,[]); title('MR Image');
subplot(1,3,2); imshow(Tumor_mask,[]); title('Tumor Mask');
subplot(1,3,3); imshow(brainmask,[]); title('Brain Mask');

mask=zeros(size(T1_image));
mask(Tumor_mask)=1;              % 1 indicates tumor, 0 is air
mask(logical(brainmask-Tumor_mask))=2;    % 2 indicates health tissues
%select input evaluation functions
eval_u = 'eval_u_step';

figure; imshow(mask,[]);
load('mask.mat')

location_tumor = [47,41];
[p] = getParam_2DHeating_Brain_Jacobian(1,1,Nx,Ny,mask,location_tumor);

sigma = 1e-2; %scaling for our nonlinear term.
F_jac = @(T,t) evalSysandJac(T,p,sigma);
%forward function which evaluates our system and Jaobian
F_nojac = @(T,t) eval_f_nonlinearSystem(T,p,sigma);
%forward function which evaluates just our system.

%Computing reference steady state solution
tol = 1e-8;
maxiters = 50;
x_ref = NewtonSolve(F_jac,zeros(size(p.A,1),1),p,sigma,tol,maxiters);

%The initial condition
x_start = zeros(size(p.A,1),1); 

%% 2) Demonstrate the TRAP solver for the entire system
dt_trap = 1;
t_stop = 300;
tolerances_newt = [1e-10,1e-10,1e0];

tic
trap_x = trap(x_start,dt_trap,t_stop,F_jac,tolerances_newt,x_ref);
trap_time = toc;

accuracy_trap = norm(x_ref - trap_x(:,end))/norm(x_ref);

fprintf('Accuracy Trap, %f, achieved in %f seconds\n'...
, accuracy_trap,trap_time)

figure
plot(trap_x(:,end))
hold on
plot(x_ref)

%% 3) Perform Model Linearization

q=10;
F_redJac = @(T,t) evalRedSysandJacNon(T, p, sigma, q);


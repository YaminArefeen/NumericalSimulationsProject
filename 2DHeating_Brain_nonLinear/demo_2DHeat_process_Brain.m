% copyright Luca Daniel, MIT 2018
% modified by Zijing 2018/10/21
clear;
close all;
addpath funcs;
load('MRI_tumor_image.mat');
Nx=64;   % larger matrix size shows higher resolution
Ny=64;
%% down-sampling the image to lower resolution
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
% 2D Heat Process
%% Generate vessel mask
% figure; imshow(T1_image,[]); title('Creat the vessel mask');
% IFH = imfreehand();
% vessel_mask = IFH.createMask();
% figure; imshow(vessel_mask,[]);
% mask(vessel_mask)=3;   % 3 indicate vessel in brain
% 
% figure; imshow(mask==3,[]);
% location_vessel_edge=ginput(2);
% location_vessel_edge=round(location_vessel_edge);
% mask(location_vessel_edge(1,2),location_vessel_edge(1,1))=4;   % 4 indicate vessel edge in brain
% mask(location_vessel_edge(2,2),location_vessel_edge(2,1))=4;   % 4 indicate vessel edge in brain
% mask(29,50)=2;
figure; imshow(mask,[]);
load('mask.mat')
%%
% location_tumor = [42,42];
location_tumor = [47,41];
% location_tumor = [13,11]; %LOWER RES LOCATION (16,16) for testing
[p] = getParam_2DHeating_Brain_Jacobian(1,1,Nx,Ny,mask,location_tumor);
%%
sigma = 1e-2; %scaling for our nonlinear term.
F = @(T) eval_f_nonlinearSystem(T,p,sigma);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Initializing Forward Euler Parameters
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
x_start = zeros(size(p.A,1),1); 
%Initial tempeatures, assuming we are starting all temperatures
%at 0, whatever our convention for 0 is.
timestep=0.001; %FE time step
t_start = 0; %FE start time
t_stop=1; %FE stop time
max_temp=10; %Maximum temperature we will allow our system to reach
visualize=20; %Flag indicating wheter we wish to visualize
step_interval_show=200;   
% the time interval of image update ( 1 updata/ N steps )
[X] = ForwardEuler2D_Brain(F,x_start,t_start,t_stop,timestep,T1_image,brainmask,max_temp,step_interval_show,visualize);
%% Checking Analytic and Finite differences Jacobian
eps = 1; %perturbation for FID jacobian.
% x_test = rand(size(x_start)); %try random inputs for testing
% x_test = zeros(size(x_start));
x_test = X;

Jf_analytic = eval_Jf_nonlinearSystem(x_test,p,sigma);
Jf_Dif = eval_J_fid(F,x_test,eps);

%calculating error based on a voxel by voxel basis
error=Jf_analytic-Jf_Dif;
disp(mean(error(:)))

%calculating error by vectorizing our Jacobians
error2 = norm(Jf_analytic(:) - Jf_Dif(:))/norm(Jf_analytic(:));
disp(error2)



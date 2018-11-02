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
T=zeros(Nx,Ny);
[F,J,p,x_start,t_start,t_stop,max_dt_FE] = getParam_2DHeating_Brain_Jacobian(T+273+37,1,1,Nx,Ny,mask,location_tumor);

Func = @(T)  getParam_2DHeating_Brain_Jacobian(T,1,1,Nx,Ny,mask,location_tumor);
%%
% eval_f = 'eval_Jf_nonLinearSystem';
eval_f = 'eval_f_nonLinearSystem';

% eval_f = 'eval_f_LinearSystem';

% figure; imagesc(p.A);
p.Nx=Nx;
p.Ny=Ny;
% test FE function
timestep=0.001;
t_stop=100000;
max_temp=10;
visualize=20;
step_interval_show=1;   % the time interval of image update ( 1 updata/ N steps )
[X] = ForwardEuler2D_Brain(eval_f,x_start,eval_u,Func,t_start,t_stop,timestep,T1_image,brainmask,max_temp,step_interval_show,visualize);
%% Jacobian Calculation
Jf_analytic = eval_Jf_LinearSystem(x_start,1,p);
dx=0.001;
N=size(p.A,1);
Jf_Dif=zeros(N);
for i=1:N
    x_new=x_start;
    x_new(i)=x_new(i)+dx;
    temp=(eval_f_LinearSystem(x_new,1,p)-eval_f_LinearSystem(x_start,1,p))./dx;
    Jf_Dif(:,i)=temp;
end
error=Jf_analytic-Jf_Dif;
disp(mean(error(:)));



function X = ForwardEuler2D_Brain(eval_f,x_start,eval_u,Func,t_start,t_stop,timestep,T1_image,brainmask,max_temp,step_interval_show,visualize)
% uses Forward Euler to simulate states model dx/dt=f(x,u,p)
% from state x_start at time t_start
% until time t_stop, with time intervals timestep
% eval_f is a string including the name of the function that evaluates f(x,u,p)
% eval_u os a string including the name of the funciton that evaluates u(t)
% 
% X = ForwardEuler(eval_f,x_start,eval_u,p,t_start,t_stop,timestep)

% copyright Luca Daniel, MIT 2018
X(:,1) = x_start;
t(1) = t_start;
if visualize
    T1_image=T1_image-min(T1_image(:));
    T1_image=T1_image./max(T1_image(:));
    figure(2)
    imshow(T1_image,[0 1],'InitialMag', 'fit'); colormap('gray');
    hold on
    Temperature_im=zeros([size(T1_image),3]);
end

num=1;
for n=1:ceil((t_stop-t_start)/timestep),
   dt = min(timestep, (t_stop-t(n)));
   t(n+1)= t(n) + dt;
%    u = feval(eval_u, t(n));
%    f = feval(eval_f, X(:,n), u, p);
    tmp=zeros(size(T1_image));
    tmp(brainmask)=X(:,n);
    f = Func(tmp);
    X(:,n+1)= X(:,n) +  dt * f;
   num=num+1;
   if mod(num,step_interval_show)==0
   if visualize
       Temperature_im(brainmask)=X(:,n+1);
       imshow(T1_image,[0 1],'InitialMag', 'fit'); colormap('gray');
       hold on
       h=imshow(Temperature_im);
       set(h, 'AlphaData', Temperature_im(:,:,1)/max_temp); title(['Time at ',num2str(t(n+1))]);
       drawnow;       
       hold on
   end
   end
end

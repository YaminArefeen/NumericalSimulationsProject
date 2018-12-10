function X = ForwardEuler2D(eval_f,x_start,eval_u,p,t_start,t_stop,timestep,visualize)
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
   visualizeResults2D(t,X,1,p.Nx,p.Ny,'.b');
end
for n=1:ceil((t_stop-t_start)/timestep),
   dt = min(timestep, (t_stop-t(n)));
   t(n+1)= t(n) + dt;
   u = feval(eval_u, t(n));
   f = feval(eval_f, X(:,n), u, p);
   X(:,n+1)= X(:,n) +  dt * f;
   if visualize
      visualizeResults2D(t,X,n+1,p.Nx,p.Ny,'.b');
   end
end
end

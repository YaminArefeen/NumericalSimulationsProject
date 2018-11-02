function u = eval_u_step(t);
% generates the value of the input at time t
% corresponding to a step of magnitude a 
%
% u = eval_u_step(t);

% copyright Luca Daniel, MIT 2018

if t <0
   u = 0;
else 
   u = 1;
end


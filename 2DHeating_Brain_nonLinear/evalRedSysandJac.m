function [f,J] ...
= evalRedSysandJac(T,p,sigma)
%Evaluates our system and its Jacobian, for the sake of how I 
%structured our trapezoidal code.

mask_brain = p.mask(p.mask>0); 
%Portions of the brain over which we simulate, will be 
%N x 1.  Essentially we are picking out portion of mask
%which just corresponds to decision voxels and not air.

vessel_indices = find(mask_brain>2);
%Indices corresponding to vessel voxels in the temperature
%vector and the Jacobian Matrix

vessel_mask = zeros(size(T(:)));
vessel_mask(vessel_indices) = 1; 
%mask all places where we have a vesseel 

J = p.A;
% for ii = 1:length(vessel_indices)
%     vidx = vessel_indices(ii);
%     J(vidx,vidx) = ...
%         J(vidx,vidx) + p.dx2.*...
%             (sigma*exp(sigma.*T(vidx)) ... 
%                 + sigma*exp(-sigma.*T(vidx)));
% end
%Our Jacobian will simply be the original nodal matrix + 
%the nonlinear term at voxels corresponding to vessels.
%Since the nonlinearity only applies on a voxel by voxel
%basis and has no influence on connections between voxels,
%we only need to add the nonlinearity to the diagonal elements
%of the Jacobian which correspond to vessels.

% f = p.A*T(:)+ vessel_mask.*...
%     (p.dx2*(exp(sigma.*T(:))-exp(-sigma.*T(:))))+p.B;
f = p.A*T(:) + p.B;

%Evaluating the linear and nonlienar portions of our system
%seperately. Only evaluate nonlinearity at vessel locations,
%based on the convention of nonlineary we've defined.

end
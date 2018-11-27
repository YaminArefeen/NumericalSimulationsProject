function [f,J] = evalRedSysandJacNon(T,p,sigma)
%Evaluates our system and its Jacobian, for the sake of how I 
%structured our trapezoidal code.

A = p.A;
B = p.B;

mask_brain = p.mask(p.mask>0); 
%Portions of the brain over which we simulate.

vessel_indices = find(mask_brain>2);
%Indices corresponding to vessel voxels in the temperature
%vector and the Jacobian Matrix

vessel_mask = zeros(size(p.A));
vessel_mask(vessel_indices) = 1;

Amat = {A; A+vessel_mask.*p.dx2.*(2*sigma)};

J = Amat{1};
for ii = 1:length(vessel_indices)
    vidx = vessel_indices(ii);
    J(vidx,vidx) = ...
        J(vidx,vidx) + p.dx2.*2*sigma*T(vidx);
end

Kmat = {0; -vessel_mask.*J*T(:)};

f = Amat{1}*T(:) + Kmat{1} + sigma*(Amat{2}*T(:) + Kmat{2}) + B;
J = Amat{1};%(f-(Amat{1}*(T(:)+ones(size(T)).*eps) + Kmat{1} + Amat{2}*(T(:)+ones(size(T)).*eps) + Kmat{2} + B))/eps;



%Our Jacobian will simply be the original nodal matrix + 
%the nonlinear term at voxels corresponding to vessels.
%Since the nonlinearity only applies on a voxel by voxel
%basis and has no influence on connections between voxels,
%we only need to add the nonlinearity to the diagonal elements
%of the Jacobian which correspond to vessels.

% f = p.A*T(:)+ vessel_mask.*...
%     (p.dx2*(exp(sigma.*T(:))-exp(-sigma.*T(:))))+p.B;

%Evaluating the linear and nonlienar portions of our system
%seperately. Only evaluate nonlinearity at vessel locations,
%based on the convention of nonlineary we've defined.



end
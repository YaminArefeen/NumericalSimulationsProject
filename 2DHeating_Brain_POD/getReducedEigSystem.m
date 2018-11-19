function [Ahat, bhat, chat] = getReducedEigSystem(A,b,c,q)
    % Get the eigen properties of the system
    [V,D] = eig(full(A));
%     diag(D)
    
    % Select the q slowest nodes
    Vq = V(:,end-q+1:end);
    size(Vq)
    
    % Adjust the system
    chat = Vq'*c;
    bhat = Vq'*b;
    Ahat = Vq'*A*Vq;

end
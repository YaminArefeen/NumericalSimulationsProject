function [Ahat, bhat, chat] = getReducedEigSystem(A,b,c,q,type_trunc)
    % Get the eigen properties of the system
    [V,D] = eig(full(A));
    
    if type_trunc == 0
        % Select smallest residues
        btild = V'*b;
        ctild = V'*c;
        checker = btild.*ctild;
        [out,idx] = sort((checker),'ascend');
        Vq = V(:,idx(1:q));
    end
    
    if type_trunc == 1
        % Select smallest quotient
        btild = V'*b;
        ctild = V'*c;
        checker = btild.*ctild./diag(D);
        [out,idx] = sort(abs(checker),'descend');
        Vq = V(:,idx(1:q));
    end
    
    if type_trunc == 2
        %Select the q slowest nodes
        Vq = V(:,end-q+1:end);
    end
    
    if type_trunc == 3
        %Select the q fastest nodes
        Vq = V(:,1:q);
    end
    
    if type_trunc == 4
        % Select half slow and half fast
        Vq(:,1:q/3) = V(:,1:q/3);
        Vq(:,q/3+1:q)= V(:,end-2*q/3+1:end);
    end
    
    % Adjust the system
    chat = Vq'*c;
    bhat = Vq'*b;
    Ahat = Vq'*A*Vq;

end
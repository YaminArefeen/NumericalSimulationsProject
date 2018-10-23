function [E,Is,alpha] = constructEIs(R,Isource,N)
    %Constructing E

    C = size(R,2); 
    S = size(Isource,2); %number of nodes, components, and sources respectively

    E = zeros(N,C); %incidence matrix
    alpha = eye(C);
    Is = zeros(N,1); %rhs current matrix

    %loop through each component of R and fill in E
    for ii = 1:C
        alpha(ii,ii) = 1/R(1,ii);   
        n1 = R(2,ii);
        n2 = R(3,ii);

        if(n1 ~= 0)
            E(n1,ii) = 1;
        end

        if(n2 ~= 0)
            E(n2,ii) = -1;
        end
    end

    %loop through each component of Isource and fill in Is
    for ii = 1:S
        current = Isource(1,ii);
        n1 = Isource(2,ii);
        n2 = Isource(3,ii);

        if(n1 ~= 0)
            Is(n1) = Is(n1) - current;       
        end

        if(n2 ~= 0)
            Is(n2) = Is(n2) + current;       
        end
    end
end
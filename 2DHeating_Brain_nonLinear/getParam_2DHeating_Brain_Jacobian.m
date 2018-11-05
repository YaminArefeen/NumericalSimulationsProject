function [p] = getParam_2DHeating_Brain_Jacobian(dx,dy,Nx,Ny,mask,location_tumor)

% Luca Daniel, MIT 2018 heatbar
% modified by Zijing 2018/10/21
% mask for brain tissues [nx ny]:  0 air, 1 tumor, 2 health tissues
% location_tumor, (x,y) locate the heat point for tumor

gamma =0.1; %thermal capacitance per unit length of the bar
km    =0.1; %thermal conductance through metal per unit length of the bar
Rair  =100; %thermal resistance to air 
RcTissue     = (1/km)*dx*dy*1; %resistance of tissue
RcTumor       = (1/km)*dx*dy*0.1;
RcVessel_air       = (1/km)*dx*dy*0.01;
RcVessel_V       = (1/km)*dx*dy*0.01;

Cstore = gamma*dx*dy; %the longer the section the larger the storage
     
%% Prepare a list that only contains real tissues (make sure A matrix does not have voltage source/air)
tissue_list = mask>0;     
N = sum(tissue_list(:));  % total number of nodes, air pixels are just ground
index_image = zeros(size(mask));
index_image(tissue_list)=1:N;    % the index of nodes
p.A = zeros(N,N);
% coupling resistors Rc between (x,y) and (x+1,y), (x,y+1), also dectect
% air pixels around (x,y)
for x = 2:Nx-1
    for y=2:Ny-1
        if (mask(x,y)>0)
            index=index_image(x,y);
            
            index_1 = index_image(x+1,y);  % 1 and 2 are using for building resistors and air detection, 
            index_2 = index_image(x,y+1);  
            index_3 = index_image(x,y-1);  % 3 and 4 are just for air detection
            index_4 = index_image(x-1,y);
            
            if(mask(x,y)==1)
                Rc = RcTumor;
            elseif (mask(x,y)==2)
                Rc = RcTissue;
            elseif (mask(x,y)==3 || mask(x,y)==4)
                Rc = RcVessel_V;
            end
            
            % build resistor connection
            if (mask(x+1,y)>0)
%             if (mask(x+1,y)>0)&&( (mask(x+1,y)<3) || ((mask(x,y)>0)&&(mask(x,y)<3)))
                Rc = resistor_lookup(mask(x,y), mask(x+1,y),RcTumor, RcTissue,RcVessel_V);
                p.A(index,index) = p.A(index,index)+(+1/Rc);
                p.A(index,index_1) = p.A(index,index_1)+(-1/Rc);
                p.A(index_1,index) = p.A(index_1,index)+(-1/Rc);
                p.A(index_1,index_1) = p.A(index_1,index_1)+(+1/Rc);
            end
            if (mask(x,y+1)>0)
%             if (mask(x,y+1)>0)&&( (mask(x,y+1)<3) || ((mask(x,y)>0)&&(mask(x,y)<3)))
                Rc = resistor_lookup(mask(x,y), mask(x,y+1),RcTumor,RcTissue,RcVessel_V);
                p.A(index,index) = p.A(index,index)+(+1/Rc);
                p.A(index,index_2) = p.A(index,index_2)+(-1/Rc);
                p.A(index_2,index) = p.A(index_2,index)+(-1/Rc);
                p.A(index_2,index_2) = p.A(index_2,index_2)+(+1/Rc);
            end
            
%             if (mask(x+1,y)>2)&&(mask(x,y)>2)  
%                 p.A(index,index) = p.A(index,index)+(+1/RcVessel_V);
%                 p.A(index,index_1) = p.A(index,index_1)+(-1/RcVessel_V);
%                 p.A(index_1,index) = p.A(index_1,index)+(-1/RcVessel_V);
%                 p.A(index_1,index_1) = p.A(index_1,index_1)+(+1/RcVessel_V);
%             end
%             if (mask(x,y+1)>2)&&(mask(x,y)>2)  
%                 p.A(index,index) = p.A(index,index)+(+1/RcVessel_V);
%                 p.A(index,index_2) = p.A(index,index_2)+(-1/RcVessel_V);
%                 p.A(index_2,index) = p.A(index_2,index)+(-1/RcVessel_V);
%                 p.A(index_2,index_2) = p.A(index_2,index_2)+(+1/RcVessel_V);
%             end
            
            if index_1==0  
                p.A(index,index) = p.A(index,index) + 1/Rair;
            end
            if index_2==0  
                p.A(index,index) = p.A(index,index) + 1/Rair;
            end
            if index_3==0  
                p.A(index,index) = p.A(index,index) + 1/Rair;
            end
            if index_4==0  
                p.A(index,index) = p.A(index,index) + 1/Rair;
            end
            
            if mask(x,y)==4  
                p.A(index,index) = p.A(index,index) + 1/RcVessel_air;
            end            
        end
    end
end
% leakage resistor Rair between edge points and air

p.B     = zeros(N,1);
idx=index_image(location_tumor(1),location_tumor(2));
p.B(idx,1)= 100;

p.A     = -p.A/Cstore; % note this will give a 1/dz^2 in A
							  % also pay attention to the negative sign
p.B     = p.B/Cstore;  % note this is important to make sure results
                       % will not depend on the number of sections N
p.mask  = mask;     
p.dx2 = dx*dy;

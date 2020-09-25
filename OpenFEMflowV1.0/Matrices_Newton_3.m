function [R, dRdU, dRdX, dRdY, dRdt] = Matrices_Newton_3(U,nE,nodes,element,nDOF,dt,CFL,BC,theta,Adjoint)

R = zeros(nDOF,1);
dRdU = sparse(nDOF,nDOF);

if Adjoint ==1
    dRdX = []; %sparse(nDOF,size(nodes,1));
    dRdY = sparse(nDOF,size(nodes,1));
    dRdt = sparse(nDOF,size(nodes,1));
else
    dRdX = [];
    dRdY = [];
    dRdt = [];
end
   

for i = nE(1):nE(2) 
    xe = [nodes(element(i,1),1) nodes(element(i,2),1)  nodes(element(i,3),1)];
    ye = [nodes(element(i,1),2) nodes(element(i,2),2)  nodes(element(i,3),2)];

    Ue = [U(4*(element(i,1)-1)+1:4*(element(i,1)-1)+4); U(4*(element(i,2)-1)+1:4*(element(i,2)-1)+4); ...
          U(4*(element(i,3)-1)+1:4*(element(i,3)-1)+4)];


    T = eye(12);    
    dTdt{1} = zeros(12);
    dTdt{2} = zeros(12);
    dTdt{3} = zeros(12);

    for n = 1:3
        c = (element(i,n)==BC(BC(:,2)==1));
        if sum(c)>=1

            t = theta(theta(:,1)==element(i,n),2);
            T(4*(n-1)+2:4*(n-1)+3,4*(n-1)+2:4*(n-1)+3) = [cos(t) sin(t); -sin(t) cos(t)];

            if Adjoint==1
               dTdt{n}(4*(n-1)+2:4*(n-1)+3,4*(n-1)+2:4*(n-1)+3) = [-sin(t) cos(t); -cos(t) -sin(t)];
            end

        end
    end

    Uexy = T'*Ue;

    if Adjoint ==0
        [Re1, dRedU1] = ElementMatrics_Newton_3(xe,ye,Uexy,dt,CFL,Adjoint);

    elseif Adjoint ==1
        [Re1, dRedU1, ~, dRedY] = ElementMatrics_Newton_3(xe,ye,Uexy,dt,CFL,Adjoint);

    end


    Re = T*Re1;
    dRedU = T*dRedU1*T';


    I = [4*(element(i,1)-1)+1 4*(element(i,1)-1)+2 4*(element(i,1)-1)+3 4*(element(i,1)-1)+4 ...
         4*(element(i,2)-1)+1 4*(element(i,2)-1)+2 4*(element(i,2)-1)+3 4*(element(i,2)-1)+4 ...
         4*(element(i,3)-1)+1 4*(element(i,3)-1)+2 4*(element(i,3)-1)+3 4*(element(i,3)-1)+4];

    R(I) = R(I) + Re;
    dRdU(I,I) = dRdU(I,I) + dRedU;

    if Adjoint ==1
        Ie = [element(i,1) element(i,2) element(i,3)];  

        dRedY = T*dRedY;      
        dRedt = [dTdt{1}*Re1 dTdt{2}*Re1 dTdt{3}*Re1] + [T*dRedU1*dTdt{1}'*Ue   T*dRedU1*dTdt{2}'*Ue   T*dRedU1*dTdt{3}'*Ue];         

        dRdY(I,Ie) = dRdY(I,Ie) + dRedY;
        dRdt(I,Ie) = dRdt(I,Ie) + dRedt;
    end
end






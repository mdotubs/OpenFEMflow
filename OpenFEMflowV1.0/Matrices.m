function [R, dRdU] = Matrices(U,nE,nodes,element,nDOF,dt,CFL,BC,theta)

IR = zeros(12*(nE(2)-nE(1)+1),1);
KR = IR;
ndoublet = 0;

IdR = zeros(12*12*(nE(2)-nE(1)+1),1);
JdR = IdR;
KdR = IdR;
ntriplets = 0 ;

for i = nE(1):nE(2) 
    xe = [nodes(element(i,1),1) nodes(element(i,2),1)  nodes(element(i,3),1)];
    ye = [nodes(element(i,1),2) nodes(element(i,2),2)  nodes(element(i,3),2)];

    Ue = [U(4*(element(i,1)-1)+1:4*(element(i,1)-1)+4); U(4*(element(i,2)-1)+1:4*(element(i,2)-1)+4); ...
          U(4*(element(i,3)-1)+1:4*(element(i,3)-1)+4)];


    T = eye(12);
    for n = 1:3
        c = (element(i,n)==BC(BC(:,2)==1));
        if sum(c)>=1

            t = theta(theta(:,1)==element(i,n),2);
            T(4*(n-1)+2:4*(n-1)+3,4*(n-1)+2:4*(n-1)+3) = [cos(t) sin(t); -sin(t) cos(t)];

        end
    end

    Uexy = T'*Ue;

    [Re, dRedU] = ElementMatrics_3(xe,ye,Uexy,dt,CFL);

    dRedU = T*dRedU*T';
    Re = T*Re;

    Indx = [4*(element(i,1)-1)+1 4*(element(i,1)-1)+2 4*(element(i,1)-1)+3 4*(element(i,1)-1)+4 ...
            4*(element(i,2)-1)+1 4*(element(i,2)-1)+2 4*(element(i,2)-1)+3 4*(element(i,2)-1)+4 ...
            4*(element(i,3)-1)+1 4*(element(i,3)-1)+2 4*(element(i,3)-1)+3 4*(element(i,3)-1)+4];
  
                      
    for krow=1:12
        ndoublet = ndoublet + 1;
        IR(ndoublet) = Indx(krow);
        KR(ndoublet) = Re(krow);
            for kcol=1:12
                ntriplets = ntriplets + 1 ;
                IdR (ntriplets) = Indx (krow) ;
                JdR (ntriplets) = Indx (kcol) ;
                KdR (ntriplets) = dRedU (krow,kcol) ;
            end
    end

end

R = sparse(IR,ones(12*(nE(2)-nE(1)+1),1),KR,nDOF,1);
dRdU = sparse(IdR,JdR,KdR,nDOF,nDOF);

    

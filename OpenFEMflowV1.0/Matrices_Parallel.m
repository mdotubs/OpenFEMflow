function [R, dRdU, dRdX, dRdY, dRdtheta] = Matrices_Parallel(U,Method)

global dt CFL Parallel Adjoint
global element nodes nDOF BC theta


nElm = size(element,1);

nCore = feature('numcores');
nP = round(nElm/nCore);

dRdU = sparse(nDOF,nDOF);
R = sparse(nDOF,1);

if Adjoint ==1
  dRdX = sparse(nDOF,size(nodes,1));
  dRdY = sparse(nDOF,size(nodes,1));
  dRdtheta = sparse(nDOF,size(nodes,1));
end

if Parallel ==1

    parfor i=1:nCore 
        if strcmpi(Method,'Picard')
            [Rp{i}, dRdUp{i}] = Matrices(U,[(i-1)*nP+1,min((i-1)*nP+nP,nElm)],nodes,element,nDOF,dt,CFL,BC,theta);

        elseif strcmpi(Method,'Newton')
            [Rp{i}, dRdUp{i}, ~, dRdYp{i}, dRdthetap{i}] = Matrices_Newton(U,[(i-1)*nP+1,min((i-1)*nP+nP,nElm)],nodes,element,nDOF,dt,CFL,BC,theta,Adjoint);
        end
    end      
        

    for i=1:nCore
        R = R +Rp{i};
        dRdU = dRdU + dRdUp{i};
        if Adjoint ==1
           dRdY = dRdY + dRdYp{i};
           dRdtheta = dRdtheta + dRdthetap{i};
        end    
    end
    clear Kp Fp Rp dRdUp
    
elseif Parallel ==0

    if strcmpi(Method,'Picard')
        [R, dRdU] = Matrices(U,[1 nElm],nodes,element,nDOF,dt,CFL,BC,theta);

    elseif strcmpi(Method,'Newton')
        [R, dRdU, dRdX, dRdY, dRdtheta] = Matrices_Newton(U,[1 nElm],nodes,element,nDOF,dt,CFL,BC,theta,Adjoint);

    end            

end



    

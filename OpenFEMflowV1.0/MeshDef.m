function [nodes2, theta, dDeltaX_dA, dDeltaY_dA, dtheta_dA] = MeshDef(nodes,element,A)

global BC Surf

nodes2 = nodes;
nNodes = size(nodes,1);
nElm = size(element,1);
DeltaX = zeros(nNodes,1);
DeltaY = zeros(nNodes,1);

dDeltaX_dA = sparse(nNodes,length(A));
dDeltaY_dA = sparse(nNodes,length(A));

%% nodes on the noundaries

B1u = Surf(Surf(:,2)==1);
B1l = Surf(Surf(:,2)==2);
B2 = BC(BC(:,2)==2);

B = [B1u; B1l; B2];       % boundary nodes

I = setdiff(1:nNodes,B);

%% wall nodes movement

Au = A(1:length(A)/2);
Al = A(length(A)/2+1:end);

[yu, dyudx, dyudA, d_dydx_dAu] = CSTairfoil(Au,nodes(B1u,1)'-nodes(B1u(1),1)); % -nodes(B1u(1),1) =+0.5
[yl, dyldx, dyldA, d_dydx_dAl] = CSTairfoil(Al,nodes(B1l,1)'-nodes(B1u(1),1)); % +0.5

DeltaY(B1u) = yu' - nodes(B1u,2);
DeltaY(B1l) = yl' - nodes(B1l,2);

theta = [B1u atan(dyudx'); B1l atan(dyldx')];

dDeltaY_dA(B1u,1:length(A)/2) = dyudA';
dDeltaY_dA(B1l,length(A)/2+1:end) = dyldA';

dtheta_dA= sparse(nNodes,length(A));

for i=1:length(A)/2
    dtheta_dA(B1u,i) = d_dydx_dAu(i,:)./(1+dyudx.^2);
    dtheta_dA(B1l,length(A)/2+i) = d_dydx_dAl(i,:)./(1+dyldx.^2);
end
dtheta_dA(B1u(1),:) = 0;    % LE

%% finding connecting nodes

for i=1:nNodes

    e = find(element==i);   
 
    clear E     
    for ne=1:length(e)
        if e(ne) <= nElm
            E(ne) = e(ne);
        elseif e(ne) > nElm && e(ne) <=2*nElm
            E(ne) = e(ne)-nElm;
        elseif e(ne) > 2*nElm
            E(ne) = e(ne)-2*nElm;
        end
    end
    
    N = [];
    for j=1:length(E)
        for n=1:3
            if isempty(find(element(E(j),n)==N,1))
                N(end+1) = element(E(j),n);
            end
        end
    end
    
    Nc{i} = setdiff(N,i);
    
end

%% Spring anaogy 

tol = 1e-3;

deltaX = DeltaX;
deltaY = DeltaY;
deltaXn = deltaX;
deltaYn = deltaY;

ddeltaX_dDeltaX = sparse(nNodes,nNodes);
ddeltaY_dDeltaY = sparse(nNodes,nNodes);

for iter=1:200
    for i=1:nNodes
        if isempty(find(i==B,1))

            k = 1./sqrt((nodes(i,1)-nodes(Nc{i},1)).^2+(nodes(i,2)-nodes(Nc{i},2)).^2);

            deltaXn(i) = sum(k.*deltaX(Nc{i}))/sum(k);
            deltaYn(i) = sum(k.*deltaY(Nc{i}))/sum(k);  
        end
    end
    
    if norm(deltaXn-deltaX)<=tol && norm(deltaYn-deltaY)<=tol
        ddeltaX_dDeltaX(i,Nc{i}) = k/sum(k);
        ddeltaY_dDeltaY(i,Nc{i}) = k/sum(k);
        break;
    else
        deltaX = deltaXn;
        deltaY = deltaYn;
    end
end

DeltaX = deltaXn;
DeltaY = deltaYn;

nodes2(:,1) = nodes(:,1) + DeltaX;
nodes2(:,2) = nodes(:,2) + DeltaY;
 
dDeltaX_dA(I,:) =   ddeltaX_dDeltaX(I,:)*dDeltaX_dA;       
dDeltaY_dA(I,:) =   ddeltaY_dDeltaY(I,:)*dDeltaY_dA;       
    
    
    
    
    
    

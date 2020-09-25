function [Xopt, Fopt, X_iter] = SQP(X0,Alpha,iter_max_sqp) 

global Sol0 t0 Options Mach Sol

%%% Optimization
%     min    Cd/Cd0
%   x,alpha
%     s.t.   Cl/Cl0 - 1 = 0
%            Cm/Cm0 - 1 <= 0   -->  Cm0 - Cm + s(1)^2 = 0
%            1 - t/t0  <= 0   -->   t0 - t  + s(2:3).^2 = 0                      

Cl0 = Sol0.Cl;
Cd0 = Sol0.Cd;
Cm0 = Sol0.Cm;

X = [X0'; Alpha; 0 ; 0 ; 0];  % three slack variables

%% SQP settings

mu = 1e-4;
tol_kkt = 1e-5;
tol_f = 1e-7;
tol_x = 1e-5;
s_max = round(log(1/tol_x)/log(2));

%%

Sol = Sol0;
H = speye(length(X));
X_iter = X';
nf = 0;

for iter = 1: iter_max_sqp
      
    if iter ==1
        [Cl, Cd, Cm, CP, Sol, dCl_dAlpha, dCd_dAlpha, dCm_dAlpha, dCl_dMach, dCd_dMach, dCm_dMach, dCl_dCST, dCd_dCST, dCm_dCST] ...
            = FEMflow(Mach, X(11), 1, Options, X(1:10), Sol);
        [t, dt_dCST] = AirfoilThickness(X(1:10));
        nf = nf+1;
    end
    % ************** KKT matrix *****************
    g = [dCd_dCST';  dCd_dAlpha; 0; 0; 0];
    
    A = [dCl_dCST/Cl0    dCl_dAlpha/Cl0  0 0 0;
         dCm_dCST/Cm0    dCm_dAlpha/Cm0 2*X(12) 0 0;
        -dt_dCST'./t0'    [0; 0]      [0 2*X(13) 0; 0 0 2*X(14)]];
    
    c = [Cl/Cl0-1; Cm/Cm0-1+X(12)^2; 1-t'./t0'+X(13:14).^2];
    
    KKT = [H   -A';
           A   sparse(length(c),length(c))];
    RHS = [-g; -c];
    
    % ************* Newton step ****************
    delta  =  KKT\RHS;
    S      =  delta(1:length(X));
    lambda =  delta(length(X)+1:end);
    
    % ************* step length ****************
    if iter==1
       P = Cd - sum(lambda.*c);

       Q  = g  - sum(meshgrid(lambda,1:length(X))'.*A)';   
       display('Iter      nF      Objfun       norm of Const.        step length         First order optimality   ' );
       display([num2str(0)  '         '  num2str(1) '           '   num2str(Cd/Cd0)  '            '  num2str(norm(c)) '               '  num2str(0) '               '  num2str(norm(Q))]);

    end
    
    a = 1;
    for s = 1:s_max
        
        dn = a*S;
        Xn = X + dn;
               
        [Cln, Cdn, Cmn, ~, Sol, dCln_dAlpha, dCdn_dAlpha, dCmn_dAlpha, ~,~,~, dCln_dCST, dCdn_dCST, dCmn_dCST] ...
            = FEMflow(Mach, Xn(11), 1, Options, Xn(1:10), Sol);
        [tn, dtn_dCST] = AirfoilThickness(Xn(1:10));
        nf = nf+1;
             
        gn = [dCdn_dCST';  dCdn_dAlpha; 0; 0; 0];
        
        An = [dCln_dCST/Cl0    dCln_dAlpha/Cl0 0 0 0 
              dCmn_dCST/Cm0    dCmn_dAlpha/Cm0 2*Xn(12) 0 0
             -dtn_dCST'./t0'    [0; 0]   [0 2*Xn(13) 0; 0 0 2*Xn(14)] ];

        cn = [Cln/Cl0-1; Cmn/Cm0-1+Xn(12)^2; 1-tn'./t0'+Xn(13:14).^2];
           
        Pn = Cdn - sum(lambda.*cn);
         
        if Pn <= P+mu*a*g'*S
            break;
        else
            a = a*0.5;
        end
    end
    
    % *************** BFGS *********************
    Qn = gn - sum(meshgrid(lambda,1:length(X))'.*An)';   
    Q  = g  - sum(meshgrid(lambda,1:length(X))'.*A)';   

    y = Qn - Q;

    if dn'*y >= 0.2*dn'*H*dn
        theta = 1;
    else
        theta = 0.8*dn'*H*dn/(dn'*H*dn-dn'*y);
    end

    eta_t = theta*y + (1-theta)*H*dn;
    H = H - H*(dn*dn')*H/(dn'*H*dn) + eta_t*eta_t'/(dn'*eta_t); 
    
    % ******************* Plotting **********************

    figure(1)
    subplot(2,2,1);
    hold on
    if iter==1
        plot(0,1, 'ko', 'MarkerFaceColor', [1 0 1]);
    end
    plot(iter,Cdn/Cd0, 'ko', 'MarkerFaceColor', [1 0 1]);
    title(num2str(Cdn/Cd0));
    xlabel('Iter');
    ylabel('F');

    subplot(2,2,2);
    hold on
    if iter==1
        plot(0,norm(c), 'ko', 'MarkerFaceColor', [1 0 1]);
    end
    plot(iter,norm(cn), 'ko', 'MarkerFaceColor', [1 0 1]);
    title(num2str(norm(cn)));
    xlabel('Iter');
    ylabel('C');

    subplot(2,2,3);
    hold on
    plot(iter,a, 'ko', 'MarkerFaceColor', [1 0 1]);
    title(num2str(a));
    xlabel('Iter');
    ylabel('Step lenght');
    
    subplot(2,2,4);
    hold on
    if iter==1
        plot(0,norm(Q), 'ko', 'MarkerFaceColor', [1 0 1]);
    end
    plot(iter,norm(Qn), 'ko', 'MarkerFaceColor', [1 0 1]);
    title(num2str(norm(Qn)));
    xlabel('Iter');
    ylabel('First Order Optimality');
    

    if round(iter/30) == iter/30
        display('Iter      nF      Objfun       norm of Const.        step length         First order optimality   ' );
    end   
    display([num2str(iter)  '         '  num2str(nf) '           '   num2str(Cdn/Cd0)  '            '  num2str(norm(cn)) '               '  num2str(a) '               '  num2str(norm(Qn))]);

    % ************ Updating *****************
     df = Cd-Cdn;
     X = Xn;  P = Pn;
     Cl = Cln; Cm = Cmn; Cd = Cdn; t = tn;
     dCl_dCST = dCln_dCST; dCd_dCST = dCdn_dCST; dCm_dCST = dCmn_dCST; dCl_dAlpha = dCln_dAlpha; dCd_dAlpha = dCdn_dAlpha; dCm_dAlpha = dCmn_dAlpha;
    
%     Sol_iter{iter} = Sol;
     X_iter(end+1,:) = X';
     
    % *********************  Convergence check **********************
    if norm(Qn) <= tol_kkt 
        disp('Firs order optimality lower than tolerance');
        break;
    elseif abs(df) <= tol_f
        disp('Cahnge in objective function is lower than tolerance');
        break;
    elseif s==s_max
        disp('Change in design variables is lower than tolerance');
        break;
    end
end

Xopt = X;
Fopt = Cdn/Cd0;


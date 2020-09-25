function [Xopt, Fopt] = RunOptimization(x0, ub, lb,options)

xLast = [];
f = [];
dfdx = [];
c = [];
dcdx = [];
ceq = [];
dceqdx = [];


fun = @objfun;
cfun = @nonlcons;


[Xopt, Fopt, exit] = fmincon(fun,x0,[],[],[],[],lb,ub,cfun,options);


    function [F, dFdX] = objfun(X)
        if  ~ isequal(X,xLast) || (nargout>1 && isempty(dfdx))
            if nargout > 1
                Adjoint = 1;
            else
                Adjoint = 0;
            end
            [f,dfdx, c, ceq, dcdx, dceqdx] = CallFEMflow(X,Adjoint);
            xLast = X;
        end
        F = f;
        dFdX = dfdx';
    end


    function [C, Ceq, dCdX, dCeqdX] = nonlcons(X)
        if  ~ isequal(X,xLast) || (nargout>2 && isempty(dcdx))
            if nargout > 2
                Adjoint = 1;
            else
                Adjoint = 0;
            end
            [f,dfdx, c, ceq, dcdx, dceqdx] = CallFEMflow(X,Adjoint);
            xLast = X;
        end
        C = c;
        Ceq = ceq;
        dCdX = dcdx';
        dCeqdX = dceqdx';
    end



end
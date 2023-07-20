function out = get_mu(Z,X,P,WN,ek,Nk,beta,fill)
%GET_MU Find the chemical potential with bisection method.
%   The bisection method is guaranteed to converge if DOS(w) >= 0.

%erreps = 1e-4;
erreps = 1e-8;

%set the initial values of the chemical potential.
muL = -10; %min(min(ek));
muR =  10; %max(max(ek));
muM = (muR+muL)/2;

fillL = get_filling(Z,X,P,ek,WN,Nk,beta,muL);
fillR = get_filling(Z,X,P,ek,WN,Nk,beta,muR);

%find mu by bisection method (linear convergence) ------------------------
fillM = get_filling(Z,X,P,ek,WN,Nk,beta,muM);
iter = 0;
loop = 1;

while (abs(fillM-fill)>erreps & loop == 1)
    iter = iter + 1;
    if (fillM > fill)
        fillR = fillM;
        muR = muM;
    elseif (fillM < fill)
        fillL = fillM;
        muL = muM;
    end
    muM = (muR+muL)/2;
    %get the new filling at the mid point.
    fillM = get_filling(Z,X,P,ek,WN,Nk,beta,muM);
    if (iter > 100)
        fprintf('\n')
        fprintf('  Cannot find the proper chemical potential!\n')
        fprintf('  Adopting mu = %g for filling = %g\n',muM,fillM)
        loop = 0;
    end
end

out = muM;
return

%find mu using Matlab fzero method (super-linear convergence) -------------
%fzero uses the algorithm that was originated by T. Dekker, which uses a
%combination of bisection, secant, and inverse quadratic interpolation
%methods. See Press et al., Numerical Recipes, Chapter 9.3 Van
%Wijngaarden-Dekker-Brent Method.
% options = optimset;
% options.TolX = erreps;
% out = fzero(@(xmu) get_filling(Z,X,P,ek,WN,Nk,beta,xmu)-fill,[muL muR],options);

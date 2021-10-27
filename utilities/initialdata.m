% set initial value in 2D
% Inputs:
% n: the number of points.
% bndcond: boundary condition, either 'Dirichlet' or 'Periodic'.
% ishomogeneous: if it is true, set it to be zero that we can go to
% inifinite accuracy.
% L_stencil: the corresponding 9 points stencil.
%
% Outputs: u,f: the given initial value.
% -------------------------------------------------------------
function [u,f] = initialdata(n,bndcond,ishomogeneous,L_stencil)

if nargin<=1
    bndcond = 'Dirichlet';
    ishomogeneous = true;
    L_stencil = [];
elseif nargin<=2
    ishomogeneous = true;
    L_stencil = [];
end

if strcmp(bndcond,'Dirichlet')
    u = zeros(n+1);
    rng(3,'v5uniform')
    u(2:end-1,2:end-1) = randn(n-1);
    if ishomogeneous
        f = zeros(n+1);
    else
        if isempty(L_stencil)
            error('The stencil needs to be given\n')
        else
            u_temp = u;
            rng(5,'v5uniform')
            u_temp(2:end-1,2:end-1) = randn(n-1);
            f = Lu(u_temp,L_stencil,1/(n-1),bndcond);
        end
    end
elseif strcmp(bndcond,'Periodic')
    rng(3,'v5uniform')
    u = randn(n);
    if ishomogeneous
        f = zeros(n);
    else
        if isempty(L_stencil)
            error('The stencil needs to be given\n')
        else
            rng(5,'v5uniform')
            u_temp = randn(n);
            f = Lu(u_temp,L_stencil,1/n,bndcond);
        end
    end
else
    error(['No such boundary condition - ' bndcond '\n'])
end

return

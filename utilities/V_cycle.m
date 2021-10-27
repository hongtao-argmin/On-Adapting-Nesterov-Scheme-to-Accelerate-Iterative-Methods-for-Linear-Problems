function u = V_cycle(u,f,h,nu1,nu2,w,L_stencil,bndcond,CycType,RelaxationMethod)    % Performs a V cycle; recursive.


for index = 1:nu1 % Perform nu1 relaxation sweeps
    switch RelaxationMethod
        case 'Jacobi'
            u = Jacobi(u,f,h,w,L_stencil,bndcond);
        case 'RB'
            u = RB(u,f,h,w,L_stencil,bndcond,false);
        case 'GS'
            u = GS(u,f,h,L_stencil,bndcond,w);
        case 'SymGS'
            u = SymmetricGS(u,f,h,L_stencil,bndcond,w);
        otherwise
            error(['No such Relaxation - ' RelaxationMethod])
    end
end
if strcmp(bndcond,'Dirichlet')
    if (length(u) <= 3)                        % Coarsest grid
        u = Jacobi(u,f,h,1,L_stencil,bndcond); % One undamped Jacobi solves exactly % because there is only one variable.
        return;
    end
elseif strcmp(bndcond,'Periodic')
    if (length(u) <= 1)                        % Coarsest grid
        for index = 1:1
            u = Jacobi(u,f,h,1,L_stencil,bndcond); % One undamped Jacobi solves exactly % because there is only one variable.
        end
        return;
    end
else
    error(['No such boundary condition - ' bndcond '\n'])
end

r = f - Lu(u,L_stencil,h,bndcond);
% Compute residual
fc = Restrict(r,bndcond);                       % Right-hand side of coarser grid
uc = zeros(size(fc));                           % Zero initial guess on coarse grid
%
for cycindex = 1:CycType
    uc = V_cycle(uc,fc,2*h,nu1,nu2,w,L_stencil,bndcond,CycType,RelaxationMethod);                % Recursive call
end
% A = ObtainExplicitMat(size(fc)-2,L_stencil,2*h,bndcond);
% uc(2:end-1,2:end-1) = reshape(A\reshape(fc(2:end-1,2:end-1),prod(size(fc)-2),1),size(fc)-2);
u = u + Prolong(uc,bndcond);                    % Interpolate and add correction

for index = 1:nu2 % Perform nu2 relaxation sweeps
    switch RelaxationMethod
        case 'Jacobi'
            u = Jacobi(u,f,h,w,L_stencil,bndcond);
        case 'RB'
            if nu1>0
                % // still uses false is faster than set true for pure MG // 
                % we may need to change the restriction.
                u = RB(u,f,h,w,L_stencil,bndcond,false); 
            else
                u = RB(u,f,h,w,L_stencil,bndcond,false);
            end
        case 'GS'
            u = GS(u,f,h,L_stencil,bndcond,w);
        case 'SymGS'
            u = SymmetricGS(u,f,h,L_stencil,bndcond,w);
        otherwise
            error(['No such Relaxation - ' RelaxationMethod])
    end
end

return

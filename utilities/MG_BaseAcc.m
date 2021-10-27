
function [u,u_set,rnorm,CPUTime] = MG_BaseAcc(L_stencil,n,nu1,nu2,n_cycles,CycType,w,bndcond,accpara,ishomogeneous,RelaxationMethod,StopTolerance)

h = 1/n;
[u,f] = initialdata(n,bndcond,ishomogeneous,L_stencil);
r = f - Lu(u,L_stencil,h,bndcond);                   % The initial residual.
rnorm = norm(r(:))/n;                                % This will store the properly normalized
%t_k = 1;
u_set{1} = u;
u_k = u;
y_k = u;
CPUTime = zeros(n_cycles+1,1);
for iters = 1:n_cycles
    StartTime = tic;
    u_k_1 = V_cycle(y_k,f,h,nu1,nu2,w,L_stencil,bndcond,CycType,RelaxationMethod);
    y_k = u_k_1 + accpara*(u_k_1-u_k);
    u_k = u_k_1;
    CPUTime(iters+1) = toc(StartTime);
    r = f - Lu(u_k,L_stencil,h,bndcond);
    rnorm = [rnorm, norm(r(:))/n];
    u_set{iters+1} = u_k;
    if rnorm(end)<StopTolerance
        break;
    end
end
u = u_k;
n_cycles = length(rnorm)-1;
CPUTime = cumsum(CPUTime);
CPUTime = CPUTime(1:n_cycles+1);
temp = rnorm(2:end) ./ rnorm(1:end-1);
fprintf('%f\n',geomean(temp(end-5:end)))

end

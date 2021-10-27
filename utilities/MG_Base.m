
function [u,u_set,rnorm,CPUTime] = MG_Base(L_stencil,n,nu1,nu2,n_cycles,CycType,w,bndcond,ishomogeneous,RelaxationMethod,StopTolerance)

h = 1/n;
[u,f] = initialdata(n,bndcond,ishomogeneous,L_stencil);
r = f - Lu(u,L_stencil,h,bndcond);                   % The initial residual.
rnorm = norm(r(:))/n;             % This will store the properly normalized

u_set{1} = u;
CPUTime = zeros(n_cycles+1,1);
for iters = 1:n_cycles
    StartTime = tic;
    u = V_cycle(u,f,h,nu1,nu2,w,L_stencil,bndcond,CycType,RelaxationMethod);
    CPUTime(iters+1) = toc(StartTime);
    r = f - Lu(u,L_stencil,h,bndcond); 
    rnorm = [rnorm, norm(r(:))/n];
    u_set{iters+1} = u;
    if rnorm(end)<StopTolerance
        break;
    end
end 
n_cycles = length(rnorm)-1;
CPUTime = cumsum(CPUTime);
CPUTime = CPUTime(1:n_cycles+1);
temp = rnorm(2:end) ./ rnorm(1:end-1);
fprintf('%f\n',geomean(temp(end-5:end)))
end

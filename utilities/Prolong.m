function u = Prolong(uc,bndcond)                  % Bi-linear interpolation
if strcmp(bndcond,'Dirichlet')
    u = zeros(length(uc)*2-1);
    u(1:2:end,1:2:end) = uc;
    u(2:2:end-1, 1:2:end) = 0.5*(u(1:2:end-2,1:2:end) + u(3:2:end,1:2:end));
    u(1:end,2:2:end-1) = 0.5*(u(1:end,1:2:end-2) + u(1:end,3:2:end)); 
elseif strcmp(bndcond,'Periodic')
    u = zeros(2*size(uc));
    u(1:2:end-1,1:2:end-1) = uc;
    % prolong except boundary; row
    u(2:2:end-2,1:2:end) = 0.5*(u(1:2:end-3,1:2:end) + u(3:2:end-1,1:2:end));
    % update boundary row
    u(end,1:2:end-1) = 0.5*(u(end-1,1:2:end-1) + u(1,1:2:end-1));
    % prolong except boundary; column
    u(1:end,2:2:end-2) = 0.5*(u(1:end,1:2:end-3) + u(1:end,3:2:end-1));
    u(1:end,end) = 0.5*(u(1:end,end-1) + u(1:end,1));
else
    error(['No such boundary condition - ' bndcond '\n'])
end

return
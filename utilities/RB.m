% implement RB in an efficient way later
function u = RB(u,f,h,w,L_stencil,bndcond,iseven)         % Red-Black Gauss-Seidel relaxation
% L_stencil is assumed to be 9 poitns.
if isempty(w)
    w = 1;
end
if nargin<=6
    iseven = false;
end

if strcmp(bndcond,'Dirichlet')
    [M,N] = size(u);
    mask = zeros(M,N);  % Create a Red-Black mask
    
    %         mask(1:2:end,1:2:end) = 1;
    %         mask(2:2:end,2:2:end) = 1;
    %
    %         % implement in more efficient way.
    %         u = u .* (1-mask) + mask .* Jacobi(u,f,h,w,L_stencil,bndcond);
    %         u = u.* mask + (1-mask) .* Jacobi(u,f,h,w,L_stencil,bndcond);
    
    mask(2:2:end-1,2:2:end-1) = 1;
    mask(3:2:end-2,3:2:end-2) = 1;
    mask_comp = zeros(M,N);
    
    % decide red-black or black-red
    if iseven
        mask_comp = mask;
        mask(2:end-1,2:end-1) = 1-mask_comp(2:end-1,2:end-1);
    else
        mask_comp(2:end-1,2:end-1) = 1-mask(2:end-1,2:end-1);
    end
    
    mask = logical(mask);
    mask_comp = logical(mask_comp);
    
    u(mask) = u(mask)+(h^2*w/L_stencil(2,2))*(f(mask) - (L_stencil(2,2)*u(mask)+L_stencil(1,2)*u(circshift(mask,-1))+...
        L_stencil(2,1)*u(circshift(mask,-1,2))+L_stencil(2,3)*u(circshift(mask,1,2))+...
        L_stencil(3,2)*u(circshift(mask,1))+L_stencil(1,1)*u(circshift(mask,[-1 -1]))+...
        L_stencil(1,3)*u(circshift(mask,[-1 1]))+L_stencil(3,1)*u(circshift(mask,[1 -1]))+...
        L_stencil(3,3)*u(circshift(mask,[1 1])))/h^2);
    
    u(mask_comp) = u(mask_comp)+(h^2*w/L_stencil(2,2))*(f(mask_comp) - (L_stencil(2,2)*u(mask_comp)+L_stencil(1,2)*u(circshift(mask_comp,-1))+...
        L_stencil(2,1)*u(circshift(mask_comp,-1,2))+L_stencil(2,3)*u(circshift(mask_comp,1,2))+...
        L_stencil(3,2)*u(circshift(mask_comp,1))+L_stencil(1,1)*u(circshift(mask_comp,[-1 -1]))+...
        L_stencil(1,3)*u(circshift(mask_comp,[-1 1]))+L_stencil(3,1)*u(circshift(mask_comp,[1 -1]))+...
        L_stencil(3,3)*u(circshift(mask_comp,[1 1])))/h^2);
    
elseif strcmp(bndcond,'Periodic')
    
    % periodic expending
    u_temp = [u(end,end) u(end,:) u(end,1); u(:,end) u u(:,1);u(1,end) u(1,:) u(1,1)];
    [M,N] = size(u_temp);
    mask = zeros(M,N);  % Create a Red-Black mask
    f_temp = zeros(M,N);
    f_temp(2:end-1,2:end-1) = f;
    mask(2:2:end-2,2:2:end-2) = 1;
    mask(3:2:end-1,3:2:end-1) = 1;
    mask_comp = zeros(M,N);
    
    if iseven
        mask_comp = mask;
        mask(2:end-1,2:end-1) = 1-mask_comp(2:end-1,2:end-1);
    else
        mask_comp(2:end-1,2:end-1) = 1-mask(2:end-1,2:end-1);
    end
    
    mask = logical(mask);
    mask_comp = logical(mask_comp);
    
    u_temp(mask) = u_temp(mask)+(h^2*w/L_stencil(2,2))*(f_temp(mask) - (L_stencil(2,2)*u_temp(mask)+L_stencil(1,2)*u_temp(circshift(mask,-1))+...
        L_stencil(2,1)*u_temp(circshift(mask,-1,2))+L_stencil(2,3)*u_temp(circshift(mask,1,2))+...
        L_stencil(3,2)*u_temp(circshift(mask,1))+L_stencil(1,1)*u_temp(circshift(mask,[-1 -1]))+...
        L_stencil(1,3)*u_temp(circshift(mask,[-1 1]))+L_stencil(3,1)*u_temp(circshift(mask,[1 -1]))+...
        L_stencil(3,3)*u_temp(circshift(mask,[1 1])))/h^2);
    % correct again boundary
    u_temp = [u_temp(end-1,end-1) u_temp(end-1,2:end-1) u_temp(end-1,2); u_temp(2:end-1,end-1) u_temp(2:end-1,2:end-1) u_temp(2:end-1,2);u_temp(2,end-1) u_temp(2,2:end-1) u_temp(2,2)];
    
    u_temp(mask_comp) = u_temp(mask_comp)+(h^2*w/L_stencil(2,2))*(f_temp(mask_comp) - (L_stencil(2,2)*u_temp(mask_comp)+L_stencil(1,2)*u_temp(circshift(mask_comp,-1))+...
        L_stencil(2,1)*u_temp(circshift(mask_comp,-1,2))+L_stencil(2,3)*u_temp(circshift(mask_comp,1,2))+...
        L_stencil(3,2)*u_temp(circshift(mask_comp,1))+L_stencil(1,1)*u_temp(circshift(mask_comp,[-1 -1]))+...
        L_stencil(1,3)*u_temp(circshift(mask_comp,[-1 1]))+L_stencil(3,1)*u_temp(circshift(mask_comp,[1 -1]))+...
        L_stencil(3,3)*u_temp(circshift(mask_comp,[1 1])))/h^2);
    
    u = u_temp(2:end-1,2:end-1);
    
else
    error(['No such boundary condition - ' bndcond '\n'])
end

end
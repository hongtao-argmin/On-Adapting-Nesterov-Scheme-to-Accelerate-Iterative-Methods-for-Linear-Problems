
% *** need to correct *** %
function u = SM_Gauss_Seidel(u,h,f,epsi,phi,nu,w)              % Jacobi relaxation with damping w
% zero boundary condition
if nargin<6
    w = 1;
end
C = cos(phi);
S = sin(phi);
u_temp = zeros(size(u)+2);
u_temp(2:end-1,2:end-1) = u;
u = u_temp;
for iters = 1:nu
    for j = 2:size(u,2)-1
        for i = 2:size(u,1)-1
            Lu_temp1 = -(...
                (epsi*C^2+S^2)*(u(i-1,j) + u(i+1,j)) + ...
                (C^2+epsi*S^2)*(u(i,j-1) + u(i,j+1)) - ...
                0.5*(1-epsi)*C*S*(u(i+1,j+1)-u(i-1,j+1)-...
                u(i+1,j-1)+u(i-1,j-1))...
                - 2*(1+epsi)*u(i,j)...
                )/h^2;
            u(i,j) = u(i,j) + w*(f(i,j) - Lu_temp1)*h^2/(2*(1+epsi)); % Assumes f = L(u) on boundaries
        end
    end
    %%%%%%% iterative from bottom
    % symmetric formulation
    for j = size(u,2)-1:-1:2
        for i = size(u,1)-1:-1:2
            Lu_temp1 = -(...
                (epsi*C^2+S^2)*(u(i-1,j) + u(i+1,j)) + ...
                (C^2+epsi*S^2)*(u(i,j-1) + u(i,j+1)) - ...
                0.5*(1-epsi)*C*S*(u(i+1,j+1)-u(i-1,j+1)-...
                u(i+1,j-1)+u(i-1,j-1))...
                - 2*(1+epsi)*u(i,j)...
                )/h^2;
            u(i,j) = u(i,j) + w*(f(i,j) - Lu_temp1)*h^2/(2*(1+epsi));
        end
    end
    u = u(2:end-1,2:end-1);
    % u = u + w*(f - L(u,h,p))./Denom_PDE(u,h,p);
end

return

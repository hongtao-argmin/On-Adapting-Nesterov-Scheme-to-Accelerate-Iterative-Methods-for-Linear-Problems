%
%
%
function u = GS(u,f,h,L_stencil,bndcond,w)
% Gauss-Seidel
% L_stencil is assumed to be 9 points.
if isempty(w)
    w = 1;
end
if strcmp(bndcond,'Dirichlet')
    [M,N] = size(u);
    
    for index_j = 2:N-1
        for index_i = 2:M-1
            
            u(index_i,index_j) = u(index_i,index_j) + w*(h^2/L_stencil(2,2))*(f(index_i,index_j) - ...
                (L_stencil(2,2)*u(index_i,index_j)+L_stencil(1,1)*u(index_i-1,index_j-1)+...
                L_stencil(1,2)*u(index_i-1,index_j)+L_stencil(1,3)*u(index_i-1,index_j+1)+L_stencil(2,1)*u(index_i,index_j-1)+...
                L_stencil(2,3)*u(index_i,index_j+1)+L_stencil(3,1)*u(index_i+1,index_j-1)+L_stencil(3,2)*u(index_i+1,index_j)+...
                L_stencil(3,3)*u(index_i+1,index_j+1))/h^2);
        end
    end
    
    
elseif strcmp(bndcond,'Periodic')
    
    
    
else
    error(['No such boundary condition - ' bndcond '\n'])
end

end
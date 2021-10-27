%
%
% b_1 and b_N: minimum and maximum eigenvalues of B which is the iteration
% matrix rather than A in our problem.
%
% c_optimal: output optimal parameters.
% r_star: predict accelerated convergence factor.
% Author: Tao Hong
% 2th Jan. 2020.
% ---------------------------------
function [c_opt,r_star] = CalculateOptimalFixedpara(b_1,b_N)

if abs(b_1)<=b_N
    % determine r
    if b_1>=-1/3*b_N
        c_opt = (2-2*sqrt(1-b_N))/b_N-1;
        r_star = 1-sqrt(1-b_N);
        r_star = real(r_star);
    elseif b_1>-b_N && b_1<0
        Omega = 2*b_1*b_N*(b_1+b_N);
        c_opt = (sqrt((4*Omega+(b_1-b_N)^2)*(b_1-b_N)^2)-(2*Omega+(b_1-b_N)^2))/(2*Omega);
        r_star = 0.5*((1+c_opt)*b_N+sqrt((1+c_opt)^2*b_N^2-4*c_opt*b_N));
        r_star = real(r_star);
    elseif b_1 == -b_N
        c_opt = 0;
        r_star = b_N;
        r_star = real(r_star);
        %     else
        %         b_N_new = -b_1;
        %         b_1_new = -b_N;
        %         if b_1_new>=-1/3*b_N_new
        %             c_opt = (2-2*sqrt(1-b_N_new))/b_N_new-1;
        %             r_star = 1-sqrt(1-b_N_new);
        %         elseif b_1_new>-b_N_new && b_1_new<0
        %             Omega = 2*b_1_new*b_N_new*(b_1_new+b_N_new);
        %             c_opt = (sqrt((4*Omega+(b_1_new-b_N_new)^2)*(b_1_new-b_N_new)^2)-(2*Omega+(b_1_new-b_N_new)^2))/(2*Omega);
        %             r_star = 0.5*((1+c_opt)*b_N_new+sqrt((1+c_opt)^2*b_N_new^2-4*c_opt*b_N_new));
        %         elseif b_1_new == -b_N_new
        %             c_opt = 0;
        %             r_star = b_N_new;
        %         end
        %         r_star = real(r_star);
    else
        error('Not define for such values of b\n')
    end
else
    if b_N<=-1/3*b_1
        c_opt = (2-2*sqrt(1-b_1))/b_1-1;
        if c_opt< -3+2*sqrt(2)
            c_opt = -3+2*sqrt(2);
        end
        %r_star = sqrt(1-b_1)-1;
        r_star = -0.5*((1+c_opt)*b_1-sqrt((1+c_opt)^2*b_1^2-4*c_opt*b_1));
        r_star = real(r_star);
    elseif  b_N>-1/3*b_1 && b_1 ~= -b_N
        Omega = 2*b_1*b_N*(b_1+b_N);
        c_opt = (sqrt((4*Omega+(b_1-b_N)^2)*(b_1-b_N)^2)-(2*Omega+(b_1-b_N)^2))/(2*Omega);
        if c_opt< -3+2*sqrt(2)
            c_opt = -3+2*sqrt(2);
        end
        r_star = 0.5*((1+c_opt)*b_N+sqrt((1+c_opt)^2*b_N^2-4*c_opt*b_N));
        r_star = real(r_star);
    elseif b_1 == -b_N
        c_opt = 0;
        r_star = b_N;
        r_star = real(r_star);
    else
        error('Not define for such values of b\n')
    end
    
end
fprintf('The optimal coffeicient and new accelerated Conv. Factor are %f,%f\n',c_opt,r_star)
return
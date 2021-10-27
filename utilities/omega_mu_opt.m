% This file is writen to calculate the optimal damping factor of a nine
% points stencil.
% Reference: Yavneh I, Olvovsky E. ``Multigrid smoothing for symmetric nine-point stencils''. Applied Mathematics and Computation, 1998, 92(2-3): 229-246.
% phi = 0;
% epsilon = 0.1;
% will be developed fully in the future.
% -------------------------------------------------------
function [omega_opt,mu_opt] = omega_mu_opt(phi,epsilon)

C = cos(phi);
S = sin(phi);
A = C^2+epsilon*S^2;
B = -(1-epsilon)*C*S;
C = epsilon*C^2+S^2;
% if A==C
%     isAC = true;
% else
%     isAC = false;
% end

D = 0;

Scal = A+C+D;
A = A/Scal;
C = C/Scal;
D = D/Scal;
disp(['Check Normalization ' num2str(A+C+D)])

% if (D+A)*(D+C)>B^2 && B==0 && A>=0&&B>=0
%     S_min = -1;
%     S_max = max(C,(D-C)/(1+A));
if (D+A)*(D+C)>B^2 && B==0
    S_min = -A-abs(D-C);
    S_max = max(A,abs(D-C)-A);
    disp('It is ellipticity')
elseif B~=0
    S_min = -1;
    S_max = sqrt(A^2+B^2);
end

omega_opt = 2/(2-(S_max+S_min));
mu_opt = (S_max-S_min)/(2-(S_max+S_min));

% fprintf('The optimal damping factor and mu is: %f, %f\n',...
%     omega_opt,mu_opt)

end

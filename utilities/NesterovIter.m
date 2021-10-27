% obtain the ACF for each symbol.
%
%

function r_final = NesterovIter(symbol,c_opt,b_1,b_N,isPrint)

if nargin<=4
    isPrint = false;
end
if isempty(c_opt)
    if isempty(b_1) || isempty(b_N)
        error('c is not given that the value of minimal and maximal eigenvalues must be given\n')
    else
        c_opt = CalculateOptimalFixedpara(b_1,b_N);
    end
end

r1 = 0.5*abs((1+c_opt)*symbol+sqrt((1+c_opt)^2*symbol^2-4*c_opt*symbol));
r2 = 0.5*abs((1+c_opt)*symbol-sqrt((1+c_opt)^2*symbol^2-4*c_opt*symbol));
r_final = max(r1,r2);

if isPrint
    fprintf('The final conv. factor for symbol %f%+fj is: %f \n',real(symbol),imag(symbol),r_final);
end

% Derive Equations - Variable grid derivatives
%
% This script derives the coefficients and expressions for computing
% the first and second finite differences for a function that is given on a
% non-uniform grid. It fits a quadratic between three points and then
% differentiates that analytically. If you set hLow==hUpp, then this will
% reduce to the standard equations for second-order accurate finite
% differences.
%
% tLow = time at lower grid point
% tMid = time at middle grid point
% tUpp = time at upper grid point
%
% hLow = tMid - tLow;
% hUpp = tUpp - tMid;
%
% xLow = x(tLow) = position at lower grid point
% xMid = x(tMid) = position at middle grid point
% xUpp = x(tUpp) = position at upper grid point
%


clc; clear;

flagCheckUniformGrid = false;
if flagCheckUniformGrid
    syms h 'real'
    hLow = h;
    hUpp = h;
else
    syms hLow hUpp 'real'
end

syms xLow xMid xUpp 'real'
syms t A B C 'real'

x = A*t^2 + B*t + C;
v = diff(x,t);
a = diff(v,t);

eqns(1) = subs(x,t,-hLow) - xLow;
eqns(2) = subs(x,t,0) - xMid;
eqns(3) = subs(x,t,hUpp) - xUpp;

vars = [A;B;C];

soln = solve(eqns,vars);

soln.A = simplify(soln.A);
soln.B = simplify(soln.B);
soln.C = simplify(soln.C);

vLow = simplify(subs(v,{'A','B','C','t'},{soln.A,soln.B,soln.C,-hLow}));
aLow = simplify(subs(a,{'A','B','C','t'},{soln.A,soln.B,soln.C,-hLow}));

vMid = simplify(subs(v,{'A','B','C','t'},{soln.A,soln.B,soln.C,0}));
aMid = simplify(subs(a,{'A','B','C','t'},{soln.A,soln.B,soln.C,0}));

vUpp = simplify(subs(v,{'A','B','C','t'},{soln.A,soln.B,soln.C,hUpp}));
aUpp = simplify(subs(a,{'A','B','C','t'},{soln.A,soln.B,soln.C,hUpp}));

disp(['vLow = ' vectorize(vLow) ';']);
disp(['vMid = ' vectorize(vMid) ';']);
disp(['vUpp = ' vectorize(vUpp) ';']);
disp(['aMid = ' vectorize(aMid) ';']);
% vLow = (hLow.^2.*xMid - hLow.^2.*xUpp - hUpp.^2.*xLow + hUpp.^2.*xMid - 2.*hLow.*hUpp.*xLow + 2.*hLow.*hUpp.*xMid)./(hLow.*hUpp.*(hLow + hUpp));
% vMid = -(hLow.^2.*xMid - hLow.^2.*xUpp + hUpp.^2.*xLow - hUpp.^2.*xMid)./(hLow.*hUpp.*(hLow + hUpp));
% vUpp = -(hLow.^2.*xMid - hLow.^2.*xUpp - hUpp.^2.*xLow + hUpp.^2.*xMid + 2.*hLow.*hUpp.*xMid - 2.*hLow.*hUpp.*xUpp)./(hLow.*hUpp.*(hLow + hUpp));
% aMid = -(2.*(hLow.*xMid - hLow.*xUpp - hUpp.*xLow + hUpp.*xMid))./(hLow.*hUpp.*(hLow + hUpp));


%%%% Gradients:
dAcc_dhLow = simplify(diff(aMid,hLow));
dAcc_dhUpp = simplify(diff(aMid,hUpp));
disp(['dAcc_dhLow = ' vectorize(dAcc_dhLow) ';']);
disp(['dAcc_dhUpp = ' vectorize(dAcc_dhUpp) ';']);
% dAcc_dhLow = (2.*(hLow.^2.*xMid - hLow.^2.*xUpp - hUpp.^2.*xLow + hUpp.^2.*xMid - 2.*hLow.*hUpp.*xLow + 2.*hLow.*hUpp.*xMid))./(hLow.^2.*hUpp.*(hLow + hUpp).^2);
% dAcc_dhUpp = (2.*(hLow.^2.*xMid - hLow.^2.*xUpp - hUpp.^2.*xLow + hUpp.^2.*xMid + 2.*hLow.*hUpp.*xMid - 2.*hLow.*hUpp.*xUpp))./(hLow.*hUpp.^2.*(hLow + hUpp).^2);


function out = tutorial2
out{1} = @init;
out{2} = @fun_eval;
out{3} = [];
out{4} = [];
out{5} = [];
out{6} = [];
out{7} = [];
out{8} = [];
out{9} = [];

% --------------------------------------------------------------------------
function dydt = fun_eval(t,kmrgd,par_r,par_a,par_b,par_c,par_d)
dydt=[par_r*kmrgd(1)*(1-kmrgd(1))-(kmrgd(1)*kmrgd(2))/(kmrgd(1)+par_a);
-par_c*kmrgd(2)+(kmrgd(1)*kmrgd(2))/(kmrgd(1)+par_a)-(par_d*kmrgd(2)^2)/(kmrgd(2)^2+par_b^2);];

% --------------------------------------------------------------------------
function [tspan,y0,options] = init
handles = feval(tutorial2);
y0=[0,0];
options = odeset('Jacobian',[],'JacobianP',[],'Hessians',[],'HessiansP',[]);
tspan = [0 10];

% --------------------------------------------------------------------------
function jac = jacobian(t,kmrgd,par_r,par_a,par_b,par_c,par_d)
% --------------------------------------------------------------------------
function jacp = jacobianp(t,kmrgd,par_r,par_a,par_b,par_c,par_d)
% --------------------------------------------------------------------------
function hess = hessians(t,kmrgd,par_r,par_a,par_b,par_c,par_d)
% --------------------------------------------------------------------------
function hessp = hessiansp(t,kmrgd,par_r,par_a,par_b,par_c,par_d)
%---------------------------------------------------------------------------
function tens3  = der3(t,kmrgd,par_r,par_a,par_b,par_c,par_d)
%---------------------------------------------------------------------------
function tens4  = der4(t,kmrgd,par_r,par_a,par_b,par_c,par_d)
%---------------------------------------------------------------------------
function tens5  = der5(t,kmrgd,par_r,par_a,par_b,par_c,par_d)

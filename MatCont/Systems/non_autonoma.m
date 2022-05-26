function out = non_autonoma
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
function dydt = fun_eval(t,kmrgd,par_a)
dydt=[par_a*sin(t)*kmrgd(1);];

% --------------------------------------------------------------------------
function [tspan,y0,options] = init
handles = feval(non_autonoma);
y0=[0];
options = odeset('Jacobian',[],'JacobianP',[],'Hessians',[],'HessiansP',[]);
tspan = [0 10];

% --------------------------------------------------------------------------
function jac = jacobian(t,kmrgd,par_a)
% --------------------------------------------------------------------------
function jacp = jacobianp(t,kmrgd,par_a)
% --------------------------------------------------------------------------
function hess = hessians(t,kmrgd,par_a)
% --------------------------------------------------------------------------
function hessp = hessiansp(t,kmrgd,par_a)
%---------------------------------------------------------------------------
function tens3  = der3(t,kmrgd,par_a)
%---------------------------------------------------------------------------
function tens4  = der4(t,kmrgd,par_a)
%---------------------------------------------------------------------------
function tens5  = der5(t,kmrgd,par_a)
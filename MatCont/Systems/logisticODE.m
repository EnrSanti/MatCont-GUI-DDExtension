function out = logisticODE
out{1} = @init;
out{2} = @fun_eval;
out{3} = [];
out{4} = [];
out{5} = [];
out{6} = [];
out{7} = [];
out{8} = [];
out{9} = [];
out{10}= @RR3;
out{11}= @r2;
out{12}= @r4;
out{13}= @c;

% --------------------------------------------------------------------------
function dydt = fun_eval(t,kmrgd,par_r)
dydt=[par_r*kmrgd(1)*(1-kmrgd(1));];

% --------------------------------------------------------------------------
function [tspan,y0,options] = init
handles = feval(logisticODE);
tspan = [0 10];

% --------------------------------------------------------------------------
function jac = jacobian(t,kmrgd,par_r)
% --------------------------------------------------------------------------
function jacp = jacobianp(t,kmrgd,par_r)
% --------------------------------------------------------------------------
function hess = hessians(t,kmrgd,par_r)
% --------------------------------------------------------------------------
function hessp = hessiansp(t,kmrgd,par_r)
%---------------------------------------------------------------------------
function tens3  = der3(t,kmrgd,par_r)
%---------------------------------------------------------------------------
function tens4  = der4(t,kmrgd,par_r)
%---------------------------------------------------------------------------
function tens5  = der5(t,kmrgd,par_r)
function userfun1=RR3(t,kmrgd,par_r)
	userfun1=0;
function userfun2=r2(t,kmrgd,par_r)
	userfun2=0;
function userfun3=r4(t,kmrgd,par_r)
	userfun3=0;
function userfun4=c(t,kmrgd,par_r)
	userfun4=0;

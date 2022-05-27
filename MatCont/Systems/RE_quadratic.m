function [out,rhs] = RE_quadratic
out{1} = @init;
out{2} = @fun_eval;
out{3} = [];
out{4} = [];
out{5} = [];
out{6} = [];
out{7} = [];
out{8} = [];
out{9} = [];
rhs{1}=@RHSre1;

% --------------------------------------------------------------------------


function dydt = fun_eval (t,state,par_gamma)
[thetaCap,wCap]=fclencurtVals();
M=10;
d1=1;
d2=0;
delayFunctions=[-3,-1];
tau_max=max(abs(delayFunctions));
ScaledNodes=UnitNodesFun()*tau_max;
ScaledDD=UnitDDFun()/tau_max;
BaryWeights=BaryWeightsFun();
yM=state(1:d2);
VM=state(d2+1:(M+1)*d2);
UM=state((d2*M+d2+1):end);
derState=kron(ScaledDD(2:end,2:end),eye(d1))*UM; %DM*state
dMDM_DDE=kron(ScaledDD(2:end,:),eye(d2));
KM = derState - kron([par_gamma./2.*(dot(commonFunctions.interpoly(+thetaCap*(-3-(-1))+-1,ScaledNodes,[0;derState(1:d1:end)],BaryWeights).*(1-commonFunctions.interpoly(+thetaCap*(-3-(-1))+-1,ScaledNodes,[0;derState(1:d1:end)],BaryWeights)),wCap)*(-1-(-3)))],ones(M,1));
dydt= [KM];

% --------------------------------------------------------------------------
function state_eq=init(M,xeq,yeq)
state_eq=[kron(ones(M,1),xeq); kron(ones(M+1,1),yeq)];

% --------------------------------------------------------------------------
function jac = jacobian(t,kmrgd,par_gamma)
% --------------------------------------------------------------------------
function jacp = jacobianp(t,kmrgd,par_gamma)
% --------------------------------------------------------------------------
function hess = hessians(t,kmrgd,par_gamma)
% --------------------------------------------------------------------------
function hessp = hessiansp(t,kmrgd,par_gamma)
%---------------------------------------------------------------------------
function tens3  = der3(t,kmrgd,par_gamma)
%---------------------------------------------------------------------------
function tens4  = der4(t,kmrgd,par_gamma)
%---------------------------------------------------------------------------
function tens5  = der5(t,kmrgd,par_gamma)

function out = UnitNodesFun
out=[0;-0.024472;-0.095492;-0.20611;-0.34549;-0.5;-0.65451;-0.79389;-0.90451;-0.97553;-1];
function out = UnitDDFun
out=[67,-81.7269,20.9443,-9.7037,5.7889,-4,3.0557,-2.5192,2.2111,-2.0502,1;20.4317,-9.9596,-14.0806,5.5055,-3.1151,2.1029,-1.5872,1.2997,-1.1363,1.0515,-0.51254;-5.2361,14.0806,-2.3416,-9.0403,4,-2.4721,1.7889,-1.4318,1.2361,-1.1363,0.55279;2.4259,-5.5055,9.0403,-0.89806,-7.1744,3.4026,-2.2301,1.7013,-1.4318,1.2997,-0.62981;-1.4472,3.1151,-4,7.1744,-0.34164,-6.4721,3.2361,-2.2301,1.7889,-1.5872,0.76393;1,-2.1029,2.4721,-3.4026,6.4721,-4.885e-15,-6.4721,3.4026,-2.4721,2.1029,-1;-0.76393,1.5872,-1.7889,2.2301,-3.2361,6.4721,0.34164,-7.1744,4,-3.1151,1.4472;0.62981,-1.2997,1.4318,-1.7013,2.2301,-3.4026,7.1744,0.89806,-9.0403,5.5055,-2.4259;-0.55279,1.1363,-1.2361,1.4318,-1.7889,2.4721,-4,9.0403,2.3416,-14.0806,5.2361;0.51254,-1.0515,1.1363,-1.2997,1.5872,-2.1029,3.1151,-5.5055,14.0806,9.9596,-20.4317;-1,2.0502,-2.2111,2.5192,-3.0557,4,-5.7889,9.7037,-20.9443,81.7269,-67];
function out = BaryWeightsFun
out=[0.5,-1,1,-1,1,-1,1,-1,1,-1,0.5];


function [thetaCap,wCap] = fclencurtVals
thetaCap=[1;0.97553;0.90451;0.79389;0.65451;0.5;0.34549;0.20611;0.095492;0.024472;0];
wCap=[0.0050505;0.04729;0.092818;0.12679;0.14961;0.15688;0.14961;0.12679;0.092818;0.04729;0.0050505];

function out = RHSre1(t,state,par_gamma)
[thetaCap,wCap]=fclencurtVals();
M=10;
d1=1;
d2=0;
delayFunctions=[-3,-1];tau_max=max(abs(delayFunctions));
ScaledNodes=UnitNodesFun()*tau_max;
ScaledDD=UnitDDFun()/tau_max;
BaryWeights=BaryWeightsFun();
yM=state(1:d2);
VM=state(d2+1:(M+1)*d2);
UM=state((d2*M+d2+1):d2*(M+1)+d1*M);
derState=kron(ScaledDD(2:end,2:end),eye(d1))*UM; %DM*state
out=par_gamma./2.*(dot(commonFunctions.interpoly(+thetaCap*(-3-(-1))+-1,ScaledNodes,[0;derState(1:d1:end)],BaryWeights).*(1-commonFunctions.interpoly(+thetaCap*(-3-(-1))+-1,ScaledNodes,[0;derState(1:d1:end)],BaryWeights)),wCap)*(-1-(-3)));



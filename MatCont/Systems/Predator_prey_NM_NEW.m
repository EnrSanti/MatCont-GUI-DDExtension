function out = Predator_prey_NM_NEW
out{1} = @init;
out{2} = @fun_eval;
out{3} = [];
out{4} = [];
out{5} = [];
out{6} = [];
out{7} = [];
out{8} = [];
out{9} = [];
out{10}= @D0;

% --------------------------------------------------------------------------


function dydt = fun_eval (t,state,par_juv_death,par_conv,par_death,par_pred,par_TAU)
M=10;
UnitQuadweights=UnitQuadweightsFun();
UnitNodes=UnitNodesFun();
UnitDD=UnitDDFun();
BaryWeights=BaryWeightsFun();
d1=0;
d2=2;
delayFunctions=[-par_TAU,-par_TAU];
tau_max=abs(min(delayFunctions));
yM=state((d1*M+1):(d1*M+d2));
VM=state((d1*M+d2+1):end);
GM = @(x) [
yM(1)*(1-yM(1))-par_pred*yM(1)*yM(2)                                                           ;
par_conv*par_pred*exp(-par_juv_death*par_TAU)*commonFunctions.interpoly(-par_TAU,tau_max*UnitNodes,[yM(1);VM(1:d2:end)],BaryWeights)*commonFunctions.interpoly(-par_TAU,tau_max*UnitNodes,[yM(2);VM(2:d2:end)],BaryWeights)-par_death*yM(2)];
KM=[];
dMDM_DDE=kron(UnitDD(2:end,:),eye(d2));
dydt= [GM(KM);(1/tau_max*dMDM_DDE)*[yM;VM]];

% --------------------------------------------------------------------------

function state_eq=init(M,xeq,yeq)
state_eq=[kron(ones(M,1),xeq); kron(ones(M+1,1),yeq)];

% --------------------------------------------------------------------------
function jac = jacobian(t,kmrgd,par_juv_death,par_conv,par_death,par_pred,par_TAU)
% --------------------------------------------------------------------------
function jacp = jacobianp(t,kmrgd,par_juv_death,par_conv,par_death,par_pred,par_TAU)
% --------------------------------------------------------------------------
function hess = hessians(t,kmrgd,par_juv_death,par_conv,par_death,par_pred,par_TAU)
% --------------------------------------------------------------------------
function hessp = hessiansp(t,kmrgd,par_juv_death,par_conv,par_death,par_pred,par_TAU)
%---------------------------------------------------------------------------
function tens3  = der3(t,kmrgd,par_juv_death,par_conv,par_death,par_pred,par_TAU)
%---------------------------------------------------------------------------
function tens4  = der4(t,kmrgd,par_juv_death,par_conv,par_death,par_pred,par_TAU)
%---------------------------------------------------------------------------
function tens5  = der5(t,kmrgd,par_juv_death,par_conv,par_death,par_pred,par_TAU)

function out = UnitQuadweightsFun
out=[0.0050505,0.04729,0.092818,0.12679,0.14961,0.15688,0.14961,0.12679,0.092818,0.04729,0.0050505];
function out = UnitNodesFun
out=[0;-0.024472;-0.095492;-0.20611;-0.34549;-0.5;-0.65451;-0.79389;-0.90451;-0.97553;-1];
function out = UnitDDFun
out=[67,-81.7269,20.9443,-9.7037,5.7889,-4,3.0557,-2.5192,2.2111,-2.0502,1;20.4317,-9.9596,-14.0806,5.5055,-3.1151,2.1029,-1.5872,1.2997,-1.1363,1.0515,-0.51254;-5.2361,14.0806,-2.3416,-9.0403,4,-2.4721,1.7889,-1.4318,1.2361,-1.1363,0.55279;2.4259,-5.5055,9.0403,-0.89806,-7.1744,3.4026,-2.2301,1.7013,-1.4318,1.2997,-0.62981;-1.4472,3.1151,-4,7.1744,-0.34164,-6.4721,3.2361,-2.2301,1.7889,-1.5872,0.76393;1,-2.1029,2.4721,-3.4026,6.4721,-4.885e-15,-6.4721,3.4026,-2.4721,2.1029,-1;-0.76393,1.5872,-1.7889,2.2301,-3.2361,6.4721,0.34164,-7.1744,4,-3.1151,1.4472;0.62981,-1.2997,1.4318,-1.7013,2.2301,-3.4026,7.1744,0.89806,-9.0403,5.5055,-2.4259;-0.55279,1.1363,-1.2361,1.4318,-1.7889,2.4721,-4,9.0403,2.3416,-14.0806,5.2361;0.51254,-1.0515,1.1363,-1.2997,1.5872,-2.1029,3.1151,-5.5055,14.0806,9.9596,-20.4317;-1,2.0502,-2.2111,2.5192,-3.0557,4,-5.7889,9.7037,-20.9443,81.7269,-67];
function out = BaryWeightsFun
out=[26214.4,-52428.8,52428.8,-52428.8,52428.8,-52428.8,52428.8,-52428.8,52428.8,-52428.8,26214.4];

function userfun1=D0(t,kmrgd,par_juv_death,par_conv,par_death,par_pred,par_TAU)
	userfun1=(-par_TAU)*(-par_TAU);

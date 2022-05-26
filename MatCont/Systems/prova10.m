function [out,rhs] = prova10
out{1} = @init;
out{2} = @fun_eval;
out{3} = [];
out{4} = [];
out{5} = [];
out{6} = [];
out{7} = [];
out{8} = [];
out{9} = [];
out{10}= @x;
rhs{1}=@RHSre1;

% --------------------------------------------------------------------------


function dydt = fun_eval (t,state,par_a)
[thetaCap,wCap]=fclencurtVals();
M=12;
d1=1;
d2=1;
delayFunctions=[-3,-3,-2];
tau_max=max(abs(delayFunctions));
ScaledNodes=UnitNodesFun()*tau_max;
ScaledDD=UnitDDFun()/tau_max;
BaryWeights=BaryWeightsFun();
yM=state(1:d2);
VM=state(d2+1:(M+1)*d2);
UM=state((d2*M+d2+1):end);
derState=kron(ScaledDD(2:end,2:end),eye(d1))*UM; %DM*state
GM = yM(1)+commonFunctions.interpoly(-3,ScaledNodes,[yM(1);VM(1:d2:end)],BaryWeights);
dMDM_DDE=kron(ScaledDD(2:end,:),eye(d2));
KM = derState - kron([dot(commonFunctions.interpoly(+thetaCap*(-2-(-3))+-3,ScaledNodes,[0;derState(1:d1:end)],BaryWeights)+commonFunctions.interpoly(+thetaCap*(-2-(-3))+-3,ScaledNodes,[yM(1);VM(1:d2:end)],BaryWeights),wCap)*(-2-(-3))],ones(M,1));
dydt= [GM;(1/tau_max*dMDM_DDE)*[yM;VM];KM];

% --------------------------------------------------------------------------
function state_eq=init(M,xeq,yeq)
state_eq=[kron(ones(M,1),xeq); kron(ones(M+1,1),yeq)];

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

function out = UnitNodesFun
out=[0;-0.017037;-0.066987;-0.14645;-0.25;-0.37059;-0.5;-0.62941;-0.75;-0.85355;-0.93301;-0.98296;-1];
function out = UnitDDFun
out=[96.3333,-117.391,29.8564,-13.6569,8,-5.3968,4,-3.1776,2.6667,-2.3431,2.1436,-2.0347,1;29.3477,-14.4195,-20.0199,7.7274,-4.2925,2.8284,-2.0706,1.633,-1.3643,1.1954,-1.0917,1.0353,-0.50867;-7.4641,20.0199,-3.4641,-12.5851,5.4641,-3.2938,2.3094,-1.778,1.4641,-1.2713,1.1547,-1.0917,0.5359;3.4142,-7.7274,12.5851,-1.4142,-9.6569,4.4614,-2.8284,2.0706,-1.6569,1.4142,-1.2713,1.1954,-0.58579;-2,4.2925,-5.4641,9.6569,-0.66667,-8.2925,4,-2.6357,2,-1.6569,1.4641,-1.3643,0.66667;1.3492,-2.8284,3.2938,-4.4614,8.2925,-0.2774,-7.7274,3.8637,-2.6357,2.0706,-1.778,1.633,-0.7944;-1,2.0706,-2.3094,2.8284,-4,7.7274,-6.6613e-15,-7.7274,4,-2.8284,2.3094,-2.0706,1;0.7944,-1.633,1.778,-2.0706,2.6357,-3.8637,7.7274,0.2774,-8.2925,4.4614,-3.2938,2.8284,-1.3492;-0.66667,1.3643,-1.4641,1.6569,-2,2.6357,-4,8.2925,0.66667,-9.6569,5.4641,-4.2925,2;0.58579,-1.1954,1.2713,-1.4142,1.6569,-2.0706,2.8284,-4.4614,9.6569,1.4142,-12.5851,7.7274,-3.4142;-0.5359,1.0917,-1.1547,1.2713,-1.4641,1.778,-2.3094,3.2938,-5.4641,12.5851,3.4641,-20.0199,7.4641;0.50867,-1.0353,1.0917,-1.1954,1.3643,-1.633,2.0706,-2.8284,4.2925,-7.7274,20.0199,14.4195,-29.3477;-1,2.0347,-2.1436,2.3431,-2.6667,3.1776,-4,5.3968,-8,13.6569,-29.8564,117.391,-96.3333];
function out = BaryWeightsFun
out=[349525.3333,-699050.6667,699050.6667,-699050.6667,699050.6667,-699050.6667,699050.6667,-699050.6667,699050.6667,-699050.6667,699050.6667,-699050.6667,349525.3333];


function [thetaCap,wCap] = fclencurtVals
thetaCap=[1;0.98296;0.93301;0.85355;0.75;0.62941;0.5;0.37059;0.25;0.14645;0.066987;0.017037;0];
wCap=[0.0034965;0.033029;0.065771;0.092382;0.11349;0.12634;0.13099;0.12634;0.11349;0.092382;0.065771;0.033029;0.0034965];

function out = RHSre1(t,state,par_a)
[thetaCap,wCap]=fclencurtVals();
M=12;
d1=1;
d2=1;
delayFunctions=[-3,-3,-2];tau_max=max(abs(delayFunctions));
ScaledNodes=UnitNodesFun()*tau_max;
ScaledDD=UnitDDFun()/tau_max;
BaryWeights=BaryWeightsFun();
yM=state(1:d2);
VM=state(d2+1:(M+1)*d2);
UM=state((d2*M+d2+1):d2*(M+1)+d1*M);
derState=kron(ScaledDD(2:end,2:end),eye(d1))*UM; %DM*state
out=dot(commonFunctions.interpoly(+thetaCap*(-2-(-3))+-3,ScaledNodes,[0;derState(1:d1:end)],BaryWeights)+commonFunctions.interpoly(+thetaCap*(-2-(-3))+-3,ScaledNodes,[yM(1);VM(1:d2:end)],BaryWeights),wCap)*(-2-(-3));


function userfun1=x(t,kmrgd,par_a)
	userfun1=kmrgd(1);

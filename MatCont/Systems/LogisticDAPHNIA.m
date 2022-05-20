function out = LogisticDAPHNIA
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


function dydt = fun_eval (t,state,par_beta,par_r,par_a,par_amax,par_k,par_gamma)
M=10;
[thetaCap,wCap]=fclencurt(20+1,0,1);
UnitQuadweights=UnitQuadweightsFun();
UnitNodes=UnitNodesFun();
UnitDD=UnitDDFun();
BaryWeights=BaryWeightsFun();
d1=1;
d2=1;
delayFunctions=[par_a,par_amax,par_a,par_amax];
tau_max=abs(min(delayFunctions));
yM=state(1:d2);
VM=state(d2+1:(M+1)*d2);
UM=state((d2*M+d2+1):end);
derState=kron(UnitDD(2:end,2:end),eye(d1))*UM; %DM*state
GM = 
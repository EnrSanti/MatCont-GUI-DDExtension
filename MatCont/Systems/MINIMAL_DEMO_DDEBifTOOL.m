function out = MINIMAL_DEMO_DDEBifTOOL
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


function dydt = fun_eval (t,state,par_TAU,par_a,par_b,par_d)
M=10;
UnitQuadweights=UnitQuadweightsFun();
UnitNodes=UnitNodesFun();
UnitDD=UnitDDFun();
BaryWeights=BaryWeightsFun();
d1=0;
d2=1;

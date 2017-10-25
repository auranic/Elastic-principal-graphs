function [TF_dynamics,t] = TF1_inhibits_TF2(parameters,tspan,initial_conditions)
%% Example of a function generating TF dynamics
%
%

[t,TF_dynamics] = ode23(@(t,y) func(t,y,parameters), tspan, initial_conditions);



end

function dydt = func(t,y,parameters)

ks1 = parameters(1,1);
ks2 = parameters(2,1);
kd1 = parameters(1,2);
kd2 = parameters(2,2);
A = parameters(2,3);
B = parameters(2,4);

dydt = [ks1-kd1*y(1);ks2*A/(B+y(1))-kd2*y(2)];

end


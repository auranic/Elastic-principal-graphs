function [TF_dynamics,t] = TF1_mutualinhibition_TF2_cofactorTF3(parameters,tspan,initial_conditions)
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
A1 = parameters(1,3);
B1 = parameters(1,4);
A2 = parameters(2,3);
B2 = parameters(2,4);

ks3 = parameters(3,1);
kd3 = parameters(3,2);
A3 = parameters(3,3);
B3 = parameters(3,4);

dydt = [-kd1*y(1)+A1*ks1/(B1+y(2))*y(3);-kd2*y(2)+ks2*A2/(B2+y(1))*y(3);-kd3*y(3)+ks3];

end




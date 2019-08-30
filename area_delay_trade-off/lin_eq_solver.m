%{
This is written by Zhiyang Ong and Andrew Mattheisen
	for EE 577B, SRAM Project - Part 2

Solver for simultaneous equations...
%}

%{
N*((0.5*Cout)/(Cin*beta))^(1/N)+N = 0
(N-1)*((0.5*Cout)/(Cin*(1-beta)))^(1/(N-1))+(N-1) = 0



syms x y alpha
[x,y] = solve(x^2*y^2, x-y/2-alpha)


eqs1 = 'x^2*y^2=1, x-y/2-alpha'
[x,y] = solve(eqs1)
%}

syms beta N Cout Cin
N=3;
%Cout=100;
%Cin=10;
[beta] = solve(N*((0.5*Cout)/(Cin*beta))^(1/N)+N, (N-1)*((0.5*Cout)/(Cin*(1-beta)))^(1/(N-1))+(N-1))

eqn='N*((0.5*Cout)/(Cin*beta))^(1/N)+N, (N-1)*((0.5*Cout)/(Cin*(1-beta)))^(1/(N-1))+(N-1)'
[b] = solve(eqn)

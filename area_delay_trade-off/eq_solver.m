%{
This is written by Zhiyang Ong (zhiyango@usc.edu;
6004 9194 12) and Andrew Mattheisen (mattheis@usc.edu;
2134 5147 11) for EE 577B, SRAM Project 
- Part 2

Solver for simultaneous equations...
%}
format long g
format compact
Cin=10%60;
Cout = 100%2176;
N=3;
j=1;
%for i=1:5;
    step=1e-3;
    beta=step;
    while beta < 1;
        % Left Hand Side (LHS) of the simultaneous equation
        a=N*((0.5*Cout)/(Cin*beta))^(1/N)+N
        % Left Hand Side (LHS) of the simultaneous equation
        b=(N-1)*((0.5*Cout)/(Cin*(1-beta)))^(1/(N-1))+(N-1)
        % Determine the difference between the LHS and RHS
        diff(j)=abs(a-b)
        b(j)=beta;
        beta=beta+step
        j=j+1;
    end
%end


%{
for i=1:length(b)
    if b(i) == 0
       b(i)=1e9 
       diff(i)=1e9
    end
end
%}


[min_diff index_md] = min(diff)
final_beta=index_md*step
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matlab file created by Andrew Mattheisen and Zhiyang Ong        %
% You can contact us at amattheisen@yahoo.com or zhiyang@ieee.org %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This matlab file determines the delay of multiple row decoders, %
% all using 2-4 predecoders and AND trees.                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp ('THIS IS OPTION 1')
%
format long g
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GENERIC PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m=10^-3;
u=10^-6;
n=10^-9;
p=10^-12;
f=10^-15;
a=10^-18;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TECHNOLOGY SPECIFIC PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lambda=.1*u;
CGate=2.396*f/u; % This is for an inverter gate, assumed to work for 
                 %   the NMOSgate, from TA's solution to part 1
CNDiff=2.244*f/u;   % nmos drain or source capacitance from TA's solution
                    % to part 1 
CPDiff=1.679*f/u;   % pmos drain or source capacitance from TA's solution
                    % to part 1 
REffPmos = 7.8648*10^3*u; % from TA's solution to part 1, problem 1e
REffNmos = 2.5832*10^3*u; % from TA's solution to part 1, problem 1e
rsheetm2 = .08; % this is in ohms per square from tech file
pinv=.5365*10^-10; % as calculated by us
tao=.1434*10^-10;  % as calculated by us
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONSTANTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gate capacitance of the inverter in the DFFPOSX1 cell connected to the
% output
Cin=4*u*CPDiff+2*u*CNDiff; % this is 1.1204e-14
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% CIRCUIT SPECIFIC PARAMETERS
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%!FILL IN SRAM Macro dimentions
NumCols=64;
NumRows=512;
%!
%!FILL IN LOGICAL EFFORT OF PATH
G=9/4; 
%!
%!FILL IN NUMBER OF STAGES
N=6;
%!
%!FILL IN PARACITIC CAPACITANCE
P=10;
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fixed Equations 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SideLoadCap=(4*lambda)*(40*lambda)*NumRows*14*a/(u*u) ...    % wire area cap
+ 2*(40*lambda)*NumRows*35*a/u;                          % + wire fringe cap
%
WordLineCap=(4*lambda)*(40*lambda)*NumCols*14*a/(u*u) ...   % wire area cap
+ 2*(40*lambda)*NumCols*35*a/u ... %                      + wire fringe cap
+ 2*CGate*NumCols*4*lambda; %                          + SRAM pass gate cap
%
% ELECTRICAL EFFORT OF PATH
H=WordLineCap/Cin; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% GUESSING - ITERATION
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%!FILL IN GUESS FOR GATE CAPACITANCE
%!Cgate_5=#####*f;
Cgate_5=6*f;
%!
%!FILL IN FANOUT OF BRANCH
%!B=#####*(#####*Cgate_5+SideLoadCap)/Cgate_5;
B=5*4*(64*Cgate_5+SideLoadCap)/Cgate_5;
%!
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% More Fixed Equations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PATH EFFORT
F=G*H*B;
%
%Should determine optimum N first, then set stage effort for that N.
%Nhat=round(log(F)/log(4)); %THIS IS COMMENTED OUT FOR DEBUGGING
Nhat=N;
%
% OPTIMAL STAGE EFFORT
fhat=F^(1/Nhat);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% Calculate Actual input Capacitances (Cin) for each stage
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% Cin = (Cout)*(g_i)/(fhat for each stage)
Cgate_actual(6)=WordLineCap*(1)/(fhat^(1+Nhat-N));
%!FILL IN STAGE EFFORT FROM LOAD TO CGATE_5
%!Cgate_5_actual=WordLineCap*(######)/(fhat^3); 
%!Use this equation: optimized using Nhat
%!   Cgate_actual(5)=WordLineCap*(6/4)/(fhat^(2+Nhat-N)); 
Cgate_actual(5)=Cgate_actual(6)*(6/4)/(fhat);
%!check to see if this matches the guess
Cgate_actual(4)=(Cgate_actual(5)*(64)+SideLoadCap)*(1)/(fhat);
Cgate_actual(3)=(Cgate_actual(4))*(6/4)/(fhat);
Cgate_actual(2)=(Cgate_actual(3)*4)*(1)/(fhat);
Cgate_actual(1)=(Cgate_actual(2)*5)*(1)/(fhat); 
%!Cgate_actual(1) should be about 11fF if the math is correct
%!
%!Update P if Nhat is different than N
P=P+(Nhat-N);
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% More Fixed Equations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
D=N*F^(1/N)+P;
Delay=N*F^(1/N)*tao+P*pinv; % this delay assumes a single cap as load
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%model load as a phi model
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%!determine the size of the resistors in the inverter that drives the load
%! select the load driving inverter's capacitance
%!CGate*total_gate_length(#####)=Cgate_actual(#####), therefore
%!total_gate_length(#####) = Cgate_actual(#####)/Cgate 
%!(1/4 is NMOS, 3/4 is PMOS)
GateLength_NMOS=round(10^7*Cgate_actual(6)/CGate*1/4)/10^7;
GateLength_PMOS=round(10^7*Cgate_actual(6)/CGate*3/4)/10^7; 
%!note - lengths are rounded to nearest .1 um
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% More Fixed Equations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
R_NMOS=REffNmos/GateLength_NMOS;
R_PMOS=REffPmos/GateLength_PMOS;
R_eff=max(R_NMOS,R_PMOS); % this is the higher of the 2 for the worst case
%
WL_C1=WordLineCap/2;
WL_C2=WordLineCap/2;
%resistance of word line is (length of wl)/(width of wl)*Rsheet
R_WL=(64*40*lambda)/(4*lambda)*rsheetm2;
%
%Use Elmore delay model to get delay of phi model
WL_PhiDelay=R_eff*WL_C1+(R_eff+R_WL)*WL_C2;
%
TOTALDELAY=Delay+WL_PhiDelay;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Print results to the prompt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(' ')
disp(' ')
disp(' ')
disp(' ')
disp(' ')
disp(' ')
disp(' ')
disp(' ')
disp(' ')
disp('--------------------------------------------------------------')
disp ('Results')
disp('--------------------------------------------------------------')
            N %= 6
disp (' ')
            Nhat %= 6
disp (' ')
            F %= 51228
disp (' ')
            fhat %=6.094
disp (' ')
            Cgate_actual %= (see next line)
%[1.13e-14, 1.38e-14, 2.10e-14, 8.56e-14, 5.73e-15, 2.33e-14]

disp (' ')
            D %=46.565
disp (' ')
            Delay %=1.0608e-9 (seconds)
            WL_PhiDelay %= 1.5665348608e-010
            TOTALDELAY %= 1.21750150552821e-009
disp (' ')
            
disp('--------------------------------------------------------------')
disp ('Error')
disp('--------------------------------------------------------------')
Guess_was=Cgate_5 %=6e-15
Result_was=Cgate_actual(5) %=5.73626513158218e-15
Error=Cgate_5-Cgate_actual(5) %=2.63734868417817e-16
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

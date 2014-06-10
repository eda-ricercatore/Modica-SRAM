%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matlab file created by Andrew Mattheisen and Zhiyang Ong        %
% You can contact us at amattheisen@yahoo.com or zhiyang@ieee.org %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This matlab file determines the delay of multiple row decoders, %
% all using 2-4 predecoders and AND trees.                        %
%                                                                 %
% This is a variation of the option1.m file.  In this variation,  %
% The number of stages is increased to the optimal number         %
% (from 6 to 8 stages). One extra buffer is inserted at the DFF   %
% output and the other drives WordLineCap                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp ('THIS IS OPTION 1A')
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
                 %   the NMOSgate
CNDiff=2.244*f/u; % This is for an nmos transistor drain or source
CPDiff=1.679*f/u; % This is for a pmos transistor drain or source
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
%Cgate_6=#####*f;
Cgate_6=91*f;
%!
%!FILL IN FANOUT OF BRANCH
%B=#####*(#####*Cgate_6+SideLoadCap)/Cgate_6;
B=5*4*(64*Cgate_6+SideLoadCap)/Cgate_6;
%
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% More Fixed Equations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PATH EFFORT
F=G*H*B;
%
%Should determine optimum N first, then set stage effort for that N.
Nhat=round(log(F)/log(4));
%Nhat=N;
%
% OPTIMAL STAGE EFFORT
fhat=round(F^(1/Nhat));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% Calculate Actual input Capacitances (Cin) for each stage
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%!FILL IN THIS ENTIRE CIRCUIT
% Cin = (Cout)*(g_i)/(fhat)
Cgate_actual(8)=(WordLineCap)*(1)/(fhat); 
Cgate_actual(7)=(Cgate_actual(8))*(1)/(fhat); 
%----------------------------------------------------------------------------%
Cgate_actual(6)=(Cgate_actual(7))*(6/4)/(fhat); %this is the value of interest
%----------------------------------------------------------------------------%
Cgate_actual(5)=(Cgate_actual(6)*64+SideLoadCap)*(1)/(fhat);
Cgate_actual(4)=(Cgate_actual(5))*(6/4)/(fhat);
Cgate_actual(3)=(Cgate_actual(4)*4)*(1)/(fhat);
Cgate_actual(2)=(Cgate_actual(3))*(1)/(fhat);
Cgate_actual(1)=(Cgate_actual(2)*2)*(1)/(fhat); % should be about 11fF
%                                                 if the math is correct
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
%! model load as a phi model
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%!determine the size of the resistors in the inverter that drives the load
%!select the load driving inverter's capacitance
%!CGate*total_gate_length(#####)=Cgate_actual(#####), therefore
%!total_gate_length(#####) = Cgate_actual(#####)/Cgate 
%!(1/4 is NMOS, 3/4 is PMOS)
GateLength_NMOS=round(10^7*Cgate_actual(8)/CGate*1/4)/10^7;
GateLength_PMOS=round(10^7*Cgate_actual(8)/CGate*3/4)/10^7; 
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
disp('--------------------------------------------------------------')
disp ('Results')
disp('--------------------------------------------------------------')
            N %= 6
disp (' ')
            Nhat %= 8
disp (' ')
            F %= 37479
disp (' ')
            fhat %=4
disp (' ')
            Cgate_actual %= (see next lines)
%[4.311e-15, 8.622e-15, 3.4488e-14, 3.4488e-14
%  9.1968e-14, 3.3288e-15, 8.8768e-15, 3.55072e-14]
disp (' ')
            D %=46.709
disp (' ')
            Delay %= 1.141e-9 (seconds)
            WL_PhiDelay %= 1.042e-10 (seconds)
            TOTALDELAY %= 1.2458e-9 (seconds)
disp (' ')
disp('--------------------------------------------------------------')
disp ('Error')
disp('--------------------------------------------------------------')
Guess_was=Cgate_6 %= 9.1e-14
Result_was=Cgate_actual(5) %= 9.1968e-14
Error=Cgate_6-Cgate_actual(5) %= -9.68000000000002e-16
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

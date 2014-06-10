%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matlab file created by Andrew Mattheisen and Zhiyang Ong        %
% You can contact us at amattheisen@yahoo.com or zhiyang@ieee.org %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This matlab file determines the delay of multiple row decoders, %
% all using 2-4 predecoders and AND trees.                        %
%                                                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp ('THIS IS OPTION 2A_real')
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
NumCols=128;
NumRows=256;
%!
%!FILL IN LOGICAL EFFORT OF PATH
G=125/64; 
%!
%!FILL IN NUMBER OF STAGES
N=8;
%!
%!FILL IN PARACITIC CAPACITANCE
P=11;
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
Cgate_7=34*f;
%!
%!FILL IN FANOUT OF BRANCH
%B=#####*(#####*Cgate_6+SideLoadCap)/Cgate_6;
B=9*8*(16*Cgate_7+SideLoadCap)/Cgate_7;
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
fhat=F^(1/Nhat);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% Calculate Actual input Capacitances (Cin) for each stage
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%!FILL IN THIS ENTIRE CIRCUIT
% Cin = (Cout)*(g_i)/(fhat)
Cgate_actual(8)=(WordLineCap)*(1)/(fhat);
if Cgate_actual(8)<3.833*f
    Cgate_actual(8)=3.833*f;
end
%----------------------------------------------------------------------------%
Cgate_actual(7)=(Cgate_actual(8))*(5/4)/(fhat); %this is the value of interest
if Cgate_actual(7)<3.833*f
    Cgate_actual(7)=3.833*f;
end
%----------------------------------------------------------------------------%
Cgate_actual(6)=(Cgate_actual(7)*16+SideLoadCap)*(1)/(fhat);
if Cgate_actual(6)<3.833*f
    Cgate_actual(6)=3.833*f;
end
Cgate_actual(5)=(Cgate_actual(6))*(5/4)/(fhat);
if Cgate_actual(5)<3.833*f
    Cgate_actual(5)=3.833*f;
end
Cgate_actual(4)=(Cgate_actual(5))*(1)/(fhat);
if Cgate_actual(4)<3.833*f
    Cgate_actual(4)=3.833*f;
end
Cgate_actual(3)=(Cgate_actual(4))*(5/4)/(fhat);
if Cgate_actual(3)<3.833*f
    Cgate_actual(3)=3.833*f;
end
Cgate_actual(2)=(Cgate_actual(3)*8)*(1)/(fhat);
if Cgate_actual(2)<3.833*f
    Cgate_actual(2)=3.833*f;
end
Cgate_check=(Cgate_actual(2)*9)*(1)/(fhat); % should be about 11fF
%                                                 if the math is correct
Cgate_actual(1)=11*f;
%!
%!Update P if Nhat is different than N
P=P+(Nhat-N);
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% More Fixed Equations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%D=N*F^(1/N)+P;
%
g=[1,1,5/4,1,5/4,1,5/4,1];
h=[Cgate_actual(2)*9/Cgate_actual(1), ... %gate 1
   Cgate_actual(3)*8/Cgate_actual(2), ... %gate 2
   Cgate_actual(4)/Cgate_actual(3), ... %gate 3
   Cgate_actual(5)/Cgate_actual(4), ... %gate 4
   Cgate_actual(6)/Cgate_actual(5), ... %gate 5
   (Cgate_actual(7)*16+SideLoadCap)/Cgate_actual(6), ... %gate 6
   Cgate_actual(8)/Cgate_actual(7), ... %gate 7
   WordLineCap/Cgate_actual(8)]; %gate 8
f=[g(1)*h(1), g(2)*h(2), g(3)*h(3), g(4)*h(4), ...
   g(5)*h(5), g(6)*h(6), g(7)*h(7), g(8)*h(8)];
D=f(1)+f(2)+f(3)+f(4)+f(5)+f(6)+f(7)+f(8)+P;
%
% Delay is the delay in seconds
Delay=(f(1)+f(2)+f(3)+f(4)+f(5)+f(6)+f(7)+f(8))*tao+P*pinv;
% this delay assumes a single cap as load
%
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
            N %= 8
disp (' ')
            Nhat %= 8
disp (' ')
            F %= 65162
disp (' ')
            fhat %=3.997
disp (' ')
            Cgate_actual %= (see next lines)
%[%----------------------------------------------------------------------------%
%[1.1e-14, 7.6714e-15, 3.833e-015, 8.475e-15
%  3.387e-14, 1.083e-13, 2.222e-14, 7.106e-14]
disp (' ')
            D %=44.023
disp (' ')
            Delay %= 1.0637e-9 (seconds)
            WL_PhiDelay %= 1.0790e-10 (seconds)
            TOTALDELAY %= 1.1716-9 (seconds)   
disp (' ')
disp('--------------------------------------------------------------')
disp ('Error')
disp('--------------------------------------------------------------')
Guess_was=Cgate_7 %= 3.4e-14
Result_was=Cgate_actual(5) %= 3.38760167673469e-014
Error=Cgate_7-Cgate_actual(5) %= 1.23983232653092e-016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

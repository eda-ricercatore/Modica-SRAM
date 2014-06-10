%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matlab file created by Andrew Mattheisen and Zhiyang Ong        %
% You can contact us at amattheisen@yahoo.com or zhiyang@ieee.org %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This matlab file determines the delay of multiple row decoders, %
% all using 2-4 predecoders and AND trees.                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp ('THIS IS OPTION 1_real')
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
rsheetm1 = .08; % this is in ohms per square from tech file
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
%model load as a pi model
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
%Use Elmore delay model to get delay of pi model
WL_PiDelay=R_eff*WL_C1+(R_eff+R_WL)*WL_C2;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% now we model the SideLoadCap as a pi model, including the input
% capacitance from the gates the SideLoad transmission line drives


%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
GateLength_NMOS_SL=round(10^7*Cgate_actual(4)/CGate*1/4)/10^7;
GateLength_PMOS_SL=round(10^7*Cgate_actual(4)/CGate*3/4)/10^7; 
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% More Fixed Equations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
R_NMOS_SL=REffNmos/GateLength_NMOS_SL;
R_PMOS_SL=REffPmos/GateLength_PMOS_SL;
R_eff_SL=max(R_NMOS_SL,R_PMOS_SL); 
% this is the higher of the 2 for the worst case
%
SL_C1=SideLoadCap/2;
SL_C2=SideLoadCap/2+Cgate_actual(5)*64;
%resistance of word line is (length of wl)/(width of wl)*Rsheet
SL_R=(512*40*lambda)/(4*lambda)*(rsheetm1);
%
%Use Elmore delay model to get delay of pi model
SL_PiDelay=R_eff_SL*SL_C1+(R_eff_SL+SL_R)*SL_C2

% (F)^(1/N)*N*tao+P*pinv
% F=G1*B1*H1 ; G1=3/2, B1=5*4=20, H1=1/2ofSideloadCapacitance
% therfore H1=.5*SideLoadCap+.5*64*Cgate_actual(5)
% P1=5 ; N=4

Delay_Gates1=((3/2)*((.5*SideLoadCap+.5*64*Cgate_actual(5))/(11*10^-15))*20)^(1/4)*4*tao+5*pinv;
temp1=((3/2)*(.5*SideLoadCap+.5*64*Cgate_actual(5))*20)^(1/4)*4*tao
temp2=5*pinv
% G2=3/2, B2=1,H2=1/2*WordLineCap/Cgate_actual(5)
% P1=3 ; N=2

Delay_Gates2=((3/2)*(1/2*WordLineCap/Cgate_actual(5))*1)^(1/2)*2*tao+3*pinv;
TOTALDELAY=Delay_Gates1 + SL_PiDelay + Delay_Gates2 + WL_PiDelay;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Print results to the prompt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(' ')
disp('--------------------------------------------------------------')
disp ('Results')
disp('--------------------------------------------------------------')
            N; %= 6
            Nhat %= 6
            F; %= 51228
            fhat; %=6.094
            Cgate_actual; %= (see next line)
            %[1.13e-14, 1.38e-14, 2.10e-14, 8.56e-14, 5.73e-15, 2.33e-14]
            D; %=46.565
disp('--------------------------------------------------------------')
disp('Total delay using capacitors for C_SL and C_WL')
disp('--------------------------------------------------------------')
            Delay %=1.0608e-9 (seconds)
           
disp('--------------------------------------------------------------')
disp('Total delay using pi models for C_SL and C_WL')
disp('--------------------------------------------------------------')
            TOTALDELAY           
            
disp('--------------------------------------------------------------')
disp('Details of calculation of total delay for pi model')
disp('--------------------------------------------------------------')
            Delay_Gates1
            SL_PiDelay
            Delay_Gates2
            WL_PiDelay 
            TOTALDELAY 
%disp (' ')
%disp('--------------------------------------------------------------')
%disp ('Error')
%disp('--------------------------------------------------------------')
Guess_was=Cgate_5; %=6e-15
Result_was=Cgate_actual(5); %=5.73626513158218e-15
Error=Cgate_5-Cgate_actual(5); %=2.63734868417817e-16








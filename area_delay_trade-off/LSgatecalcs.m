%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matlab file created by Andrew Mattheisen and Zhiyang Ong        %
% You can contact us at amattheisen@yahoo.com or zhiyang@ieee.org %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This matlab file determines the delay of the row decoder.       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
disp('-------------------------------------------------------------------')
disp ('This script calculates gate sizes for WL driver')
disp('-------------------------------------------------------------------')
format long g
format compact
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GENERIC PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m=10^-3;
u=10^-6;
n=10^-9;
p=10^-12;
f=10^-15;
a=10^-18;
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
pinv=.5499e-10; % determined by dettao.m
tao=.0778e-10; % determined by dettao.m
disp('-------------------------------------------------------------------')
disp('Calculate the width of the L-S gate of minimum size for 7 input bits.')
disp('-------------------------------------------------------------------')
disp('')
NMOSW=4;
MINPMOSW=24;
TotalWidth=0.0;
for i=1:7
    TotalWidth=TotalWidth+NMOSW+MINPMOSW*2^(i-1);
end
 TotalWidth
%side load capacitance of the minimum sized L-S decoder
Carea=TotalWidth*lambda^2*14*a/(u*u);  % wire area cap
Cfringe=TotalWidth*lambda*35*a/u; %   wire fringe cap
Cwire=Carea+Cfringe;
%  calculate input gate width of the L-S of minimum size for 7 input bits
n=7; % number of input bits
GateWidth=NMOSW*2^(n-1)+MINPMOSW*2^(n-1);
GateWidth
%input capacitance of the minimum sized L-S decoder
Cout1=GateWidth*lambda*CGate;
disp('----------------------------------------------------------------------')
disp('Note: The vertical wire capacitance is negligible compared to the gate')
disp('capacitance of the decoder')
disp('----------------------------------------------------------------------')
disp('')
disp('-------------------------------------------------------------------')
disp('logical effort for circuit 1')
disp('-------------------------------------------------------------------')
disp('')
Cin1=11*f;
H1=Cout1/Cin1;
G1=1;
B1=2;
F1=G1*B1*H1;
Nhat1=round(log(F1)/log(4))
fhat1=F1^(1/Nhat1)
% Cin = (Cout)*(g_i)/(fhat for each stage)
C1gate3=Cout1*1/fhat1
WN1gate3=round(C1gate3/CGate/4*10^6*10)/10
WP1gate3=round(C1gate3/CGate/4*3*10^6*10)/10
C1gate2=C1gate3*1/fhat1
WN1gate2=round(C1gate2/CGate/4*10^6*10)/10
WP1gate2=round(C1gate2/CGate/4*3*10^6*10)/10
disp('-------------------------------------------------------------------')
disp('logical effort for circuit 2')
disp('-------------------------------------------------------------------')
disp('')
% The input cap to the second circuit is the output capacitance of the LS
% decoder
Cin2=Cout1;
% the output cap of the second circuit is the wordline capacitance
NumCols=256;
Cout2=(4*lambda)*(40*lambda)*NumCols*14*a/(u*u) ...   % wire area cap
+ 2*(40*lambda)*NumCols*35*a/u ... %                  + wire fringe cap
+ 2*CGate*NumCols*4*lambda; %                         + SRAM pass gate cap
H2=Cout2/Cin2;
%G =(      ------ Decoder -------         )*(NAND)  
G2 = (2^6)*(1+3*((1-(1/2^7))/(1-1/2)))/(1+3)*(5/4);
%G2=(5/4);
B2=1;
F2=G2*B2*H2;
Nhat2=round(log(F2)/log(4));
%(have to have odd number of stages)
Nhat2=round(log(F2)/log(4))-mod(Nhat2,2)+1

fhat2=F2^(1/Nhat2)

% Cin = (Cout)*(g_i)/(fhat for each stage)
%alternatively
% Cout = Cin *fhat/g_i
% LS decoder
C2gate1=Cin2
WN2gate1=.4
WP2gate1=2.4
%And Gate
C2gate2=C2gate1*fhat2/((2^6)*(1+3*((1-(1/2^7))/(1-1/2)))/(1+3))
WN2gate2=round(C2gate2/CGate/5*2*10^6*10)/10
WP2gate2=round(C2gate2/CGate/5*3*10^6*10)/10
%Inverters
C2gate3=C2gate2*fhat2/(1)
WN2gate3=round(C2gate3/CGate/4*10^6*10)/10
WP2gate3=round(C2gate3/CGate/4*3*10^6*10)/10
C2gate4=C2gate3*fhat2/(1)
WN2gate4=round(C2gate4/CGate/4*10^6*10)/10
WP2gate4=round(C2gate4/CGate/4*3*10^6*10)/10
C2gate5=C2gate4*fhat2/(1)
WN2gate5=round(C2gate5/CGate/4*10^6*10)/10
WP2gate5=round(C2gate5/CGate/4*3*10^6*10)/10
disp('-------------------------------------------------------------------')
disp('Calculations for Pi model of WL')
disp('-------------------------------------------------------------------')
disp('')

Cout2
Rpiwl=rsheetm2*(40*lambda)*NumCols/(4*lambda);
Rpiwl
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%RESULTS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
-------------------------------------------------------------------
This script calculates gate sizes for WL driver
-------------------------------------------------------------------
-------------------------------------------------------------------
Calculate the width of the L-S gate of minimum size for 7 input bits.
-------------------------------------------------------------------
TotalWidth =
        3076
GateWidth =
        1792
----------------------------------------------------------------------
Note: The vertical wire capacitance is negligible compared to the gate
capacitance of the decoder
----------------------------------------------------------------------
-------------------------------------------------------------------
logical effort for circuit 1
-------------------------------------------------------------------
Nhat1 =
     3
fhat1 =
          4.27386411441149
C1gate3 =
     1.00462529576499e-013
WN1gate3 =
                      10.5
WP1gate3 =
                      31.4
C1gate2 =
     2.35062526292632e-014
WN1gate2 =
                       2.5
WP1gate2 =
                       7.4
-------------------------------------------------------------------
logical effort for circuit 2
-------------------------------------------------------------------
Nhat2 =
     5
fhat2 =
          2.83768591488344
C2gate1 =
             4.293632e-013
WN2gate1 =
                       0.4
WP2gate1 =
                       2.4
C2gate2 =
     1.09518912809823e-014
WN2gate2 =
                       1.8
WP2gate2 =
                       2.7
C2gate3 =
     3.10780276293783e-014
WN2gate3 =
                       3.2
WP2gate3 =
                       9.7
C2gate4 =
     8.81896812662451e-014
WN2gate4 =
                       9.2
WP2gate4 =
                      27.6
C2gate5 =
     2.50254616367284e-013
WN2gate5 =
                      26.1
WP2gate5 =
                      78.3
-------------------------------------------------------------------
Calculations for Pi model of WL
-------------------------------------------------------------------
Cout2 =
             5.681152e-013
Rpiwl =
                     204.8

%}

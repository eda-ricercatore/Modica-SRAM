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
disp('Calculations for optimal gate sizing on entire path')
disp('-------------------------------------------------------------------')
disp('')
Cin3=6*CGate*10^-6;
NumCols=256; 
Cout3=(4*lambda)*(40*lambda)*(NumCols+48)*14*a/(u*u) ...   % wire area cap
+ 2*(40*lambda)*(NumCols+48)*35*a/u ... %                  + wire fringe cap
+ 2*CGate*NumCols*4*lambda; %                         + SRAM pass gate cap
%note the extra 48 in the area and fringe calcs is due to routing distances


H3=Cout3/Cin3;
G3=(2^6)*(1+3*((1-(1/2^7))/(1-1/2)))/(1+3)*(5/4);
B3=2;
F3=H3*G3*B3;
Nhat3=round(log (F3)/log (4))
fhat3=F3^(1/Nhat3)
% gate capacitance ans size calculations
disp(' ')
disp('FF')
disp('-----------------------')
C3gate1=Cin3/(1)
WN3gate1=round(C3gate1/CGate/3*10^6*10)/10
WP3gate1=round(C3gate1/CGate/3*2*10^6*10)/10
disp(' ')
disp('Branched inverter (bi=2)')
disp('-----------------------')
C3gate2=C3gate1*fhat3/(1*2)
WN3gate2=round(C3gate2/CGate/4*10^6*10)/10
WP3gate2=round(C3gate2/CGate/4*3*10^6*10)/10
disp(' ')
disp('Inverter')
disp('-----------------------')
C3gate3=C3gate2*(fhat3)/(1)
WN3gate3=round(C3gate3/CGate/4*10^6*10)/10
WP3gate3=round(C3gate3/CGate/4*3*10^6*10)/10
disp(' ')
disp('Inverter')
disp('-----------------------')
C3gate4=C3gate3*fhat3/(1)
WN3gate4=round(C3gate4/CGate/4*10^6*10)/10
WP3gate4=round(C3gate4/CGate/4*3*10^6*10)/10
disp('7 bit LS Decoder')
disp('-----------------------')
C3gate5=C3gate4*fhat3/(1)
WN3gate5=round(C3gate5/CGate/7*10^6*10/64)/10
WP3gate5=round(C3gate5/CGate/7*6*10^6*10/(2^7))/10
disp(' ')
disp('AND gate')
disp('-----------------------')
C3gate6=C3gate5*fhat3/((2^6)*(1+3*((1-(1/2^7))/(1-1/2)))/(1+3))
WN3gate6=round(C3gate6/CGate/5*2*10^6*10)/10
WP3gate6=round(C3gate6/CGate/5*3*10^6*10)/10
disp(' ')
disp('Inverter')
disp('-----------------------')
C3gate7=C3gate6*(fhat3)/(5/4)
WN3gate7=round(C3gate7/CGate/4*10^6*10)/10
WP3gate7=round(C3gate7/CGate/4*3*10^6*10)/10
disp(' ')
%final load - check that it matches Cout3
disp(' ')
disp('-----------------------')
C3load=C3gate7*fhat3
disp('Should be equal to')
Cout3
disp(' ')
disp('-------------------------------------------------------------------')
disp('Calculations for Delay for path 3')
disp('-------------------------------------------------------------------')
disp('')
%delay is the sum of the gh's +P
Delay3=fhat3*7+1 %+1 is for extra inverter in non-optimal path
DelayinSec=Delay3*tao+pinv*7+pinv %.65psec!!!
disp(' ')
disp('-------------------------------------------------------------------')
disp('Calculations for Pi models')
disp('-------------------------------------------------------------------')
disp('')

CpiWL=Cout3/2
Rpiwl=rsheetm2*(40*lambda)*NumCols/(4*lambda)

CpiBL=(4*lambda)*(40*lambda)*(32*7)*14*a/(u*u)/2 ...   % wire area cap
+ (40*lambda)*(32*7)*35*a/u    %                  + wire fringe cap
RpiBL=rsheetm2*(40*lambda)*(32*7)/(4*lambda)

CpiDATA_OUT=(4*lambda)*(40*lambda)*(32*6)*14*a/(u*u)/2 ...   % wire area cap
+ (40*lambda)*(32*7)*35*a/u    %                  + wire fringe cap
RpiDATA_OUT=rsheetm2*(40*lambda)*(32*6)/(4*lambda)

CpiDATA_IN=(4*lambda)*(40*lambda)*(32*4.5)*14*a/(u*u)/2 ...   % wire area cap
+ (40*lambda)*(32*7)*35*a/u    %                  + wire fringe cap
RpiDATA_IN=rsheetm2*(40*lambda)*(32*4.5)/(4*lambda)

CpiWRITE_EN=(4*lambda)*(40*lambda)*(32*5)*14*a/(u*u)/2 ...   % wire area cap
+ (40*lambda)*(32*7)*35*a/u    %                  + wire fringe cap
RpiWRITE_EN=rsheetm2*(40*lambda)*(32*5)/(4*lambda)

CpiCLK=(4*lambda)*(40*lambda)*(32*9)*14*a/(u*u)/2 ...   % wire area cap
+ (40*lambda)*(32*7)*35*a/u    %                  + wire fringe cap
RpiCLK=rsheetm2*(40*lambda)*(32*9)/(4*lambda)

CpiROW_ADDRESS=(4*lambda)*(40*lambda)*(32*7.5)*14*a/(u*u)/2 ...   % wire area cap
+ (40*lambda)*(32*7)*35*a/u    %                  + wire fringe cap
RpiROW_ADDRESS=rsheetm2*(40*lambda)*(32*7.5)/(4*lambda)

CpiCOL_ADDRESS=(4*lambda)*(40*lambda)*(32*5)*14*a/(u*u)/2 ...   % wire area cap
+ (40*lambda)*(32*7)*35*a/u    %                  + wire fringe cap
RpiCOL_ADDRESS=rsheetm2*(40*lambda)*(32*5)/(4*lambda)


disp(' ')
disp('-------------------------------------------------------------------')
disp('Calculations for precharge control signal circuit')
disp('-------------------------------------------------------------------')
disp('')
Pre_load=3430.4*u*CGate
Pre_in=(.4+1.2)*u*CGate
HPre=Pre_load/Pre_in
FPre=HPre;
NhatPre=round(log(FPre)/log(4))
%D=gh+P
DPre=HPre*tao


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%RESULTS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{

%}

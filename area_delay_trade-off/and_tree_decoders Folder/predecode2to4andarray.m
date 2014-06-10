%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matlab file created by Andrew Mattheisen and Zhiyang Ong        %
% You can contact us at amattheisen@yahoo.com or zhiyang@ieee.org %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This matlab file determines the delay of multiple row decoders, %
% all using 2-4 predecoders and AND trees.                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
format long g


% GENERIC PARAMETERS

m=10^-3;
u=10^-6;
n=10^-9;
p=10^-12;
f=10^-15;
a=10^-18;



% TECHNOLOGY SPECIFIC PARAMETERS

lambda=.1*u;
CGate=2.396*f/u; % This is for an inverter gate
CNDiff=2.244*f/u; % This is for an nmos transistor drain or source
CPDiff=1.679*f/u; % This is for a pmos transistor drain or source



% CONSTANTS

% Gate capacitance of the inverter in the DFFPOSX1 cell connected to the
% output
Cin=4*u*CPDiff+2*u*CNDiff % this is 1.1204e-14



% CIRCUIT SPECIFIC PARAMETERS

%SRAM Macro dimentions
NumCols=128;
NumRows=256;

% LOGICAL EFFORT OF PATH
G=125/64; 

% NUMBER OF STAGES
N=8;

% PARACITIC CAPACITANCE - assume 0 for now - fix this later
P=11;

SideLoadCap=(4*lambda)*(40*lambda)*NumRows*14*a/(u*u) ...    % wire area cap
+ 2*(40*lambda)*NumRows*35*a/u;                          % + wire fringe cap

WordLineCap=(4*lambda)*(40*lambda)*NumCols*14*a/(u*u) ...   % wire area cap
+ 2*(40*lambda)*NumCols*35*a/u ... %                      + wire fringe cap
+ 2*CGate*NumCols*4*lambda; %                          + SRAM pass gate cap

% ELECTRICAL EFFORT OF PATH
H=WordLineCap/Cin; 

% BRANCHING EFFORT GUESS
Cgate_5=6.1*f;
B=2*4*(128*Cgate_5+SideLoadCap)/Cgate_5;

% PATH EFFORT
F=G*H*B;

% OPTIMAL STAGE EFFORT
Fhat=F^(1/N);

Cgate_5_actual=WordLineCap/(Fhat^3)*5/4; % check to see if this matches the guess

% Print results to the prompt
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
%SideLoadCap
%Cin
WordLineCap
%H
%B
%G
%F
%Fhat
Cgate_5
Cgate_5_actual

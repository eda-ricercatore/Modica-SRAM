%{
This is written by Zhiyang Ong (zhiyango@usc.edu;
6004 9194 12) and Andrew Mattheisen (mattheis@usc.edu;
2134 5147 11) for EE 577B, SRAM Project - Part 2, Questions
1 and 2

Performance analysis of the Lyon-Schediwy column decoder
%}
function [delay nrwl ncwl N F P cout w G H cin] = coldecoder()

% Set the accuracy of the mathematical calculations
format long g
format compact
% \lambda design rules are used for this project


global u k p f a lambda gamma c_diff_nmos c_diff_pmos
global R_eff_pmos R_eff_nmos R_eff_inv
global l_byte l_halfword l_word numbits
global c_gate_inv rc_inv nmos pmos
global R_sheet tau pinv R_nmos R_pmos cfringe carea



% Number of stages in the Lyon-Schediwy decoder
N=1;


%{
Input capacitance of the decoder = output capacitance of the
positive edge D flip-flop

Each D flip-flop is connected to an input in the row or column
decoder

The output of the D flip-flop from the OSU library is available
at the output of an inverter, which has a PMOS width of 40 
lambda and a NMOS width of 20 lambda

For a gate with widths of its NMOS and PMOS transistors 
given in terms of lambda, its output capacitance is determined
as the sum of the following: (width of NMOS in micron) *
(diffusion capacitance per micron of width for contacted NMOS
transistor) + (width of PMOS in micron) * (diffusion
capacitance per micron of width for contacted PMOS transistor)
%}
cout_dff=4*u*c_diff_pmos + 2*u*c_diff_nmos;
cin=4*cout_dff;






% Number of words considered for the width of the SRAM macro
for k=1:7
    % Start processing the number of words per row from 1
    numwds(k)=2^(k-1);
end
% Include the case for 1 row - put all words on a row
numwds=[numwds numbits/l_word];


%{
Array for the number of cells per word line in the SRAM submacro
Sizes (number of cells per word line) should be word addressable,
and cannot be an odd multiple (>1) of bytes or halfwords.
%}
%ncwl=[l_byte l_halfword];
%ncwl=[l_word l_word];
ncwl=[];
%t=[0];
for i=1:length(numwds)
	% Populate the array with other sizes
    ncwl=[ncwl (numwds(i)*l_word)];
   
%   t=[t (i*l_word)];
end

%{
For each size (number of SRAM cells) per word line in the SRAM
submacro, determine the number of required rows for that SRAM
submacro
%}
nrwl=numbits./ncwl;
% Round up the number of rows required to the nearest integer
nrwl=ceil(nrwl);






% Evaluating the delay/speed of the row decoder
for j=1:length(nrwl)
%{
    Determine the capacitance of the word line, cwl,
    and the bit line, cbl
    cwl & cbl includes area and fringing capacitance
%}
    % Fringing capacitance of cwl, cfringe
    cfringe_wl(j)=2*(40*lambda*nrwl(j))*(35*a/u);
    % Area capacitance of cwl, carea
    carea_wl(j)=(4*lambda)*(40*lambda*nrwl(j))*(14*a/(u^2));
    %cwire = 256*40*lambda*((4*lambda)*14*a/(u*u)+2*35*a/u)
    cwire_wl(j)=carea_wl(j)+cfringe_wl(j);
    
    
    
%    disp('print c_gate_inv')
%    c_gate_inv
%{
C_{word line} = Capacitance of the word line
= cwire_wl + 2*cgate*(number of cells in per row)

### IMPORTANT ASSUMPTION!!!
Assume that the gate capacitance of the inverter is the
same as that of the NMOS, since the gate capacitance is
dependent on the widths of the transistors
%}
    csram(j)=2*c_gate_inv*nrwl(j)*4*lambda;
    cwl(j)=cwire_wl(j)+csram(j);
    
    
    
    
    
    % Fringing capacitance of cbl, cfringe
    cfringe_bl(j)=2*(40*lambda*ncwl(j))*(35*a/u);
    % Area capacitance of cbl, carea
    carea_bl(j)=(4*lambda)*(40*lambda*ncwl(j))*(14*a/(u^2));
    %cwire = 256*40*lambda*((4*lambda)*14*a/(u*u)+2*35*a/u)
    cdiff_bl(j)=4*lambda*c_diff_nmos*(ncwl(j));
    cwire_bl(j)=carea_bl(j)+cfringe_bl(j);
    c_bl(j)=cwire_bl(j)+cdiff_bl(j);

%{
C_{word line} = Capacitance of the word line
= cwire_wl + 2*cgate*(number of cells in per row)

### IMPORTANT ASSUMPTION!!!
Assume that the gate capacitance of the inverter is the
same as that of the NMOS, since the gate capacitance is
dependent on the widths of the transistors
%}
%%%%%%    cbl=cwire_bl+2*cgate*nrbl(j);  
    
    
%{
Number of rows is given by nrwl(j) = output of decoding block
Number of input bits to the decoder = log2(nrwl(j)) = n
%}
    n(j)=log2(ncwl(j));
    
    
    % Size of the smallest PMOS in the Lyon-Schediwy decoder
    % w=gamma*((1 - 1/nrwl(j)) / (1 - 0.5))
    w(j)=2*gamma*(1 - 1/ncwl(j));
%{
Determine the total parasitic capacitance of the n-input
NOR gate; units are given in p_{inv}
\cite{Sutherland99}, pp. 81
p (\Sum w_d) / (1 + gamma)
%}
    P(j) = n(j)
    
    % Logical effort per input for the Lyon-Schediwy decoder
    G(j)=(2^(n(j)-1)) * ((1 + w(j)) / (1 + gamma));
    % Logical effort for gate is unity, while that per input
    % is given as above
    %%%G(j)=1
%{
Electrical effort of the decoder
= C_{out} of path / C_{in} of path
%}
    cout(j) = c_bl(j); % + cwire_bl(j);
    H(j)= cout(j) / cin;
    %H(j)= cwl(j) / cout_dff;
    
    
    
    % Determine the path effort of the Lyon-Schediwy decoder
    % F=G*H
    F(j) = G(j) * H(j);
    
    
    N(j)=ceil(log(F(j))/log(4))
    % Delay of the path
    D(j) = N(j)*(F(j)).^(1/N(j)) + P(j);
    
%{
When modeling interconnects as parasitic capacitance, the
branching factor can be determined as follows
\cite{Sutherland99}, pp. 175:
braching effort at a wire driving a gate load
= (C_{gate} + C_{wire}) /  C_{gate}
%}
end


%disp('disp c_gate_inv')
%c_gate_inv

delay=D;




%%{
nrwl
ncwl
%}


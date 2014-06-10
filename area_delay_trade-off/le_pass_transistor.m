%{
This is written by Zhiyang Ong (zhiyango@usc.edu;
6004 9194 12) and Andrew Mattheisen (mattheis@usc.edu;
2134 5147 11) for EE 577B, SRAM Project - Part 2, Questions
1 and 2

Determine the logical effort for a row decoder
@param width is the width of the pass transistor, given in meters
@param transistor_type is the type of transistor
@return le_fall is the logical effort of the pass transistor
 when its output is falling. 
@return le_rise is the logical effort of the pass transistor
 when its output is rising. 


Note that the logical effort of the pass transistor should be
calculated with a driver (driving gate to drive it, or provide
it with input signals)

Logical effort cannot be determined for the pass transistor
without a driver since it would have not input logic
%}
function [le_fall le_rise] = le_pass_transistor(width,transistor_type)
% Set the accuracy of the mathematical calculations
format long g
format compact


global u k p f a lambda gamma c_diff_nmos c_diff_pmos
global R_eff_pmos R_eff_nmos R_eff_inv
global l_byte l_halfword l_word numbits
global c_gate_inv rc_inv nmos pmos
global R_sheet tau pinv R_nmos R_pmos cfringe carea



% Effective resistance for contacted NMOS transistor
r_eff_nmos=2.5832*k*u/width;
% Effective resistance for contacted PMOS transistor
r_eff_pmos=7.8648*k*u/width;


% Printing the values of a variable in an array will not work
%disp(['width:' width '==' transistor_type '-m'])
disp('Value of width')
width
disp('Value of transistor_type')
transistor_type



sz_w=size(width)
sz_l=size(lambda)
ww=width./(4*lambda);
if strcmp(transistor_type,'nmos')
    %rc=r_eff_nmos*c_diff_nmos*width
    width;
    le_fall=(ww*R_eff_nmos*(1+1/ww)*c_diff_nmos*width)/(rc_inv);
    le_rise=(ww*R_eff_nmos*(1+2/ww)*c_diff_nmos*width)/(rc_inv);
    
    
    %le_ff=(34/4)*(1+1/34)
    %le_rr=(34/4)*(1+2/34)
    
%    le2=(0.5*R_eff_nmos*1.5*c_gate_inv*width) / (rc_inv)
else
    le_fall=(ww*R_eff_pmos*(1+2/ww)*c_diff_pmos*width) / (rc_inv);
    le_rise=(ww*R_eff_pmos*(1+1/ww)*c_diff_pmos*width) / (rc_inv);
end


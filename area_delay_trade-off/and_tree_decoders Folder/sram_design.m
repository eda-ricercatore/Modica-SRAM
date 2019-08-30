%{
This is written by Zhiyang Ong and Andrew Mattheisen for
	EE 577B, SRAM Project - Part 2, Questions 1 and 2

Performance analysis of the SRAM design


Important Notes:
#   Calculating the logical effort for a pass transistor
    requires the use of a driving gate (driver) to send/drive
    the pass transistor
#   Calculating the branching effort of a leg (which is part
    of n branches) is correct, in the way we modeled the
    different branches and the capacitance of the wire
#   The modeling and the calculation of the side load
    capacitance in the Lyon-Schediwy decoder needs to be
    checked


When using the principles/method of logical effort to
determine the sizes of the gates,
%}
% Set the accuracy of the mathematical calculations
format long g
format compact

global u k f a lambda gamma c_diff_nmos c_diff_pmos
global R_eff_pmos R_eff_nmos R_eff_inv
global l_byte l_halfword l_word numbits
global c_gate_inv rc_inv nmos pmos

% Scale for a micro, \mu 
u=1e-6;
% Scale for a kilo, k
k=1e3;
% Scale for a femto, f
f=1e-15;
% Scale for a atto, a
a=1e-18;
% Length of a \lambda
lambda = 0.1*u;



%{
Ratio of PMOS width to NMOS width in the template/unit
inverter
%}
gamma=3;





%{
Diffusion capcitance per micron of width for contacted NMOS
transistor
%}
c_diff_nmos=2.244*f/u;
%{
Diffusion capcitance per micron of width for contacted PMOS
transistor
%}
c_diff_pmos=1.679*f/u;

% Effective resistance for a PMOS transistor
R_eff_pmos=3.277*k;
% Effective resistance for a NMOS transistor
R_eff_nmos=3.229*k;

%{
If the resistances, or logical effort for a rise or fall in
the output are different, evaluate their delays separately
%}
% Effective resistance for a template (unit) inverter
R_eff_inv=0.5*(R_eff_nmos + R_eff_pmos);


% Length of a byte, in terms of bits
l_byte=8;
% Length of a halfword, in terms of bits
l_halfword=2*l_byte;
% Length of a word, in terms of bits
l_word=2*l_halfword;
% Number of bits in the SRAM macro
numbits=2^15;


%{
Effective gate capcitance per micron of width for a
template (unit) inverter
%}
c_gate_inv=2.396*f/u;
%{
Effective RC delay/product for the template (unit)
inverter
%}
rc_inv=R_eff_inv*c_gate_inv*4*lambda;




%{
Size of the pass transistors in the write circuitry that
have to be varied - Use units of meters

The pass transistors are M7, M8, M9, and M10 in the write
circuit given in the Part 1 description of the project
%}
width_pass_write=34*lambda;
% Keywords to indicated the type of transistor represented
nmos='nmos';
pmos='pmos';
% =======================================================

%{
Determine the logical effort for the rising and falling
output for a MOS pass transistor
%}
% Logical effort for a NMOS pass transistor
[le_r_pass_n le_f_pass_n]=le_pass_transistor(width_pass_write,'nmos')
% Logical effort for a PMOS pass transistor
[le_r_pass_p le_f_pass_p]=le_pass_transistor(width_pass_write,'pmos')

%{
If the logical effort for the rise and the fall in a gate
(e.g. transmission gate) or a pass transistor, do the delay
calcaulations for the rise and the fall in the output
separately

%}

% =======================================================

disp('Commence architecture exploration...')

%{
Get the delay for the row decoder based on the
Lyon-Schediwy decoder design
%}
[d_r_decoder nrows ncols] = rowdecoder()

%{
Get the delay for the row decoder based on the
Lyon-Schediwy decoder design
%}
[d_c_decoder nrows1 ncols1] = coldecoder()

[d_decoder]=max(d_r_decoder,d_c_decoder)
%{
Get the delay through the SRAM Macros
The bit and bit_bar lines drive a NMOS transistor each
Also, the charging up of data values in the SRAM cells
cannot be modeled as logical effort cannot be used to
determine the delay through logic that is dependent on
charges dropping from a certain/given level to another
%}

%{
The logic driven in the SRAM macro is a NMOS in the
selected SRAM cells
%}
le_sram_cell_f=le_f_pass_n;
le_sram_cell_r=le_r_pass_n;
% The logic driven in the column mux is a NMOS transistor
le_col_mux_f=le_f_pass_p;
le_col_mux_r=le_r_pass_p;
% The logic driven in the sense amp is a PMOS transistor
le_sense_amp_f=le_f_pass_p;
le_sense_amp_r=le_r_pass_p;

%{
Remember to take care of the dummy row and dummy column
Also, check that the diffusion capacitances at the output
of the NAND gate is considered as part of the bit and
bit_bar lines' load; the NAND gates have inputs from the
output of the row decoder, and the control signal to 
precharge the bit or bit_bar line
%}



% Model the capacitances...















disp('End architecture exploration...')










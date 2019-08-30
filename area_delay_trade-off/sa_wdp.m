%{
This is written by Zhiyang Ong and Andrew Mattheisen for
	EE 577B, SRAM Project - Part 2, Questions 1 and 2

Performance analysis of the Lyon-Schediwy row decoder
%}



% Sizing the transmission gate for the column MUX/DEMUX
format long g
clear

global lambda

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define constants for the script
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Ratio of PMOS to NMOS widths:
gamma=3;
% Magnitude for 1 micro
u=1e-6;
% Magnitude for 1 femto
f=1e-15;
% Magnitude for 1 atto
a=1e-18;
% Value for \lambda
lambda=0.1*u;
% Maximum width of NMOS in the transmission gate
max_wn=50;
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

% Wire capacitance for the bit and bit_bar lines
cwire=128*40*lambda*(2*35*f/u+(u*4*lambda)*14*a/(u*u));
% Capacitance for a word line
cbl=128*(4*lambda)*2.244*f/u+cwire;

% Resistance of the inverter in units of R
r_inv=1;
% Capacitance of the inverter in units of R
c_inv=1;
% RC product for unit inverter
rc_inv=r_inv*c_inv;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Size the transistors in a transmission gate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
index=1;
for nw=1:max_wn
    for pw=1:max_wn*gamma
        % Resistance of the transmission gate
        r_tg_r=1/(1/(2*nw) + 1/pw);
        % Capacitance of the transmission gate
        c_tg_r=pw+nw;
        % RC product of the transmission gate
        rc_tg_r(index)=r_tg_r*c_tg_r;
        
        
        % Resistance of the transmission gate
        r_tg_f=1/(1/(2*pw) + 1/nw);
        % Capacitance of the transmission gate
        c_tg_f=pw+nw;
        % RC product of the transmission gate
        rc_tg_f(index)=r_tg_f*c_tg_f;
       
        
        % Determine the max of rc_tg_f(index) & rc_tg_r(index)
        rc_tg(index)=max(rc_tg_f(index),rc_tg_r(index));
        pw_cur(index)=pw;
        nw_cur(index)=nw;
        le_tg(index)=rc_tg(index)/rc_inv;
        index=index+1;
    end
end



%{
Determine the values for the transmission gate with the
minimum logical effort
%}
for k=1:length(rc_tg)
    if rc_tg(k) == min(rc_tg)
        rctg=rc_tg(k)
        nwcur=nw_cur(k)
        pwcur=pw_cur(k)
        letg=le_tg(k)
        %break
    end
    
end



if (nwcur==1) && (pwcur==1)
    disp('Transmission gate has minimum sized transistors')
end
% Logical effort of the transmissiono gate is 4/3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determine the path delay for the path connect the SRAM
% block to the output flip-flop, and perform buffer
% insertion, if necessary

%{
Width of the PMOS pass transistors in the sense amplifier
that is connected to the bit and bit_bar lines

Sense amp does NOT use transmission gates
%}
width_sa=50*lambda
size_width_sa=size(width_sa)
% Logical effort for PMOS pass transistor
%[le_r_n le_f_n]=le_pass_transistor(width_sa,'pmos')
[le_r_n le_f_n]=le_pass_transistor(5e-6,'pmos')
% Logical effort for the pass transistor =max(le_r_n, le_f_n)
le_pt=max(le_r_n, le_f_n);


%{
Perform buffer insertion in the critical paths to minimize delay
logical effort for: col_mux sense_amp ff

Note that logical effort for col_mux, transmission gate, and 
flip-flop are the same. Since the output stage of the 
flip-flop in the OSU library is a transmission gate, I use 
the diffusion capacitances of the output stage as the output
capacitance of the flip-flop.
This "equates" the logical effort of the flip-flop with that
of the transmission gate.

Also, note that the impedance of the input/output (I/O) pads
are not considered. Else, I would model them as a load for
the output flip-flop.
%}
g=[ letg le_pt letg ];
% Logical effort of path
G=1;
for i=1:length(g)
    G=G*g(i);
end
% Determine the load capacitance at the output of the path
cload=(u*c_diff_nmos+2*u*c_diff_pmos)
%{
Use load capacitance to help determine the electrical
effort H
%}
G
H=(cload+cbl)/(8*cbl+24*lambda*c_diff_pmos)

%{
Determine the path effort for the aforementioned path
B=1, since there is no fan-out
%}
F=G*H;
N_hat=log(F)/log(4)
if N_hat < 1
   disp('Path effort is less than one') 
elseif (N_hat <= 3) && (N_hat >= 1)
   disp('Path is optimized')
else
    disp('Perform buffer insertion in the path')
end















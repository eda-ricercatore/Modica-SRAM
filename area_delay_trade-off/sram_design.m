%{
This is written by Zhiyang Ong (zhiyango@usc.edu;
6004 9194 12) and Andrew Mattheisen (mattheis@usc.edu;
2134 5147 11) for EE 577B, SRAM Project - Part 2, Questions
1 and 2

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
#   Since some of the variable names are resued between
    functions, or even reused for other purposes in this
    function/script, the workspace for MATLAB should be
    cleared so that the variables that do not contain
    correct values will be destroyed. New variables that
    are subsequently created will need to be initialized.
    Thus, the script will not use incorrect values when
    the script is rerun


Note that the input flip-flops are connected to the decoders,
which is subsequently connected to the NAND gate for the
control signals for the SRAM's word line. Next, the NAND
gate is connected to the SRAM...
Logic can be placed between them: consider buffer insertion
after the NAND gate, and a possible buffer insertion before
the row decoder



Use the principles/method of logical effort to determine
the sizes of the gates
%}
clear
% Set the accuracy of the mathematical calculations
format long g
format compact

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXPORT GLOBAL VARIABLES THAT WILL BE DEFINED SUBSEQUENTLY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global u k p f a lambda gamma c_diff_nmos c_diff_pmos
global R_eff_pmos R_eff_nmos R_eff_inv
global l_byte l_halfword l_word numbits
global c_gate_inv rc_inv nmos pmos
global R_sheet tau pinv R_nmos R_pmos cfringe carea


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DECLARE AND INITIALIZE GLOBAL VARIABLES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Scale for a micro, \mu 
u=1e-6;
% Scale for a kilo, k
k=1e3;
% Scale for a femto, f
f=1e-15;
% Scale for a atto, a
a=1e-18;
% Scale for a pico, p
p=1e-12;
% Length of a \lambda
lambda = 0.1*u;
%{
Sheet resistance of M2 from the TSMC technology file (\ohm/square)
A square is (\lambda)^2w
%}
R_sheet = 0.08;
%{
\tau is the unit for path or stage effort when estimating the
delay in a digital circuit using the principles/method of
logical effort; unit for the constant \tau is ps (picoseconds)
%}
tau = 14.34 * p;
%{
p_{inv} is the unit for the path's or stage's parasitic delay
when estimating the delay in a digital circuit using the 
principles/method of logical effort; unit for the constant 
p_{inv} is ps (picoseconds)
%}
pinv = 53.65 * p;





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
% Effective resistance for a PMOS transistor in \ohm*m
R_pmos=7.8648*k*u;
% Effective resistance for a NMOS transistor in \ohm*m
R_nmos=2.5832*k*u;


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

% Fringe capacitance per micron
cfringe = (35*a/u);
% Area capacitance per micron^2
carea = 14*a/(u^2);


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









%{
Note to authors and readers:
Use "nrows" to enumerate the number of designs that have
been considered by Zhiyang Ong.

Use src_i to indicate the selected (i^{th}) design that
have been developed by Zhiyang Ong.

Also, the designs from Andrew Mattheisen were each
separately analyzed by its own script. The results from
these scripts were collated by Zhiyang Ong and Andrew 
Mattheisen and placed/entered into arrays for delay and
area. These arrays were subsequently concatenated into the
arrays for the delays and areas of Zhiyang's designs.
Finally, a graph of delay versus area is obtained to help
determine optimal trade-off(s) for area and delay; there
may be multiple equivalent trade-offs, especially during
multi-objective optimization - see Pareto optimums.
%}

% Width of M7, M8, M9, and M10 is 34 \lambda
R_wr_nmos=R_nmos/(34*lambda);

% =======================================================


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DETERMINE THE VALUES FOR LOGICAL AND ELECTRICAL EFFORT,
% AND PARASITICS FOR EACH STAGE IN THE CRITICAL PATH
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%{
Logical effort for the 2-input NAND gate, assuming equal 
rise and fall delays
%}
le_nand=(gamma+2)/(gamma+1);



% =======================================================

disp('Commence architecture exploration...')

%{
Get the delay for the row decoder based on the
Lyon-Schediwy decoder design
%}
[d_r_decoder nrows ncols nrd frd prd cor wr gr hr cri] = rowdecoder()

%{
Get the delay for the row decoder based on the
Lyon-Schediwy decoder design
%}
[d_c_decoder nrows1 ncols1 ncd fcd pcd coc wc gc hc cci] = coldecoder()


%{
Miscellaneous: 
To determine the maximum element and the first index in
which the element occurs in the array a, try:
[b c]=max(a)
where b is the maximum element and c is its index
%}


[d_decoder]=max(d_r_decoder,d_c_decoder)















%{
Select the fastest SRAM design based on the arrangement/
architecture of the row and column decoders, and the
SRAM block of macros

This will consider for each architecture whether the
row decoder or column decoder has a larger delay, and
only uses the larger delay for determining the delay
on the critical path

Thus, for each design, determine the larger delay, and
its index for calculating the worst-case path delay in
the SRAM

Also, get the index src_i of the selected architecture
BASED SOLELY
%}
[sel_rc src_i] = min(d_decoder);

% Determine the N, F, and P of the selected decoder
if d_r_decoder(src_i) >= d_c_decoder(src_i)
%{
    The row decoder has a larger delay; use its values
    for N, F, and P; also do so for co
%}
    Ndec=nrd(src_i);
    Fdec=frd(src_i);
    Pdec=prd(src_i);
    co=cor(src_i);
    dec_w=wr(src_i);
else
%{
    The column decoder has a larger delay; use its
    values for N, F, and P
%}
    Ndec=ncd(src_i);
    Fdec=fcd(src_i);
    Pdec=pcd(src_i);
    co=coc(src_i);
    dec_w=wc(src_i);
end

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


IMPORTANT! FIX THIS IN THE DETAILED CALCULATION


ALSO, MAKE SURE THAT I DO POWER ANALYSIS,
ALONG WITH AREA CALCULATION AND VARIATION-AWARE
ANALYSIS BY VARYING THE VDD AND TEMPERATURE TILL
NOTICEABLE DEGRADATION IS SEEN IN THE PERFORMANCE 
OF THE SRAM

ALSO, IN TODAY'S CLASS (10 OCT), WE HAVE TO DO THE
MITIGATION EFFECTS FOR RADIATION FROM MICROARCHITECTURAL
LEVEL TILL LAYOUT LEVEL, AND ALSO IMPROVE PROCESS 
TECHNOLOGY

DFM TAKES CARE OF PROCESS VARIATION, BUT NOT RADIATION 
EFFECTS

MAKE IT CLOCK ALIGNED, BUT DO NOT WAIT FOR THE CLOCK
NOT BOUND TO THE CLOCK

SEE HIS NOTES DRAWN FOR JAY FOZDAR (PUT IN ACKNOWLEDGEMENTS)


ALSO, DO POWER ANALYSIS WITH NANOSIM
%}



% Model the capacitances and resistances as a pi-model













%{
Determine the value of the delay D in ps (picoseconds)

Multiply the value of the path/stage, f/F, by \tau and
the value of the parasitic delay of the path/stage,
p_{inv}

Since the parasitic of the stage/path and the path/stage
effort is measured in ps, the delay for the digital
circuit/system can be estimated in ps (picoseconds).


This implies that I will need to obtain the stage/path
effort separately from the parasitic delay for that
stage/path. Also, since the delay calculation depends on
the number of stages in the path, the values of N for
different configurations of the SRAM block is needed for
calculating delay. Thus, the functions to determine the
delay of the row and column decoders will also have to 
return the values of stage/path effort and parasitic delay.

Note that from the architecure exploration, the dimensions
for the SRAM block of macros is 128 rows * 256 columns;
d_decoder = [34.0642311319592 33.0613501080333
    32.0555818431364 31.0440203544749 33.3320807611097
    38.303428218463 43.2801252071473 63.2196124069517]
%}

%{
Ignore the column mux, since the row decoder has a larger
delay than the column decoder and the column mux.

I should implement a if-else construct to handle this...
%}

% For read operations... deocder->SRAM cell->sense amp
%{
Output capacitance of the sense amplifier
= Sum the diffusion capacitances of M5 and M1 (PMOS), 
and M3 (NMOS)
%}
csa=c_diff_pmos*((50+4)*lambda) + c_diff_nmos*(50*lambda);
% For rising output: F*g*h*b (b for sram block)
Fsram_rr=Fdec*le_sram_cell_r*le_sense_amp_r*(csa/co)*ncols(src_i)
Nsram_rr=Ndec+1+1;
Psram_rr=Pdec+1+1;
% D = N*[F^(1/N)]*\tau + P*p_{inv}
Dsram_rr=Nsram_rr*(Fsram_rr^(1/Nsram_rr))*tau + Psram_rr*pinv
Dsram_rr_tau=Nsram_rr*(Fsram_rr^(1/Nsram_rr)) + Psram_rr
% For falling output: F*g*h*b (b for sram block)
Fsram_rf=Fdec*le_sram_cell_f*le_sense_amp_f*(csa/co)*ncols(src_i);
Nsram_rf=Ndec+1+1;
Psram_rf=Pdec+1+1;
Dsram_rf=Nsram_rf*(Fsram_rf^(1/Nsram_rf))*tau + Psram_rf*pinv
Dsram_rf_tau=Nsram_rf*(Fsram_rf^(1/Nsram_rf)) + Psram_rf









% For write operations... deocder->SRAM cell
%{
Output capacitance of the write data path
= difussion capacitance for PMOS (M11 or M12, and M5 or M6)
+ difussion capacitance for NMOS (M9 or M10)
%}
cwc=c_diff_pmos*((120+50)*lambda) + c_diff_nmos*(34*lambda);
% For rising output: F*g*h*b (b for sram block)
Fsram_wr=Fdec*le_sram_cell_r*(cwc/co)*ncols(src_i)
Nsram_wr=Ndec+1;
Psram_wr=Pdec+1;
% D = N*[F^(1/N)]*\tau + P*p_{inv}
Dsram_wr=Nsram_wr*(Fsram_wr^(1/Nsram_wr))*tau + Psram_wr*pinv
Dsram_wr_tau=Nsram_wr*(Fsram_wr^(1/Nsram_wr)) + Psram_wr
% For falling output: F*g*h*b (b for sram block)
Fsram_wf=Fdec*le_sram_cell_f*(cwc/co)*ncols(src_i);
Nsram_wf=Ndec+1;
Psram_wf=Pdec+1;
Dsram_wf=Nsram_wf*(Fsram_wf^(1/Nsram_wf))*tau + Psram_wf*pinv
Dsram_wf_tau=Nsram_wf*(Fsram_wf^(1/Nsram_wf)) + Psram_wf







% Expected delay estimation for the SRAM design
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%{
IMPORTANT!!!

INSERT STUFF HERE TO CALCULATE THE AREA AND DELAY FOR THE
DIFFERENT SRAM DESIGNS


CONSIDER THE SRAM MACROS AS A HUGE BLOCK (FOR SIMPLICITY)
ESTIMATE THE DELAYS THROUGH THAT
[SIMPLICITY WAS CHOSEN SINCE THIS SHOULD BE DONE "QUICKLY"]
%}










% Dimension for w, in units of (4*lambda)
%%%%%%%%%%%%% dec_w


%{
For the selected SRAM design, the Lyon-Schediwy decoder is used
to design the row and column decoders

Also, the SRAM macros are broken into 4 quadrants of
32*64 \lambda^2
%}



%powers_all=[2.^(log2(nrows)-1)]
%{
powers_all=[2.^(log2(nrows))]
for i=1:length(powers_all)
    sum_powers(i)=sum(powers_all(1:i))
end
%}
t_log2_nrows=log2(nrows);
sum_powers=zeros([1 length(t_log2_nrows)])
sum_pow_c=zeros([1 length(t_log2_nrows)])
for i=1:length(t_log2_nrows)
    for j=1:t_log2_nrows(i)
        sum_powers(i)=sum_powers(i)+2^j;
     end
    for k=1:(max(t_log2_nrows)-t_log2_nrows(i))
        sum_pow_c(i)=sum_pow_c(i)+2^k
    end
end
sum_powers = 1+sum_powers;
sum_pow_c = 1+sum_pow_c;
%[1 3 7 15 31 63 127 255 511 1023]



% Determine powers of 2 up to 7
% powers=[1     2     4     8    16    32    64   128]
powers_r=[2.^([1:7]-1)]
%{
Length of the word line in the row decoder
= Length (horizontal) dimension of the row decoder
Sum the elements in powers to determine the length of
the wordline in the pull-up network (sum(powers_r)),
and add it to the length of the wordline in the pull-down 
network
%}
l_wl_row_dec=(sum_powers+2.*log2(nrows)).*wr*4*lambda;
l_wl_row_dec_lambda=(sum_powers+2*log2(nrows)).*wr*4
% Width of the word line in the row decoder
w_wl_row_dec=4*lambda;
w_wl_row_dec_lambda=4


%{
Height (vertical dimension) of the row decoder
= (width of wire pitch + routing pitch) * 2^n

Multiply result by 2 since there exists 2 vertical
tracks for 2 signals
%}
%l_ht_row_dec=2*lambda*(2^Ndec);
l_ht_row_dec=(4+4)*lambda*(2.^t_log2_nrows);
% Width of the signal lines running along the height
w_ht_row_dec=4*lambda;





%{
Resistance of the interconnects in the row decoder that
are connected to the inputs - signal line
%}
% Resistance R_rdec_ht of the signal line in the row decoder
% R_rdec_ht = twice the vertical height of the house
R_rdec_ht=R_sheet*(2*l_ht_row_dec.*w_ht_row_dec);
% R_rdec_ht = R_rdec_ht+length of row decoder (for worst case path);
R_rdec_ht=R_rdec_ht+R_sheet*(l_wl_row_dec.*w_wl_row_dec);
% Fringing capacitance of signal line in the row decoder
cfringe_rdec_ht=2*(2*l_ht_row_dec+l_wl_row_dec)*cfringe;
% Area capacitance of signal line in the row decoder
carea_rdec_ht=((2*l_ht_row_dec+l_wl_row_dec).*w_ht_row_dec)*carea;
% Diffusion capacitance of signal line in the row decoder
cdiff_rdec_ht=4*lambda*(c_diff_nmos*Ndec + sum(powers_r)*c_diff_pmos);
%{
Capacitance of the signal line wire in the row decoder
= area capacitance + fringing capacitance
%}
cwire_rdec_ht=carea_rdec_ht+cfringe_rdec_ht;
%{
Capacitance of the signal line in the row decoder, including
the diffusion capacitance
%}
c_rdec_ht=cwire_rdec_ht + cdiff_rdec_ht;













% Resistance R_rdec_wl of the word line in the row decoder
R_rdec_wl=R_sheet*(l_wl_row_dec.*w_wl_row_dec);
% Fringing capacitance of word line in the row decoder
cfringe_rdec_wl=(2*l_wl_row_dec)*cfringe;
% Area capacitance of word line in the row decoder
carea_rdec_wl=(l_wl_row_dec.*w_wl_row_dec)*carea;
% Diffusion capacitance of word line in the row decoder
cdiff_rdec_wl=4*lambda*(c_diff_nmos.*t_log2_nrows + sum(powers_r)*c_diff_pmos);
%{
Capacitance of the word line wire in the row decoder
= area capacitance + fringing capacitance
%}
cwire_rdec_wl=carea_rdec_wl+cfringe_rdec_wl;
%{
Capacitance of the word line in the row decoder, including
the diffusion capacitance
%}
c_rdec_wl=cwire_rdec_wl + cdiff_rdec_wl;












% Determine powers of 2 up to (3-1)
% powers=[1     2     4]
powers_c= max(t_log2_nrows)-t_log2_nrows%[2.^([1:3]-1)]

%{
Length of the word line in the column decoder
= Length (horizontal) dimension of the column decoder
Sum the elements in powers to determine the length of
the wordline in the pull-up network (sum(powers_r)),
and add it to the length of the wordline in the pull-down 
network
%}
l_wl_col_dec=(sum_pow_c+2.*powers_c).*wc*4*lambda;
l_wl_col_dec_lambda=(sum_pow_c+2.*powers_c).*wc*4
w_wl_col_dec=4*lambda;

%{
Height (vertical line) of the column decoder
= (width of wire pitch + routing pitch) * 2^n

Multiply result by 2 since there exists 2 vertical
tracks for 2 signals
%}
l_ht_col_dec=(4+4)*lambda.*[2.^powers_c]
w_ht_col_dec=4*lambda;





% Resistance R_rdec_wl of the word line in the column decoder
R_cdec_wl=R_sheet*(l_wl_col_dec.*w_wl_col_dec);
% Fringing capacitance of word line in the column decoder
cfringe_cdec_wl=(l_wl_col_dec)*cfringe;
% Area capacitance of word line in the column decoder
carea_cdec_wl=(l_wl_col_dec.*w_wl_col_dec)*carea;
% Diffusion capacitance of word line in the column decoder
cdiff_cdec_wl=4*lambda*(c_diff_nmos.*t_log2_nrows + sum_pow_c.*c_diff_pmos);
%{
Capacitance of the word line wire in the column decoder
= area capacitance + fringing capacitance
%}
cwire_cdec_wl=carea_cdec_wl+cfringe_cdec_wl;
%{
Capacitance of the word line in the column decoder, including
the diffusion capacitance
%}
c_cdec_wl=cwire_cdec_wl + cdiff_cdec_wl;





%------------------------------------------- 
%{
Parasitic resistance and capacitance for signal lines of the
column decoder
%}
% Resistance R_rdec_wl of the word line in the column decoder
R_cdec_ht=R_sheet*(l_ht_col_dec.*w_ht_col_dec);
% Fringing capacitance of word line in the column decoder
cfringe_cdec_ht=(l_ht_col_dec)*cfringe;
% Area capacitance of word line in the column decoder
carea_cdec_ht=(l_ht_col_dec.*w_ht_col_dec)*carea;
% Diffusion capacitance of word line in the column decoder
cdiff_cdec_ht=4*lambda*(c_diff_nmos.*t_log2_nrows + sum_pow_c.*c_diff_pmos);
%{
Capacitance of the word line wire in the column decoder
= area capacitance + fringing capacitance
%}
cwire_cdec_ht=carea_cdec_ht+cfringe_cdec_ht;
%{
Capacitance of the word line in the column decoder, including
the diffusion capacitance
%}
c_cdec_ht=cwire_cdec_ht + cdiff_cdec_ht;

%-------------------------------------------
%{
Developing \pi=models for the row and column decoders

For schematic representation of \pi-models, use R for resistors,
C for capacitance, and * for nodes

Schematic for the parasitics of the row and column decoders:
Use Ri and Ci for the parasitic resistance and capacitance
of the signal lines, Ro and Co for the parasitic resistance
and capacitance of the word line, and R for the input resistance
of the circuit

R = 2^(n-1)*(R_eff_pmos*width + R_eff_pmos) since the
input is connected to 2^(n-1) NMOS transistors and 2^(n-1)*width
PMOS transistors

Note that: powers_r(n)=2^(n-1)

-R--*---Ri--*------Ro---*
    |       |           |
  Ci/2  (Ci/2+Co/2)    Co/2
    |       |           |
   GND     GND         GND

%}

% Value of R, of "input" resistance to decoder
R_ip_rdec=[2.^t_log2_nrows].*(R_nmos+R_pmos*wr(length(t_log2_nrows)))

%-------------------------------------------

% +++++++++++++++++++++++++++++++++++++++++++++++++++
%{
Determine the delay through the \pi-model of the parasitic
resistance and capacitance of the row decoder
%}
tau_rdec=R_ip_rdec.*c_rdec_ht/2;
tau_rdec=tau_rdec+0.5*(R_ip_rdec+R_rdec_ht).*(c_rdec_ht+c_rdec_wl);
tau_rdec=tau_rdec+0.5*(R_ip_rdec+R_rdec_ht+R_rdec_wl).*(c_rdec_wl);

% +++++++++++++++++++++++++++++++++++++++++++++++++++
%{
Determine the delay through the \pi-model of the parasitic
resistance and capacitance of the column decoder
%}
tau_cdec=R_ip_rdec.*c_rdec_ht/2;
tau_cdec=tau_cdec+0.5*(R_ip_rdec+R_cdec_ht).*(c_cdec_ht+c_cdec_wl);
tau_cdec=tau_cdec+(R_ip_rdec+R_cdec_ht+R_cdec_wl).*(c_cdec_wl);
% +++++++++++++++++++++++++++++++++++++++++++++++++++

tau_dec=max(tau_cdec,tau_rdec)
%{
Determine the delay through the \pi-model of the parasitic
resistance and capacitance of the critical path for read
line from the input flip flop to the output flip flop
%}

%{
Number of stages in the SRAM design:
Input flip-flop;
Row/Column decoder;
SRAM block/cell;
Sense amplifer;
Output flip flop
%}
num_stages=5;


%{
Input capacitance of the SRAM is the input capacitance of
its input flip-flop, which has a transmission gate as its
first stage
Assume that the input pads of the chip have no parasitics
%}
c_ip_sram=c_gate_inv*(2+2+1+1)*u
%{
Output capacitance of the SRAM is the output capacitance 
of its output flip-flop, which is the diffusion capacitances
of the inverter at the last stage of the flip-flop
Assume that the output pads of the chip have no parasitics
%}
c_op_sram=c_diff_nmos*(2*u) + c_diff_pmos*(4*u)
%{
Electrical Effort for SRAM design
==>input capacitance of input flip-flop is 
==>the same as the output capactance of
%}
H_sram=c_op_sram/c_ip_sram
% Branching effort = 1 since there is no branching
B_sram=1


% Parasitic for a flip-flop
%%%%%%%FIX THIS
p_ff=0.25*u*(3*(2+2+1+1)+2*(2+1)+(4+2)+(4*1))/(4*lambda)
% Parasitic for a decoder = parasitic for a NOR gate
p_dec=max(prd,pcd)
% Parasitic for a SRAM cell = parasitic for an inverter
p_sram_cell=1
% Parasitic for a sense amplifier = parasitic for an inverter
p_sa=1


%{
Logical effort for the flip-flop, which is a cascade of
inverter-driven transmission gate, inverter, another 
inverter-driven transmission gate, and another inverter

Logical effort for inverter-driven transmission gate = 2
Logical effort for inverter = 1
%}
g_ff=2*1*2*1
% Logical effort for a decoder
g_dec=ones(1,8) %=max(gr,gc)
% Logical effort for SRAM cell = le_sram_cell_r
% Logical effort for sense amp
le_samp_r=le_r_pass_n

%{
Ignore logical effort of column MUX/DEMUX for simplicity,
since its values are negligible
%}

% The loop...
for z=1:length(wr)
    % Logical Effort for SRAM design
    %%%%%%%FIX THIS
    G_sram(z)=g_ff*g_dec(z)*le_sram_cell_r*le_samp_r*g_ff
    
    % Path effort for the SRAM design
    F_sram(z)=G_sram(z)*B_sram*H_sram
    % Parasitics of the SRAM design
    P_sram(z)=p_ff+p_dec(z)+p_sram_cell+p_sa+p_ff
    D_sram(z)=num_stages*(F_sram(z))^(1/num_stages)*tau+P_sram(z)*pinv
    D_sram(z)=D_sram(z)+tau_dec(z);
end
%{
256*128
64*512
128*256
128*256
%}
G_andy(1)=g_ff*(9/4)*le_sram_cell_r*le_samp_r*g_ff
G_andy(2)=g_ff*(9/4)*le_sram_cell_r*le_samp_r*g_ff
G_andy(3)=g_ff*(125/64)*le_sram_cell_r*le_samp_r*g_ff
G_andy(4)=g_ff*(125/64)*le_sram_cell_r*le_samp_r*g_ff
F_andy(1)=G_andy(1)*B_sram*H_sram*(64*5.99*f+154.8288*f)/(5.99*f)
F_andy(2)=G_andy(2)*B_sram*H_sram*(64*91.968*f+154.8288*f)/(91.968*f)
F_andy(3)=G_andy(3)*B_sram*H_sram*(16*3.833*f+154.8288*f)/(3.833*f)
F_andy(4)=G_andy(4)*B_sram*H_sram*(16*33.788*f+154.8288*f)/(33.788*f)
P_andy(1)=p_ff+10+p_sram_cell+p_sa+p_ff
P_andy(2)=p_ff+10+p_sram_cell+p_sa+p_ff
P_andy(3)=p_ff+11+p_sram_cell+p_sa+p_ff
P_andy(4)=p_ff+11+p_sram_cell+p_sa+p_ff
D_andy(1)=(num_stages+4)*(F_andy(1))^(1/(num_stages+4))*tau+P_andy(1)*pinv
D_andy(2)=(num_stages+4)*(F_andy(2))^(1/(num_stages+4))*tau+P_andy(2)*pinv
D_andy(3)=(num_stages+6)*(F_andy(3))^(1/(num_stages+6))*tau+P_andy(3)*pinv
D_andy(4)=(num_stages+6)*(F_andy(4))^(1/(num_stages+6))*tau+P_andy(4)*pinv
% =======================================================

% =======================================================



















% Start working from here!!!!
% FIX THE REST OF THE SCRIPT LATER!!!

%{
The delay for the write data path is determined by
the Elmore delay model; see schematic below


Pictorial representation of circuit:

... ---vVVV^----O----vVVV^------O-----vVVV^-----:
                |               |               |
                |               |               |
                =               =               =
                |               |               |
                GND            GND             GND

                    R_{nmos}         R_{nmos}
                    M9/M10             M7/M8
    R_{wl}  C_{wl+ndiff}    2*C_{ndiff}

Resistance of the NMOS
= R_{nmos} in \ohm*m / (total width of NMOS transistors)
Resistance of the bit and bit_bar lines
= (Number of rows*40+2*34+120)*lambda*4*lambda
{Each bit/bit_bar line cover X number of cells that measure
40 lambda by 40 lambda, and it also spans the length of 3 
transistors; two of this transistors have lengths of 34
lambda, while the other has a length of 120 lambda}
%}


% Resistance of the bit and bit_bar line for write data path
%R_wdp_bl=(nrows*40+2*34+120)*lambda*4*lambda*R_sheet;
% Capacitance of the bit and bit_bar line for write data path
%%%%%%%%%%%%%%C_wdp_bl=()*

% Resistance of the bit and bit_bar line for reading
R_read_bl=(nrows*40+50)*lambda*4*lambda*R_sheet;
% Fringing capacitance of the bit/bit_bar line for reading
cfringe_read_bl=2*cfringe*(nrows*40+50)*lambda;
% Area capacitance of the bit/bit_bar line for reading
carea_read_bl=carea*(nrows*40+50)*lambda*4*lambda;
% C_{wire} = C_{area} + C_{fringing}
cwire_read_bl=carea_read_bl+cfringe_read_bl;
% Diffusion capacitance of the bit/bit_bar line for reading
cdiff_read_bl=(nrows*40*lambda*c_diff_nmos) + 50*lambda*c_diff_pmos;
% Capacitance of the bit and bit_bar line for reading
C_read_bl=cdiff_read_bl+cwire_read_bl;

%{
Since SRAM block is connected to decoder, use its resistance

---R----*--Rrd--*
        |       |
        C       C
        |       |
       GND     GND
R=R_rdec_wl
C=C_read_bl/2;
%}
tau_read=R_rdec_wl.*C_read_bl./2;
tau_read=tau_read+(R_rdec_wl+R_read_bl).*C_read_bl./2;
D_sram
tau_read
D_sram = D_sram + tau_read;


%{
Sum the product of the transistor's length and width
in the row and column decoders
%}
sum_a=[1:8];
sum_a=sum_a-sum_a;
% Number of input bits for the column decoder of the design
num_col_bits=[0 1 2 3 4 5 6 1024]
% Number of output bits for the column decoder of the design
num_nmos_cdec=[2.^num_col_bits]
for i=1:length(wc)
    % Get the area for the pull-down network
    % 2* for Signal and Signal_bar
    % num_col_bits(i) number of transistors per output
    % 2* for the length of the transistors
    % 4* width of the minimum NMOS
    % num_nmos_cdec(i) number of rows/output
    sum_a(i)=sum_a(i)+2*2*num_col_bits(i).*num_nmos_cdec(i)*4;
    % Get the area for the pull-up network
    % 4* width of the minimum transistor that will scale by width
    % num_col_bits(i)*num_nmos_cdec(1:i)*wc(i) is the width of all PMOS
    % = \sum width of all PMOS = 2^(n-1)*width
    % 2* length of a transistor
    sum_a(i)=sum_a(i)+4*num_col_bits(i).*num_nmos_cdec(i).*wc(i)*2;
end
% Number of input bits for the row decoder of the design
num_row_bits=10 - num_col_bits
num_row_bits(8)=0
% Number of output bits for the row decoder of the design
num_nmos_rdec=[2.^num_row_bits ]
for i=1:length(wr)
    % Get the area for the pull-down network
    % 2* for Signal and Signal_bar
    % num_row_bits(i) number of transistors per output
    % 2* for the length of the transistors
    % 4* width of the minimum MOS
    % num_nmos_rdec(i) number of rows/output
    sum_a(i)=sum_a(i)+2*2*num_row_bits(i).*num_nmos_rdec(i)*4;
    % Get the area for the pull-up network
    % 4* width of the minimum transistor that will scale by width
    % num_row_bits(i)*num_nmos_rdec(1:i)*wr(i) is the width of all PMOS
    % = \sum width of all PMOS = 2^(n-1)*width
    % 2* length of a transistor
    sum_a(i)=sum_a(i)+4*num_row_bits(i).*num_nmos_rdec(i).*wr(i)*2;
end


for i=1:length(wr)
%{
Sum the product of the transistor's length and width
in the sense amplifiers and circuit in the word data path
%}  
    % Sense amplifier...
    % 32* for number of columns
    % 2* for minimum length of transistor
    % 5*50 5 MOS of 50 lambda
    % 2*4 2 PMOS of 4 lambda width
    sum_a(i) = sum_a(i) + 32*2*(50*5+4*2);
    % Write data path circuit
    % 32* for number of columns
    % 2* for minimum length of transistor
    % 4*34 4 MOS of 34 lambda width
    % 2*120 2 PMOS of 120 lambda width
    sum_a(i) = sum_a(i) + 32*2*(4*34+120*2);
    % Add the area for the column MUX/DEMUX
    % 32* 32 pass transistors; one for each column
    % 2* Minimum length for a MOS
    % 2*4 2 MOS of minimum width
    sum_a(i) = sum_a(i) + 32*2*(2*4);
    %{
        Area for the FO4 chain connecting the input flip-flops
        to the SRAM block.
    
        IMPORTANT!!!
        If the FO4 inverter chain is added to the SRAM design,
        make sure that I add this variable to sum_a(i)
    
    
        I am assuming that the FO4 inverter needed to improve
        the path delay is approximately the same as that for
        all SRAM architectures under consideration
    
        This assumption is incorrect, but it will simplify
        area estimation
    %}
    % fo4inv_chain_w=238+715+69+207+2+6+2*(55+83)+128+383+238+715;
    % sum_a(i) = sum_a(i) + fo4inv_chain_w;
    %%sum_a(i) = sum_a(i) + num_nmos_rdec(i)*2*8*(238);
    sum_a(i) = sum_a(i) + 8*(80+40);
    
    dff_w=40+20+40+20+40+20+2*(40+20)+40
    
    sum_a(i) = sum_a(i) + 42*dff_w*8;
    
    % Area=6*\sum_{all transistors} + 360*number_of_transistors
    % number_of_transistors = 202572
    sum_a(i) = 6*sum_a(i) + 360*202572;
    %sum_a(i) = sum_a(i) + max(nrows(i),ncols(i)).^2;
end

%{
Since 8 designs are considered for the SRAM block, the details
for the SRAM blocks are placed in a 1*8 array.
%}

%{
for i=1:8
    sum_area(i) = sum_a(i);
end
%}

% Calculate the final area of the SRAM design
[total_area] = area_calc(nrows, ncols, sum_a)


% Determine the area for Andy's designs...
andy_a=zeros(1,4);
for ad=1:4
    % Add the area for the sense amplifier
    andy_a(ad)=andy_a(ad)+ 32*2*(50*5+4*2);
    % Add the area for the write data path
    andy_a(ad)=andy_a(ad)+ 32*2*(4*34+120*2);
    % Add the area for the column MUX/DEMUX
    andy_a(ad)=andy_a(ad)+ 32*2*(2*4);
    % Add the area for the clock buffers...
    andy_a(ad)=andy_a(ad)+ 8*(80+40);
    % Add the area for the D flip-flops...
    dff_w=40+20+40+20+40+20+2*(40+20)+40;
    andy_a(ad)=andy_a(ad)+42*dff_w*8;
    % Obtain an overpessimistic (upper bound) estimate
    % Estimate that there will be 200,000+ transistors
    andy_a(ad)=6*andy_a(ad)+360*200000;
end
cur_f=log((9/4)*12.9*(64*5.7362651*f+154.8288*f)/(5.7362651*f))/log(4);
andy_a(1)=andy_a(1)+4*(1+cur_f+cur_f*cur_f/1.5+cur_f*cur_f*cur_f+cur_f*cur_f*cur_f*cur_f/1.5+cur_f*cur_f*cur_f*cur_f*cur_f);
cur_f=log(512*(9/4)*12.9*(64*91.968*f+154.8288*f)/(91.968*f))/log(4);
andy_a(2)=andy_a(2)+(1+cur_f+cur_f*cur_f/1.5+cur_f*cur_f*cur_f+cur_f*cur_f*cur_f*cur_f/1.5);
acur_f=log(512*(9/4)*12.9*(64*33.788*f+154.8288*f)/(33.788*f))/log(4);
ndy_a(3)=andy_a(3)+(1+cur_f+cur_f*cur_f/1.25+cur_f*cur_f*cur_f+cur_f*cur_f*cur_f*cur_f/1.25+cur_f*cur_f*cur_f*cur_f*cur_f+cur_f*cur_f*cur_f*cur_f*cur_f*cur_f);
cur_f=log(512*(9/4)*12.9*(13*3.833*f+154.8288*f)/(3.833*f))/log(4);
andy_a(4)=andy_a(4)+(1+cur_f+cur_f*cur_f/1.25+cur_f*cur_f*cur_f+cur_f*cur_f*cur_f*cur_f+cur_f*cur_f*cur_f*cur_f*cur_f);
% Determine the area usage of Andy's designs
arows=[256 64 128 128];
acols=[128 256 256 512];
[andy_area] = area_calc(arows, acols, andy_a)
% Add the area and speed for Andy's designs
D_sram = [ D_sram D_andy ]
total_area = [total_area andy_area]



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Determine the resistances and capacitances needed for the
% \pi models for the wiring and interconnects.
%
% Determine the delay of the SRAM design.
% Obtain the worst case critical path for each architecture
% trade-off.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot delay against area of the SRAM designs
%%%plot(total_area,d_decoder,'x')
%%%D_sram
plot(total_area,D_sram,'x')
title('Graph of delay versus area trade-offs for various SRAM designs')
xlabel('area \lambda^2')
ylabel('delay (s)')
disp('End architecture exploration...')
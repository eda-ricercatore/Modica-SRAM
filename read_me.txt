AndyCalcs.m determines the sizes of the FO4 inverter chain that is needed for the word line driver.

area_delay_trade-off indicates the scripts used to semi-automate the calculation/estimation of the delays for the SRAM designs, and plot the delay-area trade-off design points.
estimation_of_tau_pinv is used to estimate the value of \tau and pinv that are used in the method of logical effort to estimate delay.
sizing_gates_in_sense_amplifier contains scripts used to estimate the sizes of the transistors in the sense amplifier.

test_benches contain Nanosim scripts used for power analysis and performance computation of the SRAM design




Explanation of MATLAB scripts
The MATLAB script [sram_design.m] is used to determine the logical
effort for different architectural designs to find a good trade-off
between delay and area, while using delay minimization as the sole
objective of performance optization.

It calls the script [le_pass_transistor.m] to determine the logical
effort of the pass transistor that was intended for use in the 
column MUX/DEMUX, sense amplifier design, and the circuit for the
write data path.

It also calls the scripts [coldecoder.m] and [rowdecoder.m] to 
determine the delay through the column and row decoders for various
SRAM architectures.

The architecture of the SRAM being designed depends on the arrangement
of macros in the SRAM block. The following arrangements of SRAM
blocks are considered:
1024	rows * 32		columns = 1024	\lambda * 1024	\lambda
512		rows * 64		columns = 512	\lambda * 512	\lambda
256		rows * 128 		columns = 256	\lambda * 256	\lambda
128 	rows * 256 		columns = 256	\lambda * 256	\lambda
64 		rows * 512 		columns = 512	\lambda * 512	\lambda
32 		rows * 1024 	columns = 1024	\lambda * 1024	\lambda
16 		rows * 2048 	columns = 2048	\lambda * 2048	\lambda
1 		rows * 32768	columns = 32768	\lambda * 32768	\lambda

The arrangement of SRAM blocks determine the dimensions/sizes of the
row and column decoders. For example, a SRAM block of n rows and m
columns will have the following decoders: a a-to-n row decoder and
a b-to-m column decoder, where (a+b=10 bits) is the number of inputs
to the row and column decoders. From the size/dimension of the decoder,
it can be determined if there are any long wires/interconnects that
exceed 500 \lambda. If so, such wires/interconnects can be modeled as
a \pi-model.

In addition, from using the method of logical effort, the optimum
number of stages for this given path(s) can be determined. If the 
current number of stages in the given path is not the optimum number,
add buffers to the given path. This may subsequently require the 
use differeiation or iteration to size these gates in the given path,
and other circuits in the given path.

If a SRAM macro contains more than 32 SRAM cells in any direction or
dimension, it will need to create \pi-models for modeling the parasitic
capacitances and resistances between them.




Hence, the appropriate leftover space is:
1024 - 32 =	992		==> 992		\lambda * 1024	\lambda = 
512 - 64 = 448		==> 448		\lambda * 512	\lambda = 
256 - 128 = 128		==> 128		\lambda * 256	\lambda = 
256 - 128 = 128		==> 128		\lambda * 256	\lambda = 
512 - 64 = 448		==> 448		\lambda * 512	\lambda = 
1024 - 32 =	992		==> 992		\lambda * 1024	\lambda = 
2048 - 32 = 2032	==> 2032	\lambda * 2048	\lambda = 
32768 - 1 = 32767	==> 32767	\lambda * 32768	\lambda = 










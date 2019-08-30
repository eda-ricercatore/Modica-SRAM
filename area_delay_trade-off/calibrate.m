%{
This is written by Zhiyang Ong and Andrew Mattheisen for
	EE 577B, SRAM Project - Part 2, Questions 1 and 2

Calibrating the value of \tau, 
%}

rising_delay=[8.498e-11 1.130e-10 1.391e-10 1.634e-10];
%rising_delay=[rising_delay 8.672e-10];
%rising_delay=[rising_delay ]
falling_delay=[7.476e-11 1.094e-10 1.383e-10 1.663e-10];
%falling_delay=[falling_delay 8.048e-10];
%falling_delay=[falling_delay ]
capacitance=(2*0.95840*1e-15) * [2 4 6 8];
%h=capacitance./((2.396*1e-15/1e-6)*4*0.1*1e-6*(16))
h=capacitance./((2.396*1e-15/1e-6)*4*0.1*1e-6*(2))
avg_delay=rising_delay+falling_delay;
avg_delay = 0.5 * avg_delay




[p,err]=polyfit(h,avg_delay,1)
d_fit=polyval(p,h,err);




plot(h,avg_delay,'-^',h,d_fit,'--x')
title('Plot of delay (s) against electrical effort, h')
xlabel('electrical effort, h (unitless)')
ylabel('delay (s)')
text(5,10e-11,'delay = (14.122 ps) \cdot h + (53.045 ps)')

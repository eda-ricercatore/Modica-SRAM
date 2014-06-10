% This file calculates the delay of the unit inverter for the tsmc_018
% technology file provided for the SRAM project.

lambda=.1*(10^-6);
Cgate=2.396*(10^-15)/(10^-6);
cin=(4+12)*lambda*4*Cgate;
%{
from simulation:
 .param cload = 30f *h=2                                                             
   crisingdelay=  8.4985E-11  targ=  2.1627E-09   trig=  2.0777E-09
   cfallingdelay=  7.4758E-11  targ=  4.1541E-09   trig=  4.0794E-09

 .param cload = 60f *h=4
   crisingdelay=  1.1295E-10  targ=  2.1905E-09   trig=  2.0776E-09
   cfallingdelay=  1.0937E-10  targ=  4.1866E-09   trig=  4.0772E-09
 
 .param cload = 90f *h=6
   crisingdelay=  1.3906E-10  targ=  2.2166E-09   trig=  2.0775E-09
   cfallingdelay=  1.3826E-10  targ=  4.2153E-09   trig=  4.0770E-09

 .param cload = 120f *h=8
   crisingdelay=  1.6338E-10  targ=  2.2417E-09   trig=  2.0784E-09
   cfallingdelay=  1.6634E-10  targ=  4.2433E-09   trig=  4.0770E-09

FROM SIMULATION WITH BETTER NETLIST>>>
h=2
.param load_size = 1                                                                
   crisingdelay=  7.4254E-11  targ=  4.1888E-09   trig=  4.1145E-09
   cfallingdelay=  6.6566E-11  targ=  2.1862E-09   trig=  2.1197E-09
h=4
 .param load_size = 2
   crisingdelay=  9.0220E-11  targ=  4.2047E-09   trig=  4.1145E-09
   cfallingdelay=  8.2377E-11  targ=  2.2020E-09   trig=  2.1196E-09
h=6
 .param load_size = 3
   crisingdelay=  1.0602E-10  targ=  4.2205E-09   trig=  4.1145E-09
   cfallingdelay=  9.7521E-11  targ=  2.2171E-09   trig=  2.1196E-09
h=8
 .param load_size = 4
   crisingdelay=  1.2190E-10  targ=  4.2364E-09   trig=  4.1145E-09
   cfallingdelay=  1.1237E-10  targ=  2.2319E-09   trig=  2.1196E-09

%}
%For 2nd extra run
%Loadcap=[30e-15,   60e-15,     90e-15,    120e-15];
%dRise=[8.4985E-11, 1.1295E-10, 1.3906E-10, 1.6338E-10];
%dFall=[7.4758E-11, 1.0937E-10, 1.3826E-10, 1.6634E-10];

%Loadcap=[30e-15,   60e-15,     90e-15,    120e-15];
h=[2,4,6,8]
dRise=[7.4254E-11, 9.0220E-11, 1.0602E-10, 1.2190E-10];
dFall=[6.6566E-11, 8.2377E-11, 9.7521E-11, 1.1237E-10];




for i=1:length(dRise)
    %calculate values for delay
    dAve(i)=(dRise(i)+dFall(i))/2;

    %calculate values for h
    %h(i)=Loadcap(i)/cin;
end


% Calculate fit parameters
[p,ErrorEst] = polyfit(h,dAve,1);
% Evaluate the fit
d_fit = polyval(p,h,ErrorEst);
% Plot the data and the fit
figure(1)
plot(h,d_fit,'-',h,dAve,'+');
title('SRAM PROJECT - Determination of Tao by simulation')
xlabel('H of Inverter')
ylabel('Delay of inverter, (sec)')
% use text(<x-coord>, <y-coord>,'the text you wish to display')
p 

  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for fit_line                  %
% slope is 14.34ps              %
% y-intercept is 53.65ps        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

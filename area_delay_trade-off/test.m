powers_all=[2.^([1:10]-1)]
for i=1:length(powers_all)
    sum_powers(i)=sum(powers_all(1:i))
end
disp('Next')
sum_powers
%{
if 3>4
    disp('3 < 4')
elseif 7<8
    disp('7 is LESS than 8')
else
    disp('Err msg')
end







if 5~=12
   disp('5 is not equal to 12') 
end


if 3<4
    disp('pre test')
    error('Error found')
end
disp('Error was not processed')
%}
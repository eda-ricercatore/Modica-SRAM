% Synopsis:
% Determine the logical effort of a gate, with due consideration
% given to various CMOS technologies
%
% @param gate_type is the type of the digital logic gate
% @param cir_fam is the type of the circuit family in use
% 
function [le] = logical_effort(gate_type,cir_fam)

%{
Complementary CMOS logic gates [Weste & Harris, pp. 11]
[Rabaey & et. al., pp. 237-263]
%}
switch gate_type
    case 'nand'
        disp('nand gate found');
	case 'nor'
        disp('nor gate found');
	case 'and'
        disp('and gate found');
	case 'or'
        disp('or gate found');
	case 'not'
        disp('not gate found');
    otherwise
        % Assume that the gate 
        disp('Invalid gate type');
end


%{
Ratioed circuits: Pseudo-NMOS [Weste & Harris, pp. 327-331]
[Rabaey & et. al., pp. 263-268]
%}


%{
domino logic [Rabaey & et. al., pp. 285-301]
[Weste & Harris, pp. 327-331]
%}





% Consider DPL, CMOSTG, and CVSL for robustness
% and domino and dynamic gates for speed; also try NTL and LEAP
% See summary on pp.369




% Ka chong










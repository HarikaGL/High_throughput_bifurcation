function p = tr_move(p, type, a, b)
%TR_MOVE  Apply one transition move to the parameter vector p (1 x 14).
%   'setx', a : x = log10(k2 KTOT/(k8 PTOT)) = a, by changing k8
%   'sety', a : y = log10(k4 KTOT/(k6 PTOT)) = a, by changing k6
%   'setyk4',a: y = a, by changing k4
%   'bp',   a : k4 = k8 = k2*10^a (keeps the symmetry of Iterations 3-4)
%   'km',   a : log10(KM1/KM5) = a, by changing k5
%   'mul', a, b : p(a) = p(a)*10^b (one parameter, a = 1..14)
%   'pair', a, b: p(a) = p(a)*10^b for a pair of indices a = [i j]
KT = p(1); PT = p(2);
switch type
    case 'setx',   p(14) = p(5)*KT/PT*10^(-a);
    case 'sety',   p(11) = p(8)*KT/PT*10^(-a);
    case 'setyk4', p(8)  = p(11)*PT/KT*10^(a);
    case 'bp',     p(8) = p(5)*10^a; p(14) = p(8);
    case 'km',     KM1 = (p(4) + p(5))/p(3); p(9) = (p(10) + p(11))/(KM1*10^(-a));
    case 'mul',    p(a) = p(a)*10^b;
    case 'pair',   p(a) = p(a)*10^b;
    otherwise, error('tr_move: unknown move %s', type);
end
end

function [p, val] = tr_sweep_point(p0, typ, pj, v)
%TR_SWEEP_POINT  Parameter set at sweep value v of a two-class sweep (TR_03), and the value
%   of the swept quantity: 'x' / 'y' set to v through parameter pj; 'bp': log10(k4/k2) = v with
%   k4 = k8; 'km': log10(KM1/KM5) = v through k5; 'mul': p(pj)*10^(2v), value log10(k3/k7).
p = p0; KT = p(1); PT = p(2);
switch typ
    case 'x'
        x0 = log10(p(5)*KT/(p(14)*PT)); d = v - x0;
        if any(pj == [14 2]), p(pj) = p(pj)*10^(-d); else, p(pj) = p(pj)*10^d; end
        val = v;
    case 'y'
        y0 = log10(p(8)*KT/(p(11)*PT)); d = v - y0;
        if any(pj == [11 2]), p(pj) = p(pj)*10^(-d); else, p(pj) = p(pj)*10^d; end
        val = v;
    case 'bp'
        p(8) = p(5)*10^v; p(14) = p(8); val = v;
    case 'km'
        KM1 = (p(4) + p(5))/p(3); p(9) = (p(10) + p(11))/(KM1*10^(-v)); val = v;
    case 'mul'
        p(pj) = p(pj)*10^(2*v); val = log10(p(6)/p(12));   % shown as log10(k3/k7)
end
end

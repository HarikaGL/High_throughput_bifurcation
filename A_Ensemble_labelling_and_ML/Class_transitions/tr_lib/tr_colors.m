function col = tr_colors(what)
%TR_COLORS  Colours used in every transition figure.
%   tr_colors('main'|'detached'|'sym'|'offsym') branch colours; tr_colors(k) class k colour
%   (class_info_v3 colours, Class 7 orange, Class 8 magenta as in Fig 11).
if ischar(what)
    switch what
        case 'main',     col = [0.16 0.44 0.73];
        case 'detached', col = [0.75 0.22 0.17];
        case 'sym',      col = [0.33 0.33 0.33];
        case 'offsym',   col = [0.88 0.51 0.16];
        otherwise,       col = [0.4 0.4 0.4];
    end
else
    CI = class_info_v3(); col = 0.85*CI.colors(what, :);
    if what == 7, col = [0.93 0.62 0.05]; end
    if what == 8, col = [0.80 0.10 0.80]; end
end
end

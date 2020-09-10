% Only if flag_location_legend = 1 the function deals with the position of the
% legend
function set_fonts_lines(...
    ha, ...
    flag_location_legend)
line     = 2;
marker   = 10;
font     = 20;
fontaxes = 18;
set(0, 'ShowHiddenHandles', 'on')
hl = findobj(ha, 'Type', 'line', '-not', 'LineStyle', 'none');
set(hl,  'Linewidth', line);
set(hl,  'Markersize', marker)
hm = findobj(ha, 'Type', 'line', '-and', 'LineStyle', 'none');
set(hm, 'Markersize', marker);
set(hm, 'Linewidth', line);
he = findobj(ha, 'Type', 'errorbar');
set(he, 'Markersize', marker);
set(he, 'Linewidth', line);
set(ha, 'Fontsize', fontaxes);
if nargin >= 3
    if flag_location_legend == 1
        hLeg = legend(ha);
        set(hLeg, 'Location', 'Best');
    end
end
% hFig = gcf;
ht   = findobj(ha, 'Type', 'text');
set(ht, 'Fontsize', font);
end
function [handles] = surf_plot(trsurfs,x,y,z)
%SURF_PLOT plots the st-form surfaces in cell array trsurfs, as well as the
%original points x, y, z
%
%   INPUTS
%   trsurfs  1-by-numdivs cell array containing trimmed stform surfaces
%   x        1-by-N vector
%   y        1-by-N vector
%   z        1-by-N vector
%
%   OUTPUTS
%   handles  1-by-numdivs vector containing handles to the plotted surfaces

figure

%plot trimmed surfaces
for i=1:length(trsurfs)
    fnplt(trsurfs{i});
    hold on;
end
shading interp;
handles=get(gca,'children');

%overlay original data points
plot3(x,y,z,'wo','markerfacecolor','k');

end
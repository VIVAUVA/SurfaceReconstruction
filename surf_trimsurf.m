function [outsurf] = surf_trimsurf(insurf,xextent,yextent)
%SURF_TRIMSURF discards all values of insurf which lie outside of xextent 
%and yextent. returns outsurf in the same format as insurf
%
%   INPUTS
%   insurf stform surface
%   xextent  1-by-2 vector in the form [lowxbound upxbound]
%   yextent  1-by-2 vector in the form [lowybound upybound]
%
%   OUTPUTS
%   outsurf  stform surface (with overlap data removed)

interv=cell(1,2);
interv{1}=xextent;
interv{2}=yextent;
insurf.interv=interv;
outsurf=insurf;

end
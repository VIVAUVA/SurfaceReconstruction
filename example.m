%Surface Reconstruction Example

load data;
%data contains x,y,z coordinates for a sparse cloud of 1407 points

%calculate mosaic of thin-plate, smothing spline surfaces
trsurfs=surf_main(x,y,z);

%plot resulting surfaces along with original data
surf_plot(trsurfs,x,y,z);

%note: error between surfaces is set to present no interficial error when
%sampled to create a 256-level grayscale image; this parameter can be
%modified in surf_recalc.
function [trsurfs] = surf_main(x,y,z)
%SURF_MAIN is the main function of the surface construction tool. This
%function accepts 3-dimensional, sparse data information (x,y,z), breaks it
%into a mosaic, calculates thin-plate-spline surfaces for each mosaic
%piece, and returns a stitched meta-surface.
%
%   INPUTS
%   x  1-by-N vector
%   y  1-by-N vector
%   z  1-by-N vector
%
%   OUTPUTS
%   trsurf  1-by-N cell array, each cell contains part of the final surface
%           mosaic in stform
%
%
%
%

%check inputs



%--------------------------------------------------------------------------
%- mosaic
%--------------------------------------------------------------------------
%calculate mosaic
[segcoords xdivs ydivs trxext tryext]=surf_mosaic(x,y,z);
numdivs=xdivs*ydivs; %total number of mosaic pieces

%--------------------------------------------------------------------------
%- surface reconstruction (part I)
%--------------------------------------------------------------------------
%perform initial surface calculation (mosaic pieces overlap each other 25%)
smoothness=1; %interpolation factor (0 to 1), 1 results in zero error at
%data locations

surfs=cell(1,numdivs); %overlapping surface for each mosaic piece
trsurfs=cell(1,numdivs); %non-overlapping surface for each piece

%handle wait bar
start=tic; %timer
wait=waitbar(0,'Calculating Surfs');
barobj=findobj(wait,'type','patch');
set(barobj,'edgecolor',[63/256 79/256 186/256],'facecolor',...
    [63/256 79/256 186/256])

%calculate each surface
for i=1:numdivs
    %format data & data locations
    points=[segcoords{1,i};segcoords{2,i}];
    zdata=[segcoords{3,i}];
    
    %update wait bar
    waitbar(i/numdivs,wait,sprintf(['Calculating surfs, #%i of %i, '...
        '%i datapoints\nTotal time elapsed: %i seconds'],...
        i,numdivs,length(zdata),round(toc(start))))
    
   %thin-plate-spline surface generation
    if gt(length(zdata),3)
        surf=surf_tpaps(points,zdata,smoothness);
        surfs{i}=surf;
    else %handle empty mosaic regions
        surfs{i}=[]; %null-matrix represents a flat surface (no data)
    end
end
close(wait) %close wait bar

%--------------------------------------------------------------------------
%- surface reconstruction (part 2)
%--------------------------------------------------------------------------
%check interficial points between trimmed surfaces, if they are above a
%certain error threshold, recalculate pertinent surfaces by increasing
%their overlap and incorporating more data points (higher accuracy)

surfs=surf_recalc(x,y,z,surfs,xdivs,ydivs,start,trxext,tryext);

%calculate trimmed surfaces
for i=1:numdivs
    trsurfs{i}=surf_trimsurf(surfs{i},[trxext(1,i) trxext(2,i)],...
        [tryext(1,i) tryext(2,i)]);
end

end
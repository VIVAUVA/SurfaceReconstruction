function [segcoords xdivs ydivs trxext tryext] = surf_mosaic(x,y,z)
%SURF_MOSAIC uses (x,y) coordinates to subdivide spatial data into 
%manageable chunks--i.e. equally-sized, approximately square regions which
%contain at most 500 points (and ideally about 300 points) for optimal
%surface construction
%
%   INPUTS
%   x  1-by-N vector
%   y  1-by-N vector
%   z  1-by-N vector
%
%   OUTPUTS
%   segcoords  3-by-N cell array, nth column corresponds to nth mosaic
%              piece, columns 1,2,3 contain x,y,z coordinates respectively
%   xdivs      double, number of columns in mosaic
%   ydivs      double, number of rows in mosaic
%   trxext     2-by-N array, nth column corresponds to nth mosaic piece,
%              row 1 contains lower x bound for piece, row 2 upper x bound
%   tryext     2-by-N array, same as trxext but row 1 contains lower y
%              bound, row 2 upper y bound


%handle case where x,y coordinates do not begin at origin
xoffset=min(x);
x=x-min(x);
yoffset=min(y);
y=y-min(y);

%calculate approximate proportion of spatial distribution
xmax=max(x);
ymax=max(y);
[nprop dprop]=rat(ymax./xmax,.05);

%want about 300 points in each mosaic piece
%to ensure proper proportion of output image, number of mosaic pieces
%should be nprop*dprop*n^2
%iteratively calculate n needed
rightsize=0; %loop flag
numiter=1; %loop counter
while ~rightsize
    if eq(exist('piecedensity'),1)&&lt(max(piecedensity),450)
        %use previous value of numdivs (we overshot), flag rightsize
        numdivs=nprop.*dprop.*(numiter-2).^2;
        rightsize=1;
    else
        %try new value of numdivs (# of mosaic pieces)
        numdivs=nprop.*dprop.*numiter.^2;
    end
    
    %want ydivs/xdivs to be about ymax./xmax (via proportion of data locations)
    xdivs=round(sqrt((xmax./ymax)*(numdivs))); %num columns
    ydivs=round(sqrt((ymax./xmax)*(numdivs))); %num rows
    
    %mosaic pieces will overlap to ensure accurate construction
    %find points in each mosaic row (starting from top)
    overlap=.25; %experimentally determined overlap factor
    rsegs=cell(1,ydivs); %logical mask for row points
    for j=ydivs:-1:1
        height=ymax/ydivs; %height of each row
        lowbound=ymax*((j-1)/ydivs)-height*overlap; %lower bound for each row
        upbound=ymax*(j/ydivs)+height*overlap; %upper bound for each row
        %find mask for row points
        rsegs{ydivs+1-j}=y<upbound&y>lowbound;
    end
    
    %find points in each mosaic column (starting from left)
    csegs=cell(1,xdivs); %logical mask for column points
    for i=1:xdivs
        width=xmax/xdivs; %width of each column
        lowbound=xmax*((i-1)/xdivs)-width*overlap; %lower bound
        upbound=xmax*(i/xdivs)+width*overlap; %upper bound
        %find mask for column points
        csegs{i}=x<upbound&x>lowbound;
    end
    
    %find points in each mosaic piece (using row and column masks)
    segindices=cell(1,numdivs); %point indices for each overlapping piece
    counter=1;
    for j=1:ydivs
        for i=1:xdivs
            segindices{counter}=find(csegs{i}&rsegs{j});
            counter=counter+1;
        end
    end
    
    segcoords=cell(3,numdivs); %coordinate information for each piece
    piecedensity=zeros(1,numdivs); %number of points in each piece
    for i=1:numdivs
        currentinds=segindices{i};
        segcoords{1,i}=x(currentinds)+xoffset;
        segcoords{2,i}=y(currentinds)+yoffset;
        segcoords{3,i}=z(currentinds);
        piecedensity(i)=length(currentinds);
    end
    
    %check for optimal mean of size of displacement vectors
    if lt(max(piecedensity),500)
        rightsize=1;
    end
    
    numiter=numiter+1;
    
end

%find non-overlapping bounds for mosaic
width=xmax/xdivs; %column width
height=ymax/ydivs; %row height
xlow=0:width:(xmax-width); %column low bounds
xhigh=width:width:xmax; %column high bounds
ylow=(ymax-height):-height:0; %row low bounds
yhigh=ymax:-height:height; %row high bounds
%manipulate bounds vectors to produce extents for each mosaic piece
ylow=repmat(ylow,xdivs,1);
yhigh=repmat(yhigh,xdivs,1);
tryext=[ylow(:)';yhigh(:)']+yoffset; %y-coordinates for piece boundaries
trxext=repmat([xlow;xhigh],1,ydivs)+xoffset; %x-coordinates for piece boundaries

end


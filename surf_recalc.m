function [surfs]=surf_recalc(x,y,z,surfs,xdivs,ydivs,start,trxext,tryext)
%SURF_RECALC examines interfacial errors between surfaces in surfs and
%recalculates to greater accuracy if the error is too large
%
%   INPUTS
%   x      1-by-N vector
%   y      1-by-N vector
%   z      1-by-N vector
%   surfs  1-by-N cell array, each cell contains an stform surface
%   xdivs  double, # of columns in mosaic
%   ydivs  double, # of rows in mosaic
%   start  double, system time at start of surface calculations
%   trxext 2-by-N array, nth column corresponds to nth mosaic piece,
%          row 1 contains lower x bound for piece, row 2 upper x bound
%   tryext 2-by-N array, same as trxext but row 1 contains lower y
%          bound, row 2 upper y bound 
%   
%   OUTPUTS
%   surfs  1-by-N cell array, same format as input
%   

numdivs=xdivs*ydivs; %total # of mosaic pieces

%calculate error threshold--chosen for later grayscale image conversion
errthresh=(max(z)-min(z))./256-.01;

%determine resolution at which to surf_check each edge (heuristic)
checkres=round(4*sqrt(length(x)/numdivs));


%handle wait bar for error checking/surface recalculation
wait=waitbar(0,'Recalculating Surfs');
barobj=findobj(wait,'type','patch');
set(barobj,'edgecolor',[230/256 98/256 19/256],'facecolor',...
    [230/256 98/256 19/256]) %change waitbar color

overlapinc=.15; %overlap increase increment
recalced=zeros(5,numdivs); %first row indicates flag for recalculation,
                           %rows 2-5 store modified overlap criteria
%initialize recalc flags (only null surfaces have 0 flag)
for i=1:numdivs
    if isempty(surfs{i})
        recalced(1,i)=0;
    else
        recalced(1,i)=1;
    end
end
recalced(2:5,:)=.25; %initialize all overlaps to 25%

errors=zeros(4,numdivs); %stores interfacial errors at each edge

count=1; %iteration counter
overthresh=.75; %overlap threshold
recheck=find(recalced(1,:)); %indices of pieces to recheck
notrecalced=zeros(length(recheck)); %flags for next iteration

%Measure error at edges of each mosaic piece, recalculate surf with larger
%overlap if necessary (increase by overlapinc in pertinent directions)
%repeat until below error threshold or if any overlap exceeds the overlap
%threshold (or for a maximum of 15 iterations)
while ~eq(length(recheck),0)&&~eq(count,15)
    for q=1:length(recheck) %surf_check for null (flat) surface
        i=recheck(q); %select next index to recheck
        %update wait bar
        waitbar(q/length(recheck),wait,sprintf(['Checking surf #%i of %i for'...
            ' interficial errors, pass #%i\nTotal time elapsed: %i seconds'],...
            q,length(recheck),count,round(toc(start))))
        
        %extract extent info for this mosaic piece
        ledge=trxext(1,i);
        redge=trxext(2,i);
        aedge=tryext(2,i);
        bedge=tryext(1,i);
        
        %initialize directionalized overlap criteria
        boverlap=recalced(2,i);
        aoverlap=recalced(3,i);
        loverlap=recalced(4,i);
        roverlap=recalced(5,i);
        
        %initialize flag for recalculation
        recalc=0;
        
        %surf_check mosaic-position-dependent conditions
        if i<=(xdivs*(ydivs-1)) %surf_check if not on bottom row
            %calculate bottom-interface error
            berror=surf_check.below([surfs{i} surfs{i+xdivs}],...
                [ledge redge],bedge,checkres);
        else %on bottom row, so no pieces below
            berror=0; 
        end
        if i>=(xdivs+1) %surf_check if not on top row
            %calculate top-interface error
            aerror=surf_check.above([surfs{i} surfs{i-xdivs}],...
                [ledge redge],aedge,checkres);
        else %on top row
            aerror=0;
        end
        if mod(i,xdivs)~=1 %surf_check if not on left extremity
            %calculate left-interface error
            lerror=surf_check.left([surfs{i} surfs{i-1}],...
                [bedge aedge],ledge,checkres);
        else %on left column
            lerror=0;
        end
        if mod(i,xdivs)~=0 %surf_check if not on right extremity
            %calculate right-interface error
            rerror=surf_check.right([surfs{i} surfs{i+1}],...
                [bedge aedge],redge,checkres);
        else %on right column
            rerror=0;
        end
        
        %surf_check for recalculation
        if gt(berror,errthresh)
            recalc=1;
            %increase overlap factor in bottom direction
            boverlap=boverlap+overlapinc; 
        end
        if gt(aerror,errthresh)
            recalc=1;
            %increase overlap factor in top direction
            aoverlap=aoverlap+overlapinc;
        end
        if gt(lerror,errthresh)
            recalc=1;
            %increase overlap factor in left direction
            loverlap=loverlap+overlapinc;
        end
        if gt(rerror,errthresh)
            recalc=1;
            %increase overlap factor in right direction
            roverlap=roverlap+overlapinc;
        end
        
        
        if recalc==1 %recalculate surf with new bounds
            recalced(1,i)=1;
            
            %make overlap conditions proportionally square to prevent
            %skewing of mosaic overlap
            if ~eq(aoverlap+boverlap,loverlap+roverlap)
                if eq(max(aoverlap,boverlap),aoverlap)...
                        &&eq(max(loverlap,roverlap),roverlap)
                    aoverlap=max(aoverlap,roverlap);
                    roverlap=aoverlap;
                    boverlap=max(boverlap,loverlap);
                    loverlap=boverlap;
                elseif eq(max(aoverlap,boverlap),aoverlap)...
                        &&eq(max(loverlap,roverlap),loverlap)
                    aoverlap=max(aoverlap,loverlap);
                    loverlap=aoverlap;
                    boverlap=max(boverlap,roverlap);
                    roverlap=boverlap;
                elseif eq(max(aoverlap,boverlap),boverlap)...
                        &&eq(max(loverlap,roverlap),roverlap)
                    boverlap=max(boverlap,roverlap);
                    roverlap=boverlap;
                    aoverlap=max(aoverlap,loverlap);
                    loverlap=aoverlap;
                else
                    boverlap=max(boverlap,loverlap);
                    loverlap=boverlap;
                    aoverlap=max(aoverlap,roverlap);
                    roverlap=aoverlap;
                end
            end
            
            %store modified overlap info
            recalced(2:5,i)=[boverlap;aoverlap;loverlap;roverlap];
            %store error info
            errors(:,i)=[berror;aerror;lerror;rerror];
            
            %find points within new overlapping mosaic piece
            width=redge-ledge;
            height=aedge-bedge;
            newindices=find(lt(x,redge+roverlap*width)...
                &gt(x,ledge-loverlap*width)...
                &lt(y,aedge+aoverlap*height)...
                &gt(y,bedge-boverlap*height));
            newpoints=[x(newindices);y(newindices)];
            newvals=z(newindices);
            
            %recalculate overlapping surface with new points
            surfs{i}=surf_tpaps(newpoints,newvals,1);
        else
            %mark for removal from recalculation queue
            notrecalced(q)=1;
        end
    end
    %remove relevant entries from recalc queue
    recheck(find(notrecalced))=[];
    %re-initialize notrecalced flags
    notrecalced=zeros(length(recheck));
    count=count+1;
end
close(wait)


end


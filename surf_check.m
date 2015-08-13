%SURF_CHECK accepts two surfaces and calculates the interficial error
%between them
%
%
%   INPUTS
%   surfs     1-by-2 matrix containing two st-form surfaces
%   ybounds   1-by-2 matrix containing boundary of interface to consider,
%             only used for checking left and right surface interfaces
%   xbounds   1-by-2 matrix containing boundary of interface to consider,
%             only used for checking top and bottom surface interfaces
%   xpos      x-coordinate position reference
%   ypos      y-coordinate position reference
%   checkres  scalar specifying the number of sites at an interficial edge
%             to be checked
%   
%   OUTPUTS
%   maxerror  scalar, maximum error between surfaces of interest

classdef surf_check
    methods (Static)
        function [maxerror] = right(surfs, ybounds, xpos, checkres)
            %evaluate the boundary between surf and rsurf (right surf), return the
            %maximum error between the two. Also handle case where rsurf is
            %empty
            surf=surfs(1);
            try rsurf=surfs(2);
            catch exception %rsurf empty, assign arbitrarily large error and exit
                maxerror=1;
                return
            end
            %determine data sites to check
            lpoints=zeros(2,checkres);
            for p=1:checkres
                lpoints(2,p)=(p/checkres).*(ybounds(2)-ybounds(1))+ybounds(1);
            end
            rpoints=lpoints;
            %slightly shift sites to lie on corresponding surfaces
            lpoints(1,:)=xpos-.00000000001;
            rpoints(1,:)=xpos+.00000000001;
            %calculate the maximum error
            maxerror=max(abs(fnval(surf,lpoints)-fnval(rsurf,rpoints)));
        end
        
        function [maxerror] = left(surfs, ybounds, xpos, checkres)
            %evaluate the boundary between surf and lsurf (left surf), return the
            %max error between the two. Also handle case where lsurf is
            %empty--return -1 for error
            surf=surfs(1);
            try lsurf=surfs(2);
            catch exception %lsurf empty, assign arbitrarily large error and exit
                maxerror=1;
                return
            end
            %determine data sites to check
            lpoints=zeros(2,checkres);
            for p=1:checkres
                lpoints(2,p)=(p/checkres).*(ybounds(2)-ybounds(1))+ybounds(1);
            end
            rpoints=lpoints;
            %slightly shift sites to lie on corresponding surfaces
            lpoints(1,:)=xpos-.00000000001;
            rpoints(1,:)=xpos+.00000000001;
            %calculate the maximum error
            maxerror=max(abs(fnval(lsurf,lpoints)-fnval(surf,rpoints)));
        end
        
        function [maxerror] = below(surfs, xbounds, ypos, checkres)
            %evaluate the boundary between surf and bsurf (bottom surf), return
            %the max error between the two. Also handle case where bsurf is
            %empty
            surf=surfs(1);
            try bsurf=surfs(2);
            catch exception %bsurf empty, assign arbitrarily large error and exit
                maxerror=1;
                return
            end
            %determine data sites to check
            tpoints=zeros(2,checkres);
            for p=1:checkres
                tpoints(1,p)=(p/checkres).*(xbounds(2)-xbounds(1))+xbounds(1);
            end
            bpoints=tpoints;
            %slightly shift sites to lie on corresponding surfaces
            bpoints(2,:)=ypos-.00000000001;
            tpoints(2,:)=ypos-.00000000001;
            %calculate the maximum error
            maxerror=max(abs(fnval(bsurf,bpoints)-fnval(surf,tpoints)));
        end
        
        function [maxerror] = above(surfs, xbounds, ypos, checkres)
            %evaluate the boundary between surf and asurf (above surf),
            %return the max error between the two. Also handle case where
            %asurf is empty
            surf=surfs(1);
            try asurf=surfs(2);
            catch exception %asurf empty, assign arbitarily large error and exit
                maxerror=1;
                return
            end
            %determine data sites to check
            tpoints=zeros(2,checkres);
            for p=1:checkres
                tpoints(1,p)=(p/checkres).*(xbounds(2)-xbounds(1))+xbounds(1);
            end
            bpoints=tpoints;
            %slightly shift sites to lie on corresponding surfaces
            bpoints(2,:)=ypos-.00000000001;
            tpoints(2,:)=ypos-.00000000001;
            %calculate the maximum error
            maxerror=max(abs(fnval(surf,bpoints)-fnval(asurf,tpoints)));
        end
    end
end
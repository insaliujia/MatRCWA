classdef Triangle < PatternShape
    % This class is the subclass of PatternShape, it is used to further
    % define the shape of the pattern
    % Input example: Tri = Triangle([1,1],3,[2:4],'Si',[1,1]);
    % Tri = Triangle([centerx,centery],sidelength,[layers],material,er and ur);
    % centerx and centery the center of the rectangle shape device
    properties
        name = 'Triangle'
        SideLength
        material
        er
        ur
    end
    
    methods
        function Tri = Triangle(center,SideLength,nlayer,material)
            Tri = Tri@PatternShape(center,nlayer);
            if numel(SideLength) == 1
                Tri.SideLength = SideLength;
            else
                error('Check the radius input')
            end
                Tri.er = material.er;
                Tri.ur = material.ur;
                Tri.material = material.MaterialName;
        end
        
        function BuildPattern(Tri,Dev,varargin)
             % Purpose: make equilateral triangle pattern on the device
             % Input: center--the center of the rectangle,SideLength--side
             % length,nlayer--in which layer will the rectangle
             % be, filler--the material that will fill the hole.
             % Inputformat: Dev =Dev.RectanglePatternDevice([0.2,0.25],0.1,[3:5],[2,1]);                           
            Nx = Dev.idimension(1);
            Ny = Dev.idimension(2);
            Lx = Dev.xydimension(1);
            Ly = Dev.xydimension(2);
%             ER = Dev.ER;
%             UR = Dev.UR;
            dx = Lx/Nx;
            dy = Ly/Ny;
            SideLen = Tri.SideLength;
            h = 0.5*sqrt(3)*SideLen; 
            x0 = Tri.center(1);
            y0 = Tri.center(2);
            if abs((Lx/2-x0))+SideLen/2 > Lx/2 || h*2/3-(Ly/2-y0)>Ly/2 ...
                    ... || h/3 + (Ly/2-y0)>Ly/2
                error('The rectangle is too large!')
            end            
            nxm = floor((Lx/2-x0)/dx);
            nym = floor((Ly/2-y0)/dy);
            ny = round(h/dy); 
            ny1 = round((Ny - ny)/2)-nym; 
            ny2 = ny1 + ny - 1;
            if nargin == 2
                for ny = ny1 : ny2 
                    f = (ny - ny1)/(ny2 - ny1); 
                    nx = round(f*SideLen/Lx*Nx); 
                    nx1 = 1 + floor((Nx - nx)/2)-nxm; 
                    nx2 = nx1 + nx; 
                    Dev.ER(nx1:nx2,ny,Tri.nlayer) = Tri.er;
                    Dev.UR(nx1:nx2,ny,Tri.nlayer) = Tri.ur;
                end
            elseif nargin == 3
                for ny = ny1 : ny2 
                    f = (ny - ny1)/(ny2 - ny1); 
                    nx = round(f*SideLen/Lx*Nx); 
                    nx1 = 1 + floor((Nx - nx)/2)-nxm; 
                    nx2 = nx1 + nx; 
                    Dev.ER(nx1:nx2,ny,Tri.nlayer) = Tri.er(varargin{1},2);
                    Dev.UR(nx1:nx2,ny,Tri.nlayer) = Tri.ur(varargin{1},2);
                end
            else
                error('Check input number')
            end
                
        end
      
    end
    
end
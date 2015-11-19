classdef Rectangle < PatternShape
    % This class is the subclass of PatternShape, it is used to further
    % define the shape of the pattern
    % Input example: rec = Rectangle([1,1],[3,3],[2:4],'Si',[1,1]);
    properties
        name = 'Rectangle'
        rectxy
        material
        er
        ur
    end
    
    methods
        function Rect = Rectangle(center,rectxy,nlayer,material,varargin)            
            Rect = Rect@PatternShape(center,nlayer);
            if numel(rectxy) == 2
                Rect.rectxy = rectxy;
            else
                error('Check the rectxy input')
            end
            if nargin == 4
                Rect.er = material.er;
                Rect.ur = material.ur;
                Rect.material = material.MaterialName;
            elseif nargin == 5
                Rect.er = varargin{1}(1,1);
                Rect.ur = varargin{1}(1,2);
                Rect.material = material;
            else 
                error('There is a problem with the number of inputs');
            end                         

        end
        
        
        function BuildPattern(Rect,Dev,varargin)
            % Purpose: make rectangle pattern on the device
            % Input: object cylinder device and wavelength 
             
            x0 = Rect.center(1);
            y0 = Rect.center(2);
            rectx = Rect.rectxy(1);
            recty = Rect.rectxy(2);
            if x0-rectx/2 < 0 || y0-recty/2 < 0
                error('The rectangle is too large!')
            end
            Nx = Dev.idimension(1);
            Ny = Dev.idimension(2);
            Lx = Dev.xydimension(1);
            Ly = Dev.xydimension(2);
%             ER = Dev.ER;
%             UR = Dev.UR;            
            dx = Lx/Nx;
            dy = Ly/Ny;
            nx = ceil(rectx/(2*dx));
            nx0 = round(x0*Nx/Lx);
            nx1 = nx0-nx+1;
            nx2 = nx0+nx-1;
            ny = ceil(recty/(2*dy));
            ny0 = round(y0*Ny/Ly);
            ny1 = ny0-ny+1;
            ny2 = ny0+ny-1;
            if nargin == 2
                Dev.ER(nx1:nx2,ny1:ny2,Rect.nlayer) =  Rect.er;
                Dev.UR(nx1:nx2,ny1:ny2,Rect.nlayer) =  Rect.ur;
                
%                 for n = nx1:nx2
%                     for m = ny1:ny2
%                         Dev.ER(n,m,Rect.nlayer) =  Rect.er;
%                         Dev.UR(n,m,Rect.nlayer) =  Rect.ur;
%                     end
%                 end
            elseif nargin == 3
                Dev.ER(nx1:nx2,ny1:ny2,Rect.nlayer) = Rect.er(varargin{1},2);
                Dev.UR(nx1:nx2,ny1:ny2,Rect.nlayer) = Rect.ur(varargin{1},2);
%                 for n = nx1:nx2
%                     for m = ny1:ny2
%                         Dev.ER(n,m,Rect.nlayer) =  Rect.er(varargin{1},2);
%                         Dev.UR(n,m,Rect.nlayer) =  Rect.ur(varargin{1},2);
%                     end
%                 end
            else
                error('Check input number')
            end
                
        end
      
    end
    
end
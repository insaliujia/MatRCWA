classdef Cylinder < PatternShape
    % This class is the subclass of PatternShape, it is used to further
    % define the shape of the pattern
    % Input example: Cyl = Cylinder([1,1],3,[2:4],'Si');
    properties
        name = 'Cylinder'
        radius
        material
        er
        ur
    end
    
    methods
        function Cylin = Cylinder(center,radius,nlayer,material)
            Cylin = Cylin@PatternShape(center,nlayer);
            if numel(radius) == 1
                Cylin.radius = radius;
            else
                error('Check the radius input')
            end
            Cylin.er = material.er;
            Cylin.ur = material.ur;
            Cylin.material = material.MaterialName;
                
        end
        
        function BuildPattern(Cylin,Dev,varargin)
            % Purpose: make cylinder pattern on the device
            % Input: object cylinder device and wavelength 
                
            x0 = Cylin.center(1);
            y0 = Cylin.center(2);

            r = Cylin.radius;
            if x0-r <= 0 || y0-r <= 0
                error('The radius is too large!')
            end
            Nx = Dev.idimension(1);
            Ny = Dev.idimension(2);
            Lx = Dev.xydimension(1);
            Ly = Dev.xydimension(2);
%             ER = Dev.ER;
%             UR = Dev.UR;
            dx = Lx/Nx;
            dy = Ly/Ny;
            nx = round(r/dx);
            nx0 = round(x0*Nx/Lx);
            nx1 = nx0 - nx;
            nx2 = nx0 + nx;

            if nargin == 2
                for n = nx1:nx2
                    ny1 = real(-sqrt(r^2 - (n*dx-x0)^2)+y0);
                    ny2 = real(sqrt(r^2 - (n*dx-x0)^2)+y0);
                    ny1 = round(ny1/dy);
                    ny2 = round(ny2/dy);
                    Dev.ER(n,ny1:ny2,Cylin.nlayer)=Cylin.er;
                    Dev.UR(n,ny1:ny2,Cylin.nlayer)=Cylin.ur;
                end
            elseif nargin == 3
                for n = nx1:nx2
                    ny1 = real(-sqrt(r^2 - (n*dx-x0)^2)+y0);
                    ny2 = real(sqrt(r^2 - (n*dx-x0)^2)+y0);
                    ny1 = round(ny1/dy);
                    ny2 = round(ny2/dy);
                    Dev.ER(n,ny1:ny2,Cylin.nlayer)=Cylin.er(varargin{1},2);
                    Dev.UR(n,ny1:ny2,Cylin.nlayer)=Cylin.ur(varargin{1},2);
                end
            else
                error('Check input number')
            end
                
        end
      
    end
    
end
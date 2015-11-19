classdef Ellipse < PatternShape
    % This class is the subclass of PatternShape, it is used to further
    % define the shape of the pattern
    % Input example: Ellip = Ellipse([0.5,0.5],[0.2,0.3],[2:4],'Si',[1,1]);
    % from left to right is center of the pattern, radius y and radius y, in which
    % layer the pattern will be,the name of the material
    properties
        name = 'Ellipse'
        radius
        material
        er
        ur
    end
    
    methods
        function Ellip = Ellipse(center,radius,nlayer,material)
            Ellip = Ellip@PatternShape(center,nlayer);
            if numel(radius) == 2
                Ellip.radius = radius;
            else
                error('Check the radius input')
            end
            Ellip.er = material.er;
            Ellip.ur = material.ur;
            Ellip.material = material.MaterialName;
                
        end
        
        function BuildPattern(Ellip,Dev,varargin)
            % Purpose: make cylinder pattern on the device
            % Input: object cylinder device and wavelength 
                
            x0 = Ellip.center(1);
            y0 = Ellip.center(2);

            r1 = Ellip.radius(1,1);
            r2 = Ellip.radius(1,2);
            Nx = Dev.idimension(1);
            Ny = Dev.idimension(2);
            Lx = Dev.xydimension(1);
            Ly = Dev.xydimension(2);            
            if x0-r1 <= 0 || Lx-x0-r1 <= 0 || y0-r2 <= 0 || Ly-y0-r2 <=0
                error('The radius is too large!')
            end

%             ER = Dev.ER;
%             UR = Dev.UR;
            dx = Lx/Nx;
            dy = Ly/Ny;
            nx = round(r1/dx);
            ny = round(r2/dy);
            nx0 = round(x0*Nx/Lx);
            nx1 = nx0 - nx;
            nx2 = nx0 + nx;

            if nargin == 2
                for n = nx1:nx2
                    ny1 = real(-sqrt(1-((n*dx-x0)^2)/r1^2)*r2+y0);
                    ny2 = real(sqrt(1-((n*dx-x0)^2)/r1^2)*r2+y0);
                    ny1 = round(ny1/dy);
                    ny2 = round(ny2/dy);
                    Dev.ER(n,ny1:ny2,Ellip.nlayer)=Ellip.er;
                    Dev.UR(n,ny1:ny2,Ellip.nlayer)=Ellip.ur;
                end
            elseif nargin == 3
                for n = nx1:nx2
                    ny1 = real(-sqrt(1-((n*dx-x0)^2)/r1^2)*r2+y0);
                    ny2 = real(sqrt(1-((n*dx-x0)^2)/r1^2)*r2+y0);
                    ny1 = round(ny1/dy);
                    ny2 = round(ny2/dy);
                    Dev.ER(n,ny1:ny2,Ellip.nlayer)=Ellip.er(varargin{1},2);
                    Dev.UR(n,ny1:ny2,Ellip.nlayer)=Ellip.ur(varargin{1},2);
                end
            else
                error('Check input number')
            end
                
        end
      
    end
    
end
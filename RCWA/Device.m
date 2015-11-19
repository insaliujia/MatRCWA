classdef Device < handle
    % device is the class for constucting the whole device
    % functions includ: AddLayer, BuildDevice, CylinderPatternDevice,
    % RectanglePatternDevice,ConvDevice, ShowLayer
    
    properties (SetAccess = protected)
        xydimension             % dimension of the device format [x,y]
        length                  % length of each layer
        material                % material of each layer
        idimension              % imaginary dimension of the layer format [Nx,Ny],
                                % Nx:number of point along x in real-space grid Ny:number 
                                % of point along y in real-space grid
        ilayer                  % number of the layer
        PQR                     % number of spatial harmonics along x and y (HAS TO BE ODD NUM)
        Pattern@cell            % patterns on the substract
        URC                     % caculated convolution matrix of permittivity of the device
        ERC                     % caculated convolution matrix of permeability of the device
    end
    properties
        ER                      % constructed permittivity of the device
        UR                      % constructed permeability of the device
    end

    methods
        function Dev = Device(xydimension,idimension,PQR)
        % Purpose: Define the reflective and tranmission region of the device and period
        % Input: permeability and permittivity in reflection and transmission region, and period
        % The xydimension and idimension should be proportional
            if numel(xydimension) == 2 && numel(idimension) == 2
                Dev.xydimension = xydimension;
                Dev.idimension = idimension;
            else
                error('Check the device dimension')
            end
            Dev.PQR = PQR;
        end
        
        function AddLayer(Dev,material,length,ilayer)
        % Purpose define properties of each layer for the device
            Dev.material{end+1} = material;
            Dev.length(end+1) = length;
            Dev.ilayer(end+1) = ilayer;
        end
        
        function AddPattern(Dev,shape,center,size,nlayer,material)
            % Purpose: add different paterns on specific layers

            switch shape
                case 'Cylinder'
                    radius = size;
                    Dev.Pattern = [Dev.Pattern,{Cylinder(center,radius,nlayer,material)}];
                case 'Ellipse'
                    radius = size;
                    Dev.Pattern = [Dev.Pattern,{Ellipse(center,radius,nlayer,material)}];
                case 'Rectangle'
                    rectxy = size;
                    Dev.Pattern = [Dev.Pattern,{Rectangle(center,rectxy,nlayer,material)}];
                case 'Triangle'
                    SideLength = size;
                    Dev.Pattern = [Dev.Pattern,{Triangle(center,SideLength,nlayer,material)}];
                case 'Pyramid'
                    rectxy = size;
                    Dev.Pattern = [Dev.Pattern,{Pyramid(center,rectxy,nlayer,material)}];
                otherwise
                error('Please Check your spell.');
            end
        end
        
        function BuildLayer(Dev,varargin)
            % Purpose: assign er and ur to different layers
            % Input: object
            % Output: constructed er and ur
            % INITIALIZE LAYERS TO er AND ur
            if nargin == 2
                Dev.ER = ones(Dev.idimension(1),Dev.idimension(2),sum(Dev.ilayer));
                Dev.UR = ones(Dev.idimension(1),Dev.idimension(2),sum(Dev.ilayer));
                low = 1;
                for n = 1:numel(Dev.ilayer)
                    high = sum(Dev.ilayer(1:n));
                    if numel(Dev.material{n}.er) == 1 && numel(Dev.material{n}.ur) == 1
                        Dev.ER(:,:,low:high) = Dev.material{n}.er ; 
                        Dev.UR(:,:,low:high) = Dev.material{n}.ur ; 
                    else
                        Dev.ER(:,:,low:high) = Dev.material{n}.er(varargin{1},2) ;
                        Dev.UR(:,:,low:high) = Dev.material{n}.ur(varargin{1},2) ; 
                    end
                    low = high+1;
                end
            elseif nargin == 1
                Dev.ER = ones(Dev.idimension(1),Dev.idimension(2),sum(Dev.ilayer));
                Dev.UR = ones(Dev.idimension(1),Dev.idimension(2),sum(Dev.ilayer));
                low = 1;
                for n = 1:numel(Dev.ilayer)
                    high = sum(Dev.ilayer(1:n));
                    Dev.ER(:,:,low:high) = Dev.material{n}.er ; 
                    Dev.UR(:,:,low:high) = Dev.material{n}.ur ; 
                    low = high+1;    
                end
            end
        end
        
        function BuildPattern(Dev,varargin)
            % Purpose: build different patterns on the device
            % Usage: The function will take the pattern object out and compute the
            % ER and UR in its object and transfer the data back. The
            % function here will handle the er and ur 
            for n = 1:numel(Dev.Pattern)
                if numel(Dev.Pattern{n}.er) == 1
                    BuildPattern(Dev.Pattern{n},Dev);
                else
                    BuildPattern(Dev.Pattern{n},Dev,varargin{1});
                end
            end     
        end
        
        function BlurDevice(Dev,BlurCoef)
            % Purpose: blur the inhomogeneous layer of the device only ER
            layers = [];
            for n = 1:numel(Dev.Pattern)
                layers = [layers,Dev.Pattern{1,n}.nlayer];    
            end
            for i = layers
                Dev.ER(:,:,i) = Blur(Dev.ER(:,:,i),BlurCoef);
            end

        end
        
        function ShrinkDevice(Dev,ShrinkCoef)
            % Purpose: Shrink each layer of the device the effect including
            % blur and shrink the device
            Blockdimensionx = (Dev.idimension(1) - mod(Dev.idimension(1), ShrinkCoef))/ShrinkCoef;
            Blockdimensiony = (Dev.idimension(2) - mod(Dev.idimension(2), ShrinkCoef))/ShrinkCoef;
            BlockER = ones(Blockdimensionx,Blockdimensiony,sum(Dev.ilayer));
%             BlockUR = ones(Blockdimensionx,Blockdimensiony,sum(Dev.ilayer));
            layers = [];
            for n = 1:numel(Dev.Pattern)
                layers = [layers,Dev.Pattern{1,n}.nlayer];    
            end
            for i = layers
                Dev.ER(:,:,i) = Blur(Dev.ER(:,:,i),ShrinkCoef);
%                 Dev.UR(:,:,i) = Blur(Dev.UR(:,:,i),ShrinkCoef);
                BlockER(:,:,i) = BlockMean(Dev.ER(:,:,i),ShrinkCoef);
%                 BlockUR(:,:,i) = BlockMean(Dev.UR(:,:,i),ShrinkCoef);
            end            
            Dev.ER = BlockER;
%             Dev.UR = BlockUR;           
        end            
            
        
        
      
 
         function ConvDevice(Dev)
             % Purpose: caculate the convolution matrix of the device
             NH  = prod(Dev.PQR);                   %total number of spatial harmonics
             Dev.URC = ones(NH,NH,sum(Dev.ilayer));
             Dev.ERC = ones(NH,NH,sum(Dev.ilayer));
             for n = 1 : sum(Dev.ilayer)
                 Dev.URC(:,:,n) = convmat(Dev.UR(:,:,n),Dev.PQR);            % 1 layer convolution matrices for ur
                 Dev.ERC(:,:,n) = convmat(Dev.ER(:,:,n),Dev.PQR);            % 1 layer convolution matrices for ur 
             end
             
         end
        
         function ShowLayer(Dev,nlayer)
             % Purpose: show certain layer of builded device
             % Input: nlayer--the layer which will be shown
             figure 
             subplot(121);
             imagesc(real(Dev.ER(:,:,nlayer)'));
             title(['ER image of layer: ', num2str(nlayer)])
             axis equal;
             set(gca,'XAxisLocation','top');
             colorbar;
             subplot(122);
             imagesc(real(Dev.UR(:,:,nlayer)'));
%              xlabel('x (\mum)');
%              ylabel('y (\mum)');
             title(['UR image of layer: ', num2str(nlayer)])
             axis equal;
%              set(gca,'YDir','reverse');
%              set(gca,'XAxisLocation','top');
%              colormap(flipud(autumn));
             colorbar;
         end
         
         function ShowConvLayer(Dev,nlayer)
             % Purpose: show certain layer of builded device
             % Input: nlayer--the layer which will be shown
             figure 
             subplot(121);
             imagesc(real(Dev.ERC(:,:,nlayer)'));
             title(['Convolution ER image of layer: ', num2str(nlayer)])
             axis equal;
             colorbar;
             subplot(122);
             imagesc(real(Dev.URC(:,:,nlayer)'));
             xlabel('x (\mum)');
             ylabel('y (\mum)');
             title(['Convolution UR image of layer: ', num2str(nlayer)])
             axis equal;
             colorbar;
         end

         
         function numlayer = numlayer(Dev)
            % Purpose; cacualte the number of layers
            numlayer = numel(Dev.length);
         end        
            
         function LengthLayer = GetLengthLayer(Dev)
             % Purpose: get the length of each layer
             low = 1;
             for n = 1:numel(Dev.length)
                 high = sum(Dev.ilayer(1:n));
                 LengthLayer(low:high) = Dev.length(n)/Dev.ilayer(n);
                 low = high+1; 
             end
         end
         

         
%          function ShowDevice(Dev)
%              Lx = Dev.xydimension(1);
%              Ly = Dev.xydimension(2);
%              Lz = sum(Dev.length);
%              % Vertices
%              VertexData = [Lx*ones(8,1),Ly*ones(8,1),Lz*ones(8,1)]...
%                  .*[0,0,0;
%                  1,0,0;
%                  0,1,0;
%                  0,0,1;
%                  1,1,0;
%                  0,1,1;
%                  1,0,1;
%                  1,1,1];
%              % Patches
%              Index_Patch = ...
%                  [1,2,5,3;
%                  1,3,6,4;
%                  1,4,7,2;
%                  4,7,8,6;
%                  2,5,8,7;
%                  3,6,8,5];
%              n_pat = 6;
%              for i_pat=1:n_pat
%                  % Patches data
%                  PatchData_X(:,i_pat) = VertexData(Index_Patch(i_pat,:),1);
%                  PatchData_Y(:,i_pat) = VertexData(Index_Patch(i_pat,:),2);
%                  PatchData_Z(:,i_pat) = VertexData(Index_Patch(i_pat,:),3);
%              end
%              % Draw patches
%              figure(1);
%              h = patch(PatchData_X,PatchData_Y,PatchData_Z,'y');
%              set(h,'FaceLighting','phong','EdgeLighting','phong');
%              % Axes settings
%              xlabel('x','FontSize',14);
%              ylabel('y','FontSize',14);
%              zlabel('z','FontSize',14);
%              set(gca,'FontSize',14);
%              axis vis3d equal;
%              view([-37.5,30]);
%              camlight;
%              grid on;
%              xlim([-0.15,0.35]);
%              ylim([-0.2,0.3]);
%              zlim([-0.1,0.4]);
%     end
         

         
%         function PyramidePatternDevice(Dev,center,SideLength,nlayer,varargin)
%             if nargin == 5
%                 fillerer = varargin{1}(1,1);
%                 fillerur = varargin{1}(1,2);
%             elseif nargin == 6
%                 fillerer = varargin{1}.er(varargin{2},2);
%                 fillerur = varargin{1}.ur(varargin{2},2);
%             else 
%                 error('There is a problem with the number of inputs');
%             end            
%             layernum = numel(nlayer);
%             gap = SideLength/(layernum-1);
%             p = 0;
%             
%             for i =  nlayer
%                 RectanglePatternDevice(Dev,center,[SideLength-p*gap,SideLength-p*gap],i,[fillerer,fillerur])
%                 p = p+1; 
%             end
%         end
        
%        
           
         
%           function CylinderPatternDevice(Dev,center,radius,nlayer,varargin)
%             % Purpose: make cylinder pattern on the device
%             % Input: center--center of the hole, radius--radius of the
%             % hole,nlayer--in which layer will the hole be, filler--the
%             % material that will fill the hole [er,ur] or name of the
%             % material and its conresponding wavelength
%             % Input format: Dev = Dev.CylinderPatternDevice([0.15,0.15],0.1,[1:3],7);
%             if nargin == 5
%                 fillerer = varargin{1}(1,1);
%                 fillerur = varargin{1}(1,2);
%             elseif nargin == 6
%                 fillerer = varargin{1}.er(varargin{2},2);
%                 fillerur = varargin{1}.ur(varargin{2},2);
%             else 
%                 error('There is a problem with the number of inputs');
%             end
%                 
%             x0 = center(1);
%             y0 = center(2);
% 
%             r = radius;
%             if x0-r <= 0 || y0-r <= 0
%                 error('The radius is too large!')
%             end
%             Nx = Dev.idimension(1);
%             Ny = Dev.idimension(2);
%             Lx = Dev.xydimension(1);
%             Ly = Dev.xydimension(2);
%             dx = Lx/Nx;
%             dy = Ly/Ny;
%             nx = round(r/dx);
%             nx0 = round(x0*Nx/Lx);
%             nx1 = nx0 - nx;
%             nx2 = nx0 + nx;
%             for n = nx1:nx2
%                 ny1 = real(-sqrt(r^2 - (n*dx-x0)^2)+y0);
%                 ny2 = real(sqrt(r^2 - (n*dx-x0)^2)+y0);
%                 ny1 = round(ny1/dy);
%                 ny2 = round(ny2/dy);
%                 Dev.ER(n,ny1:ny2,nlayer)=fillerer;
%                 Dev.UR(n,ny1:ny2,nlayer)=fillerur;
%             end                  
%         end
%         
%         function RectanglePatternDevice(Dev,center,rectxy,nlayer,varargin)
%              % Purpose: make rectangle pattern on the device
%              % Input: center--the center of the rectangle,rectx--width in
%              % the x direction,recty--width in y direction, nlayer--in which layer will the rectangle
%              % be, filler--the material that will fill the hole.
%              % Inputformat: Dev =Dev.RectanglePatternDevice([0.2,0.25],[0.3,0.3],[0.19,0.2],[3:5],[2,1]);
%             if nargin == 5
%                 fillerer = varargin{1}(1,1);
%                 fillerur = varargin{1}(1,2);
%             elseif nargin == 6
%                 fillerer = varargin{1}.er(varargin{2},2);
%                 fillerur = varargin{1}.ur(varargin{2},2);
%             else 
%                 error('There is a problem with the number of inputs');
%             end             
%             x0 = center(1);
%             y0 = center(2);
%             rectx = rectxy(1);
%             recty = rectxy(2);
%             if x0-rectx/2 < 0 || y0-recty/2 < 0
%                 error('The rectangle is too large!')
%             end
%             Nx = Dev.idimension(1);
%             Ny = Dev.idimension(2);
%             Lx = Dev.xydimension(1);
%             Ly = Dev.xydimension(2);
%             dx = Lx/Nx;
%             dy = Ly/Ny;
%             nx = ceil(rectx/(2*dx));
%             nx0 = round(x0*Nx/Lx);
%             nx1 = nx0-nx+1;
%             nx2 = nx0+nx-1;
%             ny = ceil(recty/(2*dy));
%             ny0 = round(y0*Ny/Ly);
%             ny1 = ny0-ny+1;
%             ny2 = ny0+ny-1;            
%             for n = nx1:nx2
%                 for m = ny1:ny2
%                     Dev.ER(n,m,nlayer) = fillerer;
%                     Dev.UR(n,m,nlayer) = fillerur;
%                 end
%             end
%         end
%         
%         function TrianglePatternDevice(Dev,center,SideLength,nlayer,varargin)
%              % Purpose: make equilateral triangle pattern on the device
%              % Input: center--the center of the rectangle,SideLength--side
%              % length,nlayer--in which layer will the rectangle
%              % be, filler--the material that will fill the hole.
%              % Inputformat: Dev =Dev.RectanglePatternDevice([0.2,0.25],0.1,[3:5],[2,1]);
%             if nargin == 5
%                 fillerer = varargin{1}(1,1);
%                 fillerur = varargin{1}(1,2);
%             elseif nargin == 6
%                 fillerer = varargin{1}.er(varargin{2},2);
%                 fillerur = varargin{1}.ur(varargin{2},2);
%             else 
%                 error('There is a problem with the number of inputs');
%             end                           
%             Nx = Dev.idimension(1);
%             Ny = Dev.idimension(2);
%             Lx = Dev.xydimension(1);
%             Ly = Dev.xydimension(2);
%             dx = Lx/Nx;
%             dy = Ly/Ny;
%             h = 0.5*sqrt(3)*SideLength; 
%             x0 = center(1);
%             y0 = center(2);
%             if abs((Lx/2-x0))+SideLength/2 > Lx/2 || h*2/3-(Ly/2-y0)>Ly/2 ...
%                     ... || h/3 + (Ly/2-y0)>Ly/2
%                 error('The rectangle is too large!')
%             end            
%             nxm = floor((Lx/2-x0)/dx);
%             nym = floor((Ly/2-y0)/dy);
%             ny = round(h/dy); 
%             ny1 = round((Ny - ny)/2)+nym; 
%             ny2 = ny1 + ny - 1; 
%             for ny = ny1 : ny2 
%                 f = (ny - ny1)/(ny2 - ny1); 
%                 nx = round(f*SideLength/Lx*Nx); 
%                 nx1 = 1 + floor((Nx - nx)/2)-nxm; 
%                 nx2 = nx1 + nx; 
%                 Dev.ER(nx1:nx2,ny,nlayer) = fillerer;
%                 Dev.UR(nx1:nx2,ny,nlayer) = fillerur;
%             end
%         end

    end
   
end
      
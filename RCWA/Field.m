classdef Field < handle
    % This object will record the data which will be used to caculate
    % the field inside the device
    % Functions will be used to caculate field inside the device
    properties 
        wavenumber
        L
        xydimension
        W_V_ref
        W_V_trn
        W_V
        LAM
        Kx
        Ky
        k0
        URC
        ERC
        C_input
    end
    properties (SetAccess = protected)
        CoordXYZ
        GridXY
        LayerZ
        GridZ
        Ex
        Ey
        Ez
        Hx
        Hy
        Hz
    end
    
    methods
        function GetPointField(Field,CoordinateXYZ)
            % This function get the parameters which will be caculated
            % later 
            if numel(CoordinateXYZ) == 3
                Field.CoordXYZ(1) = CoordinateXYZ(1);
                Field.CoordXYZ(2) = CoordinateXYZ(2);
                Field.CoordXYZ(3) = CoordinateXYZ(3);
            else
                error('Check the number of the Coordinate')
            end
        end
        
        function GetLayerField(Field,XYGird,ZLayer)
            % This function get the layer parameters which will be caculated
            if numel(XYGird) == 2
                Field.GridXY(1) = XYGird(1);
                Field.GridXY(2) = XYGird(2);
            else
                error('Check the number in the X and Y direction');
            end
            Field.LayerZ     = ZLayer;
        end
        
        function GetGridField(Field,XYGird,Z_Gird)
            % This function get the grid parameters which will be caculated
            if numel(XYGird) == 2
                Field.GridXY(1) = XYGird(1);
                Field.GridXY(2) = XYGird(2);
            else
                error('Check the number in the X and Y direction')
            end
            Field.GridZ = Z_Gird;
        end
        
        function CaculateLayerField(Field)
            % This function caculate the field in grid
            [numLayer,Z_position] = FindLayer(Field,Field.LayerZ);
            NH = length(Field.Kx);
            Z = zeros(2*NH,2*NH);
%             C1 = Field.W_V{1}\Field.W_V_ref*Field.C_input;          % define input C at the first layer
            S_U = Field.W_V_ref*Field.C_input;                       % defint the input C in the first layer 
            if numLayer ~= 1
                for n = 1 : numLayer-1
                    Middle = diag(Field.LAM{n})*Field.k0*Field.L(n);
                    Elamda1 = diag(exp(-Middle));       % cacualte diangal lamda
                    Elamda2 = diag(exp( Middle)); 
                    Elamda_Diag = [Elamda1,Z;Z,Elamda2];
                    C = Field.W_V{n}\S_U;          % define input C at the each layer
                    S_U = Field.W_V{n}*Elamda_Diag*C;    % input for the first layer
                end
            end
            Middle =diag(Field.LAM{numLayer})*Field.k0*Z_position;
            Elamda1 = diag(exp(-Middle));
            Elamda2 = diag(exp(Middle));
            Elamda_Diag = [Elamda1,Z;Z,Elamda2];
            C = Field.W_V{numLayer}\S_U;
            S_U = Field.W_V{numLayer}*Elamda_Diag*C;
            S_x = S_U(1:NH,1);
            S_y = S_U(NH+1:2*NH,1);
            U_x = S_U(2*NH+1:3*NH,1);
            U_y = S_U(3*NH+1:4*NH,1);
            S_z = 1i*Field.URC(:,:,numLayer)\(Field.Kx*U_y - Field.Ky*U_x);
            U_z = 1i*Field.ERC(:,:,numLayer)\(Field.Kx*S_y - Field.Ky*S_x);
            
            % cacuate E and H at this layer
            Nx = Field.GridXY(1);
            Ny = Field.GridXY(1);
            Lx = Field.xydimension(1);
            Ly = Field.xydimension(2);
            % initilize grid
            dx = Lx/Nx;                        %grid resolution along x
            dy = Ly/Ny;                        %grid resolution along y
            xa = [0:Nx-1]*dx;                      %x axis array
            xa = xa - mean(xa);                    %center x axis at zero
            ya = [0:Ny-1]*dy;                      %y axis vector
            ya = ya - mean(ya);                    %center y axis at zero
            % initilize E and H
            IniZ = zeros(Nx,Ny);
            Field.Ex = IniZ;
            Field.Ey = IniZ;
            Field.Ez = IniZ;
            Field.Hx = IniZ;
            Field.Hy = IniZ;
            Field.Hz = IniZ;
            
            for m = 1:numel(ya)
                for n = 1:numel(xa)
                    Phase = exp(-1i*(diag(Field.Kx)*xa(n)+diag(Field.Ky)*ya(m)));
                    Field.Ex(n,m) = sum(S_x.*Phase);
                    Field.Ey(n,m) = sum(S_y.*Phase);
                    Field.Ez(n,m) = sum(S_z.*Phase);
                    Field.Hx(n,m) = sum(U_x.*Phase);
                    Field.Hy(n,m) = sum(U_y.*Phase);
                    Field.Hz(n,m) = sum(U_z.*Phase);
                end
            end 
        end
        
        function CaculateGridField(Field)
            % This function caculate field in a grid
            % since here we take count of the first reflective layer the
            % device will start from the second layer 
            % initilize num_position_Layer
            num_layer = numel(Field.GridZ);
            num_position_Layer = zeros(num_layer,2);
            for i = 1:num_layer
                [num_position_Layer(i,1),num_position_Layer(i,2)] = FindLayer(Field,Field.GridZ(i));
            end
            NH = length(Field.Kx);
            Z = zeros(2*NH,2*NH);
            S_U_C = Field.W_V_ref*Field.C_input;                       % defint the input C in the first layer
            % initilize S_U
            C_input_layer = cell(1,max(num_position_Layer(:,1)));
            % caculate the input C for each layer
            C_input_layer{1} = S_U_C;
            % caculate input in every layer in which will caculate the field
            for n = 1 : max(num_position_Layer(:,1))-1 
                Middle = diag(Field.LAM{n})*Field.k0*Field.L(n);
                Elamda1 = diag(exp(-Middle));       % cacualte diangal lamda
                Elamda2 = diag(exp( Middle));
                Elamda_Diag = [Elamda1,Z;Z,Elamda2];
                C = Field.W_V{n}\C_input_layer{n};          % define input C at the each layer
                C_input_layer{n+1} = Field.W_V{n}*Elamda_Diag*C;    % input for the first layer
            end
            
            %initilize each layer
            S_U = cell(1,num_layer);
            S_x = cell(1,num_layer);
            S_y = cell(1,num_layer);
            S_z = cell(1,num_layer);
            U_x = cell(1,num_layer);
            U_y = cell(1,num_layer);
            U_z = cell(1,num_layer);
            % caculate field in each layer
            for m = 1 : num_layer
                caculated_layer = num_position_Layer(m,1);
                Middle = diag(Field.LAM{caculated_layer})*Field.k0*num_position_Layer(m,2);
                Elamda1 = diag(exp(-Middle));       % cacualte diangal lamda
                Elamda2 = diag(exp( Middle));
                Elamda_Diag = [Elamda1,Z;Z,Elamda2];
                C = Field.W_V{caculated_layer}\C_input_layer{caculated_layer};          % define input C at the each layer
                S_U{m} = Field.W_V{caculated_layer}*Elamda_Diag*C;    % input for the first layer
                S_x{m} = S_U{m}(1:NH,1);
                S_y{m} = S_U{m}(NH+1:2*NH,1);
                U_x{m} = S_U{m}(2*NH+1:3*NH,1);
                U_y{m} = S_U{m}(3*NH+1:4*NH,1);
                S_z{m} = 1i*Field.URC(:,:,caculated_layer)\(Field.Kx*U_y{m} - Field.Ky*U_x{m});
                U_z{m} = 1i*Field.ERC(:,:,caculated_layer)\(Field.Kx*S_y{m} - Field.Ky*S_x{m});
            end
            clear S_U  C_input_layer;
            
            % cacuate E and H at this layer
            Nx = Field.GridXY(1);
            Ny = Field.GridXY(1);
            Lx = Field.xydimension(1);
            Ly = Field.xydimension(2);
            % initilize grid
            xa = linspace(-Lx/2,Lx/2,Nx+1);
            ya = linspace(-Ly/2,Ly/2,Ny+1);
%             [xxa,yya] = meshgrid(xa,ya);
%             dx = Lx/Nx;                        %grid resolution along x
%             dy = Ly/Ny;                        %grid resolution along y
%             xa = [0:Nx-1]*dx;                      %x axis array
%             xa = xa - mean(xa);                    %center x axis at zero
%             ya = [0:Ny-1]*dy;                      %y axis vector
%             ya = ya - mean(ya);                    %center y axis at zero
            % initilize E and H
            IniZ = zeros(Nx,Ny,num_layer);
            Ex_ = IniZ;
            Ey_ = IniZ;
            Ez_ = IniZ;
            Hx_ = IniZ;
            Hy_ = IniZ;
            Hz_ = IniZ;        
            Kx_ = diag(Field.Kx);
            Ky_ = diag(Field.Ky);
            % cauclate E and H field         
            for i = 1:num_layer
                for m = 1:numel(ya)
                    for n = 1:numel(xa)
                        Phase = exp(-1i*(Kx_*xa(n)+Ky_*ya(m)));
                        Ex_(n,m,i) = sum(S_x{i}.*Phase);
                        Ey_(n,m,i) = sum(S_y{i}.*Phase);
                        Ez_(n,m,i) = sum(S_z{i}.*Phase);
                        Hx_(n,m,i) = sum(U_x{i}.*Phase);
                        Hy_(n,m,i) = sum(U_y{i}.*Phase);
                        Hz_(n,m,i) = sum(U_z{i}.*Phase);
                    end
                end
            end
            Field.Ex = Ex_;
            Field.Ey = Ey_;
            Field.Ez = Ez_;
            Field.Hx = Hx_;
            Field.Hy = Hy_;
            Field.Hz = Hz_;
            clear Ex_ Ey_ Ez_ Hx_ Hy_ Hz_ IniZ
        end
            

        
        function [numLayer,Z_position] = FindLayer(Field,depth)
            % This function will find in which layer the number is
            if depth < 0
                error('Depth should be a positive number')
            elseif depth > sum(Field.L)
                error('Depth should not pass the depth of the device')
            end
            z_depth = depth;
            z_depth_device = Field.L(1);
            numLayer = 1;
            Z_position = z_depth;
            while z_depth > z_depth_device
                numLayer = numLayer + 1;
                Z_position = z_depth - z_depth_device;
                z_depth_device = z_depth_device + Field.L(numLayer);
            end
            
        end
        
        function ShowLayerField(Field,varargin)
            % Purpose: this function show field in certain layer
            if nargin == 1
                figure
                subplot(231)
                imagesc(real(Field.Ex));
                title(['Ex at z layer ', num2str(Field.LayerZ)]);
                subplot(232)
                imagesc(real(Field.Ey));
                title(['Ey at z layer', num2str(Field.LayerZ)]);
                subplot(233)
                imagesc(real(Field.Ez));
                title(['Ez at z layer ', num2str(Field.LayerZ)]);
                subplot(234)
                imagesc(real(Field.Hx));
                title(['Hx at z layer ', num2str(Field.LayerZ)]);
                subplot(235)
                imagesc(real(Field.Hy));
                title(['Hy at z layer ', num2str(Field.LayerZ)]);
                subplot(236)
                imagesc(real(Field.Hz));
                title(['Hz at z layer ', num2str(Field.LayerZ)]);
            elseif nargin == 2
                figure
                switch varargin{1}
                    case 'Ex' 
                        imagesc(real(Field.Ex));
                        title(['Ex at z layer ', num2str(Field.LayerZ)]);
                    case 'Ey'
                        imagesc(real(Field.Ey));
                        title(['Ey at z layer ', num2str(Field.LayerZ)]);
                    case 'Ez'
                        imagesc(real(Field.Ez));
                        title(['Ez at z layer ', num2str(Field.LayerZ)]);
                    case 'Hx'
                        imagesc(real(Field.Hx));
                        title(['Hx at z layer ', num2str(Field.LayerZ)]);
                    case 'Hy'
                        imagesc(real(Field.Hy));
                        title(['Hy at z layer ', num2str(Field.LayerZ)]);
                    case 'Hz'
                        imagesc(real(Field.Hz));
                        title(['Hz at z layer ', num2str(Field.LayerZ)]);
                    case 'E'
                        subplot(131)
                        imagesc(real(Field.Ex));
                        title(['Ex at z layer ', num2str(Field.LayerZ)]);
                        axis square;
                        subplot(132)
                        imagesc(real(Field.Ey));
                        title(['Ey at z layer ', num2str(Field.LayerZ)]);
                        axis square;
                        subplot(133)
                        imagesc(real(Field.Ez));
                        title(['Ez at z layer ', num2str(Field.LayerZ)]);
                        axis square;
                    case 'H'
                        subplot(131)
                        imagesc(real(Field.Hx));
                        title(['Hx at z layer ', num2str(Field.LayerZ)]);
                        axis square;
                        subplot(132)
                        imagesc(real(Field.Hy));
                        title(['Hy at z layer ', num2str(Field.LayerZ)]);
                        axis square;
                        subplot(133)
                        imagesc(real(Field.Hz));
                        title(['Hz at z layer ', num2str(Field.LayerZ)]);
                        axis square;
                end
            else
                error('Too many input')
            end
    
        end
        
        function ShowGridField(Field,varargin)
            if nargin == 1
                figure
                subplot(231)
                h = slice(real(Field.Ex),[],[],1:size(Field.Ex,3));
                set(h, 'EdgeColor','none', 'FaceColor','interp')
                alpha(0.5)
                title(['Ex at grid from ', num2str(Field.GridZ(1)), ' to ',num2str(Field.GridZ(end))]);
                subplot(232)
                h = slice(real(Field.Ey),[],[],1:size(Field.Ey,3));
                set(h, 'EdgeColor','none', 'FaceColor','interp')
                alpha(0.5)
                title(['Ey at grid from ', num2str(Field.GridZ(1)), ' to ',num2str(Field.GridZ(end))]);
                subplot(233)
                h = slice(real(Field.Ez),[],[],1:size(Field.Ez,3));
                set(h, 'EdgeColor','none', 'FaceColor','interp')
                alpha(0.5)
                title(['Ez at grid from ', num2str(Field.GridZ(1)), ' to ',num2str(Field.GridZ(end))]);
                subplot(234)
                h = slice(real(Field.Hx),[],[],1:size(Field.Hx,3));
                set(h, 'EdgeColor','none', 'FaceColor','interp')
                alpha(0.5)
                title(['Hx at grid from ', num2str(Field.GridZ(1)), ' to ',num2str(Field.GridZ(end))]);
                subplot(235)
                h = slice(real(Field.Hy),[],[],1:size(Field.Hy,3));
                set(h, 'EdgeColor','none', 'FaceColor','interp')
                alpha(0.5)
                title(['Hy at grid from ', num2str(Field.GridZ(1)), ' to ',num2str(Field.GridZ(end))]);
                subplot(236)
                h = slice(real(Field.Hz),[],[],1:size(Field.Hz,3));
                set(h, 'EdgeColor','none', 'FaceColor','interp')
                alpha(0.5)
                title(['Hz at grid from ', num2str(Field.GridZ(1)), ' to ',num2str(Field.GridZ(end))]);
            elseif nargin == 2
                figure
                switch varargin{1}
                    case 'Ex'
                        h = slice(real(Field.Ex),[],[],1:size(Field.Ex,3));
                        set(h, 'EdgeColor','none', 'FaceColor','interp')
                        alpha(0.5)
                        title(['Ex at grid from ', num2str(Field.GridZ(1)), ' to ',num2str(Field.GridZ(end))]);
                    case 'Ey'
                        h = slice(real(Field.Ey),[],[],1:size(Field.Ey,3));
                        set(h, 'EdgeColor','none', 'FaceColor','interp')
                        alpha(0.5)
                        title(['Ey at grid from ', num2str(Field.GridZ(1)), ' to ',num2str(Field.GridZ(end))]);
                    case 'Ez'
                        h = slice(real(Field.Ez),[],[],1:size(Field.Ez,3));
                        set(h, 'EdgeColor','none', 'FaceColor','interp')
                        alpha(0.5)
                        title(['Ez at grid from ', num2str(Field.GridZ(1)), ' to ',num2str(Field.GridZ(end))]);
                    case 'Hx'
                        h = slice(real(Field.Hx),[],[],1:size(Field.Hx,3));
                        set(h, 'EdgeColor','none', 'FaceColor','interp')
                        alpha(0.5)
                        title(['Hx at grid from ', num2str(Field.GridZ(1)), ' to ',num2str(Field.GridZ(end))]);
                    case 'Hy'
                        h = slice(real(Field.Hy),[],[],1:size(Field.Hy,3));
                        set(h, 'EdgeColor','none', 'FaceColor','interp')
                        alpha(0.5)
                        title(['Hy at grid from ', num2str(Field.GridZ(1)), ' to ',num2str(Field.GridZ(end))]);
                    case 'Hz'
                        h = slice(real(Field.Hz),[],[],1:size(Field.Hz,3));
                        set(h, 'EdgeColor','none', 'FaceColor','interp')
                        alpha(0.5)
                        title(['Hz at grid from ', num2str(Field.GridZ(1)), ' to ',num2str(Field.GridZ(end))]);
                    case 'E'
                        subplot(131)
                        h = slice(real(Field.Ex),[],[],1:size(Field.Ex,3));
                        set(h, 'EdgeColor','none', 'FaceColor','interp')
                        alpha(0.5)
                        title(['Ex at grid from ', num2str(Field.GridZ(1)), ' to ',num2str(Field.GridZ(end))]);
                        axis square;
                        subplot(132)
                        h = slice(real(Field.Ey),[],[],1:size(Field.Ey,3));
                        set(h, 'EdgeColor','none', 'FaceColor','interp')
                        alpha(0.5)
                        title(['Ey at grid from ', num2str(Field.GridZ(1)), ' to ',num2str(Field.GridZ(end))]);
                        axis square;
                        subplot(133)
                        h = slice(real(Field.Ez),[],[],1:size(Field.Ez,3));
                        set(h, 'EdgeColor','none', 'FaceColor','interp')
                        alpha(0.5)
                        title(['Ez at grid from ', num2str(Field.GridZ(1)), ' to ',num2str(Field.GridZ(end))]);
                        axis square;
                    case 'H'
                        subplot(131)
                        h = slice(real(Field.Hx),[],[],1:size(Field.Hx,3));
                        set(h, 'EdgeColor','none', 'FaceColor','interp')
                        alpha(0.5)
                        title(['Hx at grid from ', num2str(Field.GridZ(1)), ' to ',num2str(Field.GridZ(end))]);
                        axis square;
                        subplot(132)
                        h = slice(real(Field.Hy),[],[],1:size(Field.Hy,3));
                        set(h, 'EdgeColor','none', 'FaceColor','interp')
                        alpha(0.5)
                        title(['Hy at grid from ', num2str(Field.GridZ(1)), ' to ',num2str(Field.GridZ(end))]);
                        axis square;
                        subplot(133)
                        h = slice(real(Field.Hz),[],[],1:size(Field.Hz,3));
                        set(h, 'EdgeColor','none', 'FaceColor','interp')
                        alpha(0.5)
                        title(['Hz at grid from ', num2str(Field.GridZ(1)), ' to ',num2str(Field.GridZ(end))]);
                        axis square;
                end
            else
                error('Too many input')
            end     
        end
        
    end
end









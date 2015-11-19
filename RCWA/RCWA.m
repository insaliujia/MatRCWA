classdef RCWA < handle
    properties(SetAccess=protected)
        ShowProcess             % show the process of the simulation
        BlurCoef = 0          % blur device to accelerate the caculation
        ShrinkCoef = 0         % blur and shrink device to accelerate the caculation
        dispersion = 1          % define dispersion
        RecordField = 0         % reconstructe field inside device 
        RecordDifOrder = 0      % record diffrection order for every wavelength
        referur                 % permeability and permittivity in reflection region
        trnerur                 % permeability and permittivity in transmission region
        R
        T
        Ref_order               % reflection in each order
        Trn_order               % transmission in each order
        source                  % source of the simulation
        device                  % device of the simulation
    end
    
    properties 
        field                   % field of the simualtion
        
    end
    
    properties (Constant,Hidden)
        % Define the units
        micrometers = 1;
        nanometers  = 1 / 1000;
        centimeter = 10000;
        meter = 1000000;
        radians     = 1;
        degrees = pi/180;
        sur_nor = [0; 0; -1];                     % define surface normal
    end
    
    methods
        function ObjRCWA = RCWA(referur,trnerur,ShowProcess)
            % Purpose: Define object RCWA
            % Input: the number of spatial harmonics along x and y (HAS TO BE ODD NUM)
            ObjRCWA.referur = referur;
            ObjRCWA.trnerur = trnerur;
            ObjRCWA.ShowProcess = ShowProcess;
        end
        
        function ObjRCWA = UseBlurEffect(ObjRCWA,BlurCoef)
            % Purpose: open the blur effect to accelerate caculation
            % Input: blur coefficient--How many times the device will be shrinked
            if numel(BlurCoef) == 1
                ObjRCWA.BlurCoef = BlurCoef;
            else
                error('BlurCoef should be a single number')
            end
        end
        
        function ObjRCWA = UseShrinkEffect(ObjRCWA,ShrinkCoef)
            % Purpose: open the blur effect to accelerate caculation
            % Input: blur coefficient--How many times the device will be shrinked
            if numel(ShrinkCoef) == 1
                ObjRCWA.ShrinkCoef = ShrinkCoef;
            else
                error('ShrinkCoef should be a single number')
            end            
            
        end
        
        function ObjRCWA = Dispersion(ObjRCWA,dispersion)
            % Purpose: define the dispersion of materials
            ObjRCWA.dispersion = dispersion;
        end
            
        
        function ObjRCWA = ConstructField(ObjRCWA)
            % Purpose: start record field inside device
            ObjRCWA.RecordField = 1;
%             % initialize field object
%             ObjRCWA.field = 
        end
        
        function ObjRCWA = RecordDiffOrder(ObjRCWA)
           % Purpose: record diffrection order
           ObjRCWA.RecordDifOrder = 1;
        end
        
        function RCWARun(ObjRCWA,source,device,field)
            if ObjRCWA.ShowProcess == 1
                h = waitbar(0,'1','Name','RCWA Caculating...',...
                    'CreateCancelBtn',...
                    'setappdata(gcbf,''canceling'',1)');
                setappdata(h,'canceling',0)
            end
            
            ObjRCWA.source = source;
            ObjRCWA.device = device;
            if  ObjRCWA.dispersion == 0             % deal with dispersion of the mateiral
                BuildLayer(ObjRCWA.device);
                BuildPattern(ObjRCWA.device);
                if ObjRCWA.BlurCoef ~= 0         %open the blur effect
                    BlurDevice(ObjRCWA.device,ObjRCWA.BlurCoef)
                end
                if ObjRCWA.ShrinkCoef ~= 0       %open the shrink effect
                    ShrinkDevice(ObjRCWA.device,ObjRCWA.ShrinkCoef)
                end
                ConvDevice(ObjRCWA.device);
            end
            if ObjRCWA.RecordField == 1             % initialize Filed object to be prepared to record field
                LayerNum = sum(ObjRCWA.device.ilayer);
                ObjRCWA.field = field;
                ObjRCWA.field.W_V = cell(1,LayerNum);
                ObjRCWA.field.LAM = cell(1,LayerNum);
            end
            
            if ObjRCWA.RecordDifOrder == 1             % initialize matrix to record diffrection for each order
%                 x_k = size(device.PQR);
                Ref_order_ = zeros(device.PQR(1),device.PQR(2),source.snum);            % record refection in each order
                Trn_order_ = Ref_order_;            % record transmission in each order
            end
%             if ObjRCWA.RecordField == 1             % initialize Filed object to be prepared to record field
%                 LayerNum = sum(ObjRCWA.device.ilayer);
%                 ObjRCWA.field = Field;
%                 for wavenum = 1: NLAM
%                     ObjRCWA.field(wavenum).wavenumber = wavenum;
%                     ObjRCWA.field(wavenum).W_V = cell(1,LayerNum);
%                     ObjRCWA.field(wavenum).LAM = cell(1,LayerNum);
%                 end
%             end
            NLAM = source.snum;                  %determine how many simulations
            for nlam = 1 : NLAM
                if ObjRCWA.dispersion == 1
                    BuildLayer(ObjRCWA.device,nlam);
                    BuildPattern(ObjRCWA.device,nlam);
                    if ObjRCWA.BlurCoef ~= 0         %open the blur effect
                        BlurDevice(ObjRCWA.device,ObjRCWA.BlurCoef)
                    end
                    if ObjRCWA.ShrinkCoef ~= 0       %open the shrink effect
                        ShrinkDevice(ObjRCWA.device,ObjRCWA.ShrinkCoef)
                    end
                    ConvDevice(ObjRCWA.device);
                end
                % Make patterns in the layer input:[center],radius,[in which ilayer],[er, ur]
                [Ref,Trn] = RCWAer(ObjRCWA,source,ObjRCWA.device,nlam);
                if ObjRCWA.RecordDifOrder == 1
                    Ref_order_(:,:,nlam) = reshape(Ref,device.PQR);            % record refection in each order
                    Trn_order_(:,:,nlam) = reshape(Trn,device.PQR);            % record transmission in each order
                end
                Ref = sum(Ref(:));
                Trn = sum(Trn(:));
                ObjRCWA.R(nlam) = 100*Ref;
                ObjRCWA.T(nlam) = 100*Trn;
                % caculate field in device 
                if ObjRCWA.RecordField == 1
                    if ObjRCWA.field.CoordXYZ ~= 0
                        ObjRCWA.field.CaculatePointField
                    elseif ObjRCWA.field.LayerZ ~= 0
                        ObjRCWA.field.CaculateLayerField
                    elseif sum(ObjRCWA.field.GridZ) ~= 0
                        ObjRCWA.field.CaculateGridField;
                    end
                end
                if ObjRCWA.ShowProcess == 1
                    waitbar(nlam/NLAM,h,'I am working, please don''t disturb me ...');
                    if getappdata(h,'canceling')
                        delete(h)
                        break
                    end
                end
            end
            if ObjRCWA.RecordDifOrder == 1 
                ObjRCWA.Ref_order = Ref_order_;
                ObjRCWA.Trn_order = Trn_order_;
                clear Ref_order_ Trn_order_
            end
            if  ObjRCWA.ShowProcess == 1
                delete(h)
            end
            
        end
        
        
        function PlotRT(ObjRCWA)
            % Create figure
            figure1 = figure('Name','Reflection and Transmission','NumberTitle','off');
            
            % Create axes
            axes1 = axes('Parent',figure1,'FontWeight','demi','FontSize',14);
            box(axes1,'on');
            hold(axes1,'all');
            
            plot(ObjRCWA.source.wavelength/ObjRCWA.nanometers,ObjRCWA.R,'-r','LineWidth',2); hold on;
            plot(ObjRCWA.source.wavelength/ObjRCWA.nanometers,ObjRCWA.T,'-b','LineWidth',2);
            plot(ObjRCWA.source.wavelength/ObjRCWA.nanometers,100-(ObjRCWA.R+ObjRCWA.T),'-k','LineWidth',2); hold off;
            legend('Reflectance', 'Transmittance', 'Conservation')
            axis([min(ObjRCWA.source.wavelength)/ObjRCWA.nanometers max(ObjRCWA.source.wavelength)/ObjRCWA.nanometers 0 105]);
            xlabel('Wavelength (nm)','FontWeight','demi','FontSize',12);
            ylabel('%   ','Rotation',0,'FontWeight','demi','FontSize',12);
            title('SPECTRAL RESPONSE','FontWeight','bold','FontSize',14);
        end
        
        function PlotR(ObjRCWA)
            figure1 = figure('Name','Reflection','NumberTitle','off');
            % Create axes
            axes1 = axes('Parent',figure1,'FontWeight','demi','FontSize',14);
            box(axes1,'on');
            hold(axes1,'all');
            
            plot(ObjRCWA.source.wavelength/ObjRCWA.nanometers,ObjRCWA.R,'-r','LineWidth',2);
            legend('Reflectance')
            axis([min(ObjRCWA.source.wavelength)/ObjRCWA.nanometers max(ObjRCWA.source.wavelength)/ObjRCWA.nanometers 0 105]);
            xlabel('Wavelength (nm)','FontWeight','demi','FontSize',12);
            ylabel('%   ','Rotation',0,'FontWeight','demi','FontSize',12);
            title('Reflection','FontWeight','bold','FontSize',14);
        end
           
        function PlotT(ObjRCWA)
            figure1 = figure('Name','Transmission','NumberTitle','off');
            % Create axes
            axes1 = axes('Parent',figure1,'FontWeight','demi','FontSize',14);
            box(axes1,'on');
            hold(axes1,'all');
            
            plot(ObjRCWA.source.wavelength/ObjRCWA.nanometers,ObjRCWA.T,'-b','LineWidth',2);
            legend('Transmittance')
            axis([min(ObjRCWA.source.wavelength)/ObjRCWA.nanometers max(ObjRCWA.source.wavelength)/ObjRCWA.nanometers 0 105]);
            xlabel('Wavelength (nm)','FontWeight','demi','FontSize',12);
            ylabel('%   ','Rotation',0,'FontWeight','demi','FontSize',12);
            title('Transmission','FontWeight','bold','FontSize',14);
        end
        
        function PlotA(ObjRCWA)
            figure1 = figure('Name','Absorption','NumberTitle','off');
            % Create axes
            axes1 = axes('Parent',figure1,'FontWeight','demi','FontSize',14);
            box(axes1,'on');
            hold(axes1,'all');
            plot(ObjRCWA.source.wavelength/ObjRCWA.nanometers,100-(ObjRCWA.R+ObjRCWA.T),'-k','LineWidth',2);
            legend('Absorption')
            axis([min(ObjRCWA.source.wavelength)/ObjRCWA.nanometers max(ObjRCWA.source.wavelength)/ObjRCWA.nanometers 0 105]);
            xlabel('Wavelength (nm)','FontWeight','demi','FontSize',12);
            ylabel('%   ','Rotation',0,'FontWeight','demi','FontSize',12);
            title('Absorption','FontWeight','bold','FontSize',14);
        end
        
        function SaveData(ObjRCWA,name)
            filename = strcat(name,'.mat');
            save(filename);
        end
                       
    end
            
end
    
        
classdef Material < handle
    % Manipulate er and ur for the material of caculation
    % Function inculde Plotn--plot the reflective index, Selectn--select
    % the reflective index which will be used.
    properties
        MaterialName@char
        n
        er
        ur
        

    end
    
    methods 
        function Mat = Material(material,varargin)
            Mat.MaterialName = material;
            if nargin == 1
                filename = strcat(Mat.MaterialName,'.txt');
                Mat.n = dlmread(filename);
                re_n = Mat.n(:,2);
                im_n = Mat.n(:,3);
                Mat.n(:,2) = conj(complex(re_n,im_n));
                Mat.n(:,3) = [];
                Mat.er = horzcat(Mat.n(:,1),Mat.n(:,2).^2);   %permittivity of device
                Mat.ur = horzcat(Mat.n(:,1),ones(numel(Mat.n)/2,1)); 
            
            
            elseif nargin == 2
                if numel(varargin{1}) == 1
                    Mat.n = varargin{1};
                    Mat.er = Mat.n^2;
                    Mat.ur = 1;
                elseif numel(varargin{1}) == 2
                Mat.er = varargin{1}(1,1);
                Mat.ur = varargin{1}(1,2);
                Mat.n = sqrt(Mat.er*Mat.ur);
                end
                
            else
                error('Too many input!')
            end

        end
        
                
        
        function Plotn(Mat)
            % Purpose: plot reflective index
            figure('name',Mat.MaterialName);
            real_n = real(Mat.n(:,2));
            img_n  = -imag(Mat.n(:,2));
            x_label = Mat.n(:,1);
            plot(x_label*1000,real_n,'-b'); 
            hold on;
            plot(x_label*1000,img_n,'-r'); 
            hold off;
            xlabel('Wavelength (nm)');
            ylabel('Reflectivity');
            title('Real and imaginary refractive index');
            legend('real part','imaginary part');
        end
        
        function Fitn(Mat,index,varargin)
            % Purpose: fit and caculate reflective index of materials
            % input: Mat--material, source--source of simulation,
            % index--fitting index
            % if the fitting is not coordinate with your data, please use
            % more advanced fitting function
            
                p_re = polyfit(Mat.n(:,1),real(Mat.n(:,2)),index);
                p_im = polyfit(Mat.n(:,1),imag(Mat.n(:,2)),index);
                x_wavelength = round(Mat.n(1,1)*1000)/1000:0.001:round(Mat.n(end,1)*1000)/1000;
                f_re = polyval(p_re,x_wavelength);
                f_im = polyval(p_im,x_wavelength);
                if nargin == 2
                    figure('name',Mat.MaterialName);
                    subplot(121);
                    plot(Mat.n(:,1),real(Mat.n(:,2)),'o');
                    hold on
                    plot(x_wavelength,f_re,'r');
                    hold off;
                    xlabel('Wavelength (um)');
                    ylabel('index');
                    title('Real refractive index');
                    subplot(122);
                    plot(Mat.n(:,1),imag(Mat.n(:,2)),'o');
                    hold on
                    plot(x_wavelength,f_im,'r');
                    hold off;
                    xlabel('Wavelength (um)');
                    ylabel('index');
                    title('imaginary refractive index');
                elseif nargin == 3
                    Mat.n(1:numel(x_wavelength),1) = x_wavelength;
                    Mat.n(:,2) = conj(complex(f_re,f_im));
                    Mat.er = horzcat(Mat.n(:,1),Mat.n(:,2).^2);   %permittivity of device
                    Mat.ur = horzcat(Mat.n(:,1),ones(numel(Mat.n)/2,1)); 
                end
                
        end
        
        
        
        
        function Selectn(Mat,source)
            % Purpose: choose reflective index for further usage
            if numel(source.wavelength) > 1
                p = (Mat.n(2,1)-Mat.n(1,1));
                q = source.wavelength(2)-source.wavelength(1);
                if p > q 
                    error('The wavelength of source don''t match the step of reflectivity')
                end
            end
            [row_l,col_l] = find(Mat.n(:,1)== smin(source));
            [row_h,col_h] = find(Mat.n(:,1)== smax(source));
            step = round(q*1000);
            Mat.n = Mat.n(row_l:step:row_h,:);  %reflective index of silicon
            Mat.er = Mat.er(row_l:step:row_h,:);   %permittivity of device
            Mat.ur = Mat.ur(row_l:step:row_h,:);
        end
        
        function num_n = numn(Mat)
            % Purpose: get the number of the reflectivity
            num_n = numel(Mat.n);
        end
         
        
    end

    
    
end


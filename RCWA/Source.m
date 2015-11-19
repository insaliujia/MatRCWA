classdef Source < handle
    % input format S = Source([300:1:400],[30,60],[1,0])
    % source is the the class for incident source
    % Properties for the class
        
    properties (SetAccess = immutable)
        wavelength          % wavelength of the source, unit is nanometer
        angle               % incident angle, unite is radians
        polarization        % polarization of the incident, no unite
    end
    methods
        % Get the input from users
        function s = Source(wavelength,angle,polarization)
        % Purpose: constructor for the source
        % Input: wavelength -- the wavelength range of incident
        %        angle      -- the angle of incident, constructed as array
        %        polarization -- polarization of incident, constructed as
        %        array
                 s.wavelength = wavelength *  RCWA.nanometers;
                 s.angle      = angle * RCWA.degrees;
                 s.polarization = polarization;
        end
        
       
        function disp(s)
        % Usage: disp(s)
        % Purpose: Disply the data 
        
            fprintf('Wavelength is from %s micrometer to %s micrometer\n',num2str(s.wavelength(1)),num2str(s.wavelength(end)));
            fprintf('Incident angle theta is %s radians and phi is %s radians\n', num2str(s.angle(1)), num2str(s.angle(end)));
            fprintf('TE and TM is %s and %s\n', num2str(s.polarization(1)),num2str(s.polarization(2)));
        end
        
        
        function snum = snum(s)
        % Purpose: Get the num of the wavelength
            snum =  length(s.wavelength);
        end
        
        function smin = smin(s)
        % Purpose: Get the minimum wavelength 
            smin = min(s.wavelength);
        end
        
        function smax = smax(s)
        % Purpose: Get the maximum wavelength
            smax = max(s.wavelength);
        end
        
        function sk = sk(s,wavelengthnum)
        % Purpose: Get the wave vector conresponding to the wavelength
        % INPUT: the number of the wavelength
        % OUTPUT: k vector in this wavelength
            sk = 2*pi/s.wavelength(wavelengthnum);
        end
        
        
        function skinc = skinc(s,n)
        % Purpose: get the source vector in the reflective region
        % INPUT: n is the reflective of the place that the incident is in
        % OUTPUT: Source vector k
            skinc = n*[sin(s.angle(1))*cos(s.angle(2));sin(s.angle(1))*sin(s.angle(2));cos(s.angle(1))];
                
        end  
        
        function sP = sP(s,wavelengthnum,n)
            % caculate vector along polarizations
            % Input: the number of the wavelength
            if s.angle(1) == 0
                a_te = [0;1;0];
            else
                a_te = cross((s.sk(wavelengthnum)*s.skinc(n)),RCWA.sur_nor)/norm(cross((s.sk(wavelengthnum)*s.skinc(n)),RCWA.sur_nor));
            end
            a_tm = cross(a_te,(s.sk(wavelengthnum)*s.skinc(n)))/norm(cross(a_te,(s.sk(wavelengthnum)*s.skinc(n))));

            % Composite polarization vector 
            sP = s.polarization(1)*a_te + s.polarization(2)*a_tm;
            sP = sP/norm(sP);
            
        end
            
     end
        
end



        
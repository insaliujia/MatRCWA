classdef PatternShape 
    % Shape is a class to define pattern parameters of the structures on the
    % device
    properties
        center
        nlayer
    end
    
    methods
        
        function Shape = PatternShape(center,nlayer)
            if numel(center) == 2
                Shape.center = center;
            end
            Shape.nlayer = nlayer;
        end  
    end
end
classdef Material
    properties
        E = 30e9;
        v = 0.3;
        C0 = 10;
        CU = 0;
        ro = 2500;
        fsmax = 3e6;
        fsmin = -4e7;
        ftmax = 8e8;
        sigmaRandom = 0;
    end

    properties (Dependent)
        G
    end

    methods
        function G = get.G(obj)
            G = obj.E / (2 * (1 + obj.v));
        end
    end
end
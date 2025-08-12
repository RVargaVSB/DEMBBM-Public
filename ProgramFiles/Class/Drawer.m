classdef Drawer < handle
    properties
        stepCount = 1000000;
        deforScale = 100;

        basicPlots = [1,1,1] % [convergece+deformation, Element stress, Contact stress]
        stressPlot = [1,2,3,4] % [s_x, s_y, s_z, t_xy, t_xz, t_zy]


        f1
        f2
        f3

        Time=[];
        maxU =[];
        minU =[];
        maxA = [];
        minA = [];
        maxV =[];
        minV = [];
        unboundParticles =[];
        transformationActualization = [];
        maxSigmaB = zeros(3,0);





        stepIndex = 0;
        
        stepLength = [];
    end

    methods
        function obj = Record(obj, problem)
            obj.stepIndex = obj.stepIndex + 1;
            obj.stepLength(end+1) = problem.currentState.deltaTime;

            obj.Time(end+1) = problem.solver.timestep;
            obj.maxU(end+1) = max(max(problem.elements.displacement(1:3,:)));
            obj.minU(end+1) = min(min(problem.elements.displacement(1:3,:)));

            obj.maxA(end+1) = max(max(problem.elements.acceleration(1:3,:)));
            obj.minA(end+1) = min(min(problem.elements.acceleration(1:3,:)));

            obj.maxV(end+1) = max(max(problem.elements.velocity(1:3,:)));
            obj.minV(end+1) = min(min(problem.elements.velocity(1:3,:)));

            obj.unboundParticles(end+1) = sum(problem.elements.bounded.tensionVector);
    
            if problem.currentState.updateTrasform == true
                obj.transformationActualization(end+1) = 1;
            else
                obj.transformationActualization(end+1) = 0;
            end

            obj.maxSigmaB(:,end+1) = [max(max(problem.result.sigContact([1,5],:))); min(min(problem.result.sigContact([2,6],:))); max(max(problem.result.sigContact([3,4,7,8],:)))];

            if mod(problem.solver.stepCount,obj.stepCount) == 0
                DrawNow(problem)
            end


        end

        function obj = DrawStart(obj)
             if obj.basicPlots(1)==1
                 obj.f1=figure;
            end
            if obj.basicPlots(2)==1
                 obj.f2=figure;
            end
            if obj.basicPlots(3)==1
                 obj.f3=figure;
            end
        end
    end
    

end

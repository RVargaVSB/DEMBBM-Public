classdef Solver < handle
    properties
        timestep = 1e-04;
        maxTime = 0.02;
        maxStep = inf;
        solverType = "cg+ichol"
        solverOpt = 1e-6 

        limDeg = 1;

        gamma = 0.5;
        beta = 0.25;

        alphaC = 1e-4;
        betaC = 1e-3;

        ram = 5;

        currentTime = 0;
        stepCount = 0;
    end

    properties (Dependent)
        limDelt
    end
    
    methods

    function val = get.limDelt(obj)
        val = cosd(obj.limDeg);
    end

    function problem = Solve(obj,problem)
        
        problem.elements.bounded = BoundedContactsElementAssambler(problem);
        problem.elements.unbounded = UnboundParticlePreparation(problem);
        problem = BoundedMatrixAssembler(problem);

        problem.drawer.DrawStart();
        
        while obj.currentTime < obj.maxTime && obj.stepCount < obj.maxStep
              obj.stepCount = obj.stepCount + 1;
              obj.currentTime = obj.currentTime + problem.currentState.deltaTime;

              problem.elements.externalforces = ExternalForcesEvaluation(problem);
            
              problem = ContactCalculation(problem);

              if ~isempty(problem.elements.unbounded.parisA)
                 [problem] = UnboundStiffnessAssembler(problem); 
              end

              problem=DEMBBMSolverIteration(problem);

              problem.drawer.Record(problem);



        end
        
        DrawNow(problem)

    end

    end
end


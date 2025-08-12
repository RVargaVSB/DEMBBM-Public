function [obj] = CurrentState(problem)
    
        obj.deltaTime = problem.solver.timestep;
        obj.shortStepsCount = 0;
        obj.isShortStep = false;  

end


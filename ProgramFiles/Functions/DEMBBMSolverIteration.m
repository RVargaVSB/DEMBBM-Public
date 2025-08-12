function [problem] = DEMBBMSolverIteration(problem)

    % Transformation matrices
    problem=TransformationMatrixAssembler(problem);

    problem=UnboundForcesCalculation(problem);
    
    problem=FullMatrixAssambler(problem);

    problem=NewMarkSolver(problem);
    

end


function problem = TensionProblem(cords,mainradius,neckRatio,forceX)

    elements = RectangularShape(cords, mainradius);
    elements = ApplyNecking(elements, cords, mainradius, neckRatio);
    problem.geometrySet = elements;

    %% Boundary conditions
    
    dirichlets{1}.position = [find(problem.geometrySet.node(1,:)==0)];
    dirichlets{1}.dofs=[1 1 1 1 1 1];

    problem.dirichlets=Dirichlets(dirichlets,problem);

    %% Force conditions

       problem.forceSetting.type = "time";
       problem.forceSetting.startTime = 0;
       problem.forceSetting.endTime = 0.01;
       problem.force{1}.value = forceX;
       problem.force{1}.direction = 1;
       problem.forceSetting.unload = 0;
       problem.forceSetting.stopOnCrack = 1;
       problem.forceSetting.stopOnCrackCount = 10;

       problem.force{1}.position = [find(problem.geometrySet.node(1,:)==cords(1))];

end


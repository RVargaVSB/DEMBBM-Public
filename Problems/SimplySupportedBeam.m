function problem = SimplySupportedBeam(cords,radius,loadPos,timeModel,loadVector,unloadOptions,stopOnCrack)

    %%%%%%Geometry set

    problem.geometrySet = RectangularShape(cords, radius);

    %%%%%% Boundary conditions
    
    dirichlets{1}.position = [find(and(problem.geometrySet.node(1,:)==0,problem.geometrySet.node(3,:)==0))];
    dirichlets{1}.dofs=[1 0 1 0 0 0];

    dirichlets{2}.position=[find(and(problem.geometrySet.node(1,:)==cords(1),problem.geometrySet.node(3,:)==0))];  %as Rows
    dirichlets{2}.dofs=[1 0 1 0 0 0];

    problem.dirichlets=Dirichlets(dirichlets,problem);

    %%%%%% Force conditions

    if loadPos=="middle"
        problem.force{1}.position = [find(abs(problem.geometrySet.node(1,:) - cords(1)/2) <= radius ...
                             & problem.geometrySet.node(3,:) == cords(3))];
    else
        problem.force{1}.position = find(problem.geometrySet.node(3,:) == cords(3));
    end

    if timeModel=="time"
       problem.forceSetting.type = "time";
       problem.forceSetting.startTime = loadVector(3);
       problem.forceSetting.endTime = loadVector(4);
       problem.forceSetting.stopOnCrack = stopOnCrack;
       problem.force{1}.value = loadVector(1);
       problem.force{1}.direction = loadVector(2);
       problem.forceSetting.stopOnCrackCount = 1;
    else
       problem.forceSetting.type = "step";
       problem.forceSetting.startStep = loadVector(3);
       problem.forceSetting.endStep = loadVector(4);
       problem.force{1}.value = loadVector(1);
       problem.force{1}.direction = loadVector(2); 
       problem.forceSetting.stopOnCrack = stopOnCrack;
       problem.forceSetting.stopOnCrackCount = 1;
    end

    if unloadOptions(1) ~= -1
       problem.forceSetting.unload = 1;
       problem.forceSetting.unloadPause = unloadOptions(1);
       problem.forceSetting.unloadMultiplikator = unloadOptions(2);
    else
       problem.forceSetting.unload = 0;
    end

end

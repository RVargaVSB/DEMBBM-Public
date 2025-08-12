function [externalForces] = ExternalForcesEvaluation(problem)

    externalForces = zeros(size(problem.elements.externalforces));

    for i = 1 : length(problem.force)
            
         if problem.forceSetting.type == "time"
            
            curnTime = problem.solver.currentTime;
            maxTime = problem.forceSetting.endTime;
            startTime = problem.forceSetting.startTime;
    
            if startTime <= curnTime
               forceValue = problem.force{i}.value / length(problem.force{i}.position);
               
               dir = problem.force{i}.direction;
               pos = problem.force{i}.position;

               if curnTime <= maxTime
                   if and(sum(problem.elements.unboundParticles) >= problem.forceSetting.stopOnCrackCount,problem.forceSetting.stopOnCrack==1)
                      externalForces = problem.elements.externalforces;
                   else
                      externalForces(dir,pos) = ones(1,length(pos))*forceValue * (curnTime - startTime) / (maxTime - startTime);
                   end
               else
                   if problem.forceSetting.unload == 1
                     if and(curnTime >= maxTime+problem.forceSetting.unloadPause, curnTime <= maxTime+problem.forceSetting.unloadPause+(startTime-maxTime)/problem.forceSetting.unloadMultiplikator)
                        unloadDuration = (startTime - maxTime) / problem.forceSetting.unloadMultiplikator;
                        unloadStart = maxTime + unloadPause;
                        externalForces(dir,pos) = ones(1,length(pos))*(forceValue * (1 - (curnTime - unloadStart) / unloadDuration));
                     else
                        externalForces = problem.elements.externalforces;
                     end
                        
                    else
                     externalForces(dir,pos) = ones(1,length(pos))*forceValue;
                   end
                end
            end
        end
end
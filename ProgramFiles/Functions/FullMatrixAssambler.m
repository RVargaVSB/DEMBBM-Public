function [problem] = FullMatrixAssambler(problem)

   if problem.currentState.updateTrasform == true
    if ~isempty(problem.elements.unbounded.parisA)

        Ku=problem.elements.unbounded.Localization*problem.elements.unbounded.T'*problem.elements.unbounded.Kdiag*problem.elements.unbounded.T*problem.elements.unbounded.Localization';
        Kb=problem.elements.bounded.Localization*problem.elements.bounded.T'*problem.elements.bounded.Kdiag*problem.elements.bounded.T*problem.elements.bounded.Localization';
        problem.K=Kb+Ku;

    else
        problem.K=problem.elements.bounded.Localization*problem.elements.bounded.T'*problem.elements.bounded.Kdiag*problem.elements.bounded.T*problem.elements.bounded.Localization';
    end
   end

end


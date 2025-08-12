function [problem]=UnboundForcesCalculation(problem)

    if ~isempty(problem.elements.unbounded.parisA)
        
        gapVec=zeros(12,length(problem.elements.unbounded.parisA));
        gapVec(1,:)=-problem.elements.unbounded.gap/2;
        gapVec(7,:)=+problem.elements.unbounded.gap/2;
        gapVec=gapVec(:);

        uA=problem.elements.displacement(:,problem.elements.unbounded.parisA);
        uB=problem.elements.displacement(:,problem.elements.unbounded.parisB);
        u=[uA;uB];        

        % Fn=problem.unbounded.Kdiag*problem.unbounded.T*u(:);
        Fn=problem.elements.unbounded.Kdiag*gapVec(:);

        problem.elements.unbounded.intposXA(problem.elements.unbounded.exist==0)=problem.elements.node(1,problem.elements.unbounded.parisA(problem.elements.unbounded.exist==0));
        problem.elements.unbounded.intposYA(problem.elements.unbounded.exist==0)=problem.elements.node(2,problem.elements.unbounded.parisA(problem.elements.unbounded.exist==0));
        problem.elements.unbounded.intposZA(problem.elements.unbounded.exist==0)=problem.elements.node(3,problem.elements.unbounded.parisA(problem.elements.unbounded.exist==0));
       
        problem.elements.unbounded.intposXB(problem.elements.unbounded.exist==0)=problem.elements.node(1,problem.elements.unbounded.parisB(problem.elements.unbounded.exist==0));
        problem.elements.unbounded.intposYB(problem.elements.unbounded.exist==0)=problem.elements.node(2,problem.elements.unbounded.parisB(problem.elements.unbounded.exist==0));
        problem.elements.unbounded.intposZB(problem.elements.unbounded.exist==0)=problem.elements.node(3,problem.elements.unbounded.parisB(problem.elements.unbounded.exist==0));

        Funbound=problem.elements.unbounded.T'*Fn;

        % Funbound=problem.unbounded.T*Fn;

        problem.elements.unbounded.ForceXA=Funbound(1:12:end);
        problem.elements.unbounded.ForceYA=Funbound(2:12:end);
        problem.elements.unbounded.ForceZA=Funbound(3:12:end);

        problem.elements.unbounded.ForceXB=Funbound(7:12:end);
        problem.elements.unbounded.ForceYB=Funbound(8:12:end);
        problem.elements.unbounded.ForceZB=Funbound(9:12:end);

    end

end
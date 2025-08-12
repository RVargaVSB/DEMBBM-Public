function [problem] = NewMarkSolver(problem)

    u=problem.elements.displacement;

    if ~isempty(problem.elements.unbounded.parisA) 
        Fu=zeros(size(problem.elements.unbounded.Kdiag,2),1);

        Fu(1:12:end)=problem.elements.unbounded.ForceXA;
        Fu(2:12:end)=problem.elements.unbounded.ForceYA;
        Fu(3:12:end)=problem.elements.unbounded.ForceZA;
        Fu(7:12:end)=problem.elements.unbounded.ForceXB;
        Fu(8:12:end)=problem.elements.unbounded.ForceYB;
        Fu(9:12:end)=problem.elements.unbounded.ForceZB;
        
        Fbound=MatrixMultiplicationBackOrder({problem.elements.bounded.Localization,problem.elements.bounded.T',problem.elements.bounded.Kdiag,problem.elements.bounded.T,problem.elements.bounded.Localization'},u(:));
        % Fbound=problem.elements.bounded.Localization*problem.elements.bounded.T'*problem.elements.bounded.Kdiag*problem.elements.bounded.T*problem.elements.bounded.Localization'*u(:);
        Funbound=problem.elements.unbounded.Localization*Fu;
        problem.elements.internalforces=Fbound+Funbound;
    else
         % problem.elements.internalforces=problem.elements.bounded.Localization*problem.elements.bounded.T'*problem.elements.bounded.Kdiag*problem.elements.bounded.T*problem.elements.bounded.Localization'*u(:);
         problem.elements.internalforces=MatrixMultiplicationBackOrder({problem.elements.bounded.Localization,problem.elements.bounded.T',problem.elements.bounded.Kdiag,problem.elements.bounded.T,problem.elements.bounded.Localization'},u(:));
    end

   
    [problem,a,v,du]=SolverFunction(problem,problem.solver.solverType);

    % problem.C=problem.mass_matrix*problem.solver.alfaC+problem.K*problem.solver.betaC;
    % problem.A=problem.K*problem.time.step^2*problem.solver.beta+problem.C*problem.time.step*problem.solver.gama+problem.mass_matrix;
    % problem.A=(problem.A+problem.A')/2;

    % vel=problem.elements.velocity(:);
    % axel=problem.elements.acceleration(:);
    % Fext=problem.elements.externalforces(:);
    % Fint=problem.elements.internalforces(:);
    % F=Fext-Fint;

    % as=problem.time.step*vel+1/2*(I-2*problem.solver.beta*I)*axel*problem.time.step^2;
    % bs=vel+(I-problem.solver.gama*I)*problem.time.step*axel;


    % a=zeros(dof,1);
    
    % B=F-problem.K*as-problem.C*bs;

    % if isempty(problem.LL)
    %      problem.LL=chol(problem.A(~problem.boundarylogvec,~problem.boundarylogvec),"lower");
    % end

    %if isempty(problem.unbounded.parisA)
    %    y=problem.LL\B(~problem.boundarylogvec);
    %    a(~problem.boundarylogvec)=problem.LL' \ y;
    %else
    %a(~problem.boundarylogvec)=problem.A(~problem.boundarylogvec,~problem.boundarylogvec)\B(~problem.boundarylogvec); 
    %end  

    u=problem.elements.displacement+reshape(du,6,[]);

    fu = MatrixMultiplicationBackOrder({problem.elements.bounded.Kdiag,problem.elements.bounded.T,problem.elements.bounded.Localization'},u(:));
    % fu=problem.elements.bounded.Kdiag*problem.elements.bounded.T*problem.elements.bounded.Localization'*u(:);
    fu(1:12:end)=-fu(1:12:end);

    sg=problem.elements.bounded.Kbcross*fu;
    sg=reshape(sg,12,[]);

    sgval=[sg(1,:)+(sg(5,:).^2+sg(6,:).^2).^(1/2); -(sg(1,:)-(sg(5,:).^2+sg(6,:).^2).^(1/2)); (sg(2,:).^2+sg(3,:).^2).^(1/2); abs(sg(4,:))
           sg(7,:)+(sg(11,:).^2+sg(12,:).^2).^(1/2);-(sg(7,:)-(sg(11,:).^2+sg(12,:).^2).^(1/2)); (sg(8,:).^2+sg(9,:).^2).^(1/2); abs(sg(10,:))];

    crack=sgval>problem.elements.bounded.sigLim;
    [~,cpos]=find(crack);
    
    problem.elements.bounded.tensionVectorPrev = problem.elements.bounded.tensionVector;
    problem.elements.bounded.tensionVector(cpos)=0;
       
    problem.elements.unboundParticles(problem.elements.bounded.pairs(1,(problem.elements.bounded.tensionVector==0)))=1;
    problem.elements.unboundParticles(problem.elements.bounded.pairs(2,(problem.elements.bounded.tensionVector==0)))=1;

    % fu=problem.elements.bounded.T'*problem.elements.bounded.Kdiag*problem.elements.bounded.T*problem.elements.bounded.Localization'*u(:);
    fu = MatrixMultiplicationBackOrder({problem.elements.bounded.T',problem.elements.bounded.Kdiag,problem.elements.bounded.T,problem.elements.bounded.Localization'},u(:));

    fu=reshape(-fu,12,[]);

    sij1=problem.elements.bounded.sij1;
    sij2=problem.elements.bounded.sij2;

    Sigvalb1=[fu(1,:).*sij1(1,:);fu(2,:).*sij1(2,:);fu(3,:).*sij1(3,:);(fu(1,:).*sij1(2,:)+fu(2,:).*sij1(1,:))/2
             (fu(3,:).*sij1(2,:)+fu(2,:).*sij1(3,:))/2
             (fu(1,:).*sij1(3,:)+fu(3,:).*sij1(1,:))/2];
    Sigvalb2=[fu(7,:).*sij2(1,:);fu(8,:).*sij2(2,:);fu(9,:).*sij2(3,:);(fu(7,:).*sij2(2,:)+fu(8,:).*sij2(1,:))/2
             (fu(9,:).*sij2(2,:)+fu(8,:).*sij2(3,:))/2
             (fu(7,:).*sij2(3,:)+fu(9,:).*sij2(1,:))/2];

    

    problem.result.nodeContact=(problem.elements.node(:,problem.elements.bounded.pairs(1,:))-(problem.elements.node(:,problem.elements.bounded.pairs(1,:))-problem.elements.node(:,problem.elements.bounded.pairs(2,:)))./2).*problem.elements.bounded.tensionVector;
    problem.result.nodeCrack=((problem.elements.node(:,problem.elements.bounded.pairs(1,:))-(problem.elements.node(:,problem.elements.bounded.pairs(1,:))-problem.elements.node(:,problem.elements.bounded.pairs(2,:)))./2).*(problem.elements.bounded.tensionVector-1)).*-1;
    problem.result.sigContact=diag([1 -1 1 1 1 -1 1 1])*sgval;



    if ~isempty(problem.elements.unbounded.parisA)
        dFu = MatrixMultiplicationBackOrder({problem.elements.unbounded.T',problem.elements.unbounded.Kdiag,problem.elements.unbounded.T,problem.elements.unbounded.Localization'},du(:));
        % dFu=problem.elements.unbounded.T'*problem.elements.unbounded.Kdiag*problem.elements.unbounded.T*problem.elements.unbounded.Localization'*du(:);
        problem.elements.unbounded.ForceXA=problem.elements.unbounded.ForceXA+dFu(1:12:end);
        problem.elements.unbounded.ForceYA=problem.elements.unbounded.ForceYA+dFu(2:12:end);
        problem.elements.unbounded.ForceZA=problem.elements.unbounded.ForceZA+dFu(3:12:end);

        problem.elements.unbounded.ForceXB=problem.elements.unbounded.ForceXB+dFu(7:12:end);
        problem.elements.unbounded.ForceYB=problem.elements.unbounded.ForceYB+dFu(8:12:end);
        problem.elements.unbounded.ForceZB=problem.elements.unbounded.ForceZB+dFu(9:12:end);

        Fu=zeros(size(problem.elements.unbounded.Kdiag,2),1);

        Fu(1:12:end)=problem.elements.unbounded.ForceXA;
        Fu(2:12:end)=problem.elements.unbounded.ForceYA;
        Fu(3:12:end)=problem.elements.unbounded.ForceZA;
        Fu(7:12:end)=problem.elements.unbounded.ForceXB;
        Fu(8:12:end)=problem.elements.unbounded.ForceYB;
        Fu(9:12:end)=problem.elements.unbounded.ForceZB;

        fu=reshape(-Fu,12,[]);

        sij1=problem.elements.unbounded.sij1;
        sij2=problem.elements.unbounded.sij2;

        Sigvalu1=[fu(1,:).*sij1(1,:);fu(2,:).*sij1(2,:);fu(3,:).*sij1(3,:);(fu(1,:).*sij1(2,:)+fu(2,:).*sij1(1,:))/2
                 (fu(3,:).*sij1(2,:)+fu(2,:).*sij1(3,:))/2
                 (fu(1,:).*sij1(3,:)+fu(3,:).*sij1(1,:))/2];
        Sigvalu2=[fu(7,:).*sij2(1,:);fu(8,:).*sij2(2,:);fu(9,:).*sij2(3,:);(fu(7,:).*sij2(2,:)+fu(8,:).*sij2(1,:))/2
                 (fu(9,:).*sij2(2,:)+fu(8,:).*sij2(3,:))/2
                 (fu(7,:).*sij2(3,:)+fu(9,:).*sij2(1,:))/2];
            
        n=size(problem.elements.bounded.pairs(1,:),2)+size(problem.elements.unbounded.parisA,2);
        sig=[Sigvalu1,Sigvalb1,Sigvalu2,Sigvalb2];
        row=[1; 2; 3; 4; 5; 6];
        row=repmat(row,1,2*n);
        col=[repmat(problem.elements.unbounded.parisA,6,1),repmat(problem.elements.bounded.pairs(1,:),6,1),repmat(problem.elements.unbounded.parisB,6,1),repmat(problem.elements.bounded.pairs(2,:),6,1)];
        sig=full(sparse(row,col,sig));
     else
        n=size(problem.elements.bounded.pairs,2);
        sig=[Sigvalb1,Sigvalb2];
        row=[1; 2; 3; 4; 5; 6];
        row=repmat(row,1,2*n);
        col=[repmat(problem.elements.bounded.pairs(1,:),6,1),repmat(problem.elements.bounded.pairs(2,:),6,1)];
        sig=full(sparse(row,col,sig));
    end

    V=repmat(problem.elements.volume,6,1);
    problem.elements.stress=sig./V;
    
    problem.elements.displacement=problem.elements.displacement+reshape(du,6,[]);
    problem.elements.node=problem.elements.nodestart+problem.elements.displacement(1:3,:);
        
    problem.elements.velocity=reshape(v,6,[]);

    problem.elements.acceleration=reshape(a,6,[]);

end


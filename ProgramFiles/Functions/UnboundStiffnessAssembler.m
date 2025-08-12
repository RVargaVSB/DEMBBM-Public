function [problem] = UnboundStiffnessAssembler(problem)

        LambdaVec=zeros(1,size(problem.elements.unbounded.parisA,1));
        L=problem.elements.radius(problem.elements.unbounded.parisA)+problem.elements.radius(problem.elements.unbounded.parisB);
        LambdacontPosVec=problem.elements.radius(problem.elements.unbounded.parisA)<=problem.elements.radius(problem.elements.unbounded.parisB);
        LambdaVec(LambdacontPosVec)=problem.elements.lambdas(problem.elements.unbounded.parisA(LambdacontPosVec));
        LambdaVec(~LambdacontPosVec)=problem.elements.lambdas(problem.elements.unbounded.parisB(~LambdacontPosVec));
        Rb=LambdaVec.*min(problem.elements.radius(problem.elements.unbounded.parisA),problem.elements.radius(problem.elements.unbounded.parisB));
        
        Ab=pi*Rb.^2;

        E = problem.material{1}.E;
        C = problem.material{1}.C0;

        K1_1=E*Ab./L;
        K2_2=ones(size(K1_1))*C;
        K3_3=ones(size(K1_1))*C;

    n_dof=size(K1_1,2)*12;
    P=1:12:n_dof;
    rKc=[P,P,P+1,P+1,P+2,P+2,P+6,P+6,P+7,P+7,P+8,P+8,n_dof];
    cKc=[P,P+6,P+1,P+7,P+2,P+8,P,P+6,P+1,P+7,P+2,P+8,n_dof];
    vKc=[K1_1,-K1_1,K2_2,-K2_2,K3_3,-K3_3,-K1_1,K1_1,-K2_2,K2_2,-K3_3,K3_3,0];
    
    problem.elements.unbounded.Kdiag=sparse(rKc,cKc,vKc);

    u.ca=problem.elements.unbounded.parisA;
    u.cb=problem.elements.unbounded.parisB;
    
    rL=[u.ca*6-5; u.ca*6-4; u.ca*6-3;u.ca*6-2;u.ca*6-1; u.ca*6 ; u.cb*6-5; u.cb*6-4; u.cb*6-3; u.cb*6-2; u.cb*6-1; u.cb*6];
    cL=[P; P+1; P+2; P+3; P+4; P+5; P+6; P+7; P+8; P+9; P+10; P+11];
    vL=true(size(cL));
    
    rm=size(problem.elements.node,2);
    
    problem.elements.unbounded.Localization=sparse(rL,cL,vL,rm*6,length(cL(:)));

    


end


function [problem] = TransformationMatrixAssembler(problem)

       problem.currentState.updateTrasform = true;

    if problem.solver.stepCount == 1
       problem.currentState.updateTrasform = true;
       problem.currentState.LastTransform = 0;
    else

        x1=problem.elements.node(:,problem.elements.bounded.pairs(1,:));
        x2=problem.elements.node(:,problem.elements.bounded.pairs(2,:));
    
        v=(x2-x1);
        normv=v./vecnorm(v);

        A = problem.currentState.LastTransformVector;
        B = (v-(vecnorm(v)-(problem.elements.radius(problem.elements.bounded.pairs(1,:))+problem.elements.radius(problem.elements.bounded.pairs(2,:)))).*normv)/2;

       if and(any( acos( sum(A .* B, 1) ./ (vecnorm(A) .* vecnorm(B)) ) < deg2rad(problem.solver.limDeg) )...
               ,  isequal([problem.elements.unbounded.parisA; problem.elements.unbounded.parisB], problem.currentState.unboundedPairs) )

           problem.currentState.updateTrasform = false;
           problem.currentState.LastTransform = problem.currentState.LastTransform+1;

       end

    end

    if problem.currentState.updateTrasform == true

    x1=problem.elements.node(:,problem.elements.bounded.pairs(1,:));
    x2=problem.elements.node(:,problem.elements.bounded.pairs(2,:));

    v=(x2-x1);
    normv=v./vecnorm(v);

    problem.elements.bounded.sij=(v-(vecnorm(v)-(problem.elements.radius(problem.elements.bounded.pairs(1,:))+problem.elements.radius(problem.elements.bounded.pairs(2,:)))).*normv)/2;
    problem.elements.bounded.sij1=(normv).*problem.elements.radius(problem.elements.bounded.pairs(1,:));
    problem.elements.bounded.sij2=(-normv).*problem.elements.radius(problem.elements.bounded.pairs(2,:));
    
    problem.currentState.LastTransformVector =  problem.elements.bounded.sij;

    AXS=[zeros(1,length(normv(1,:)));-normv(3,:);normv(2,:)];

    normAXS=AXS./vecnorm(AXS);

    cosT=[normv(1,:)];
    sinT=(normv(2,:).^2+normv(3,:).^2).^(1/2);

    T=atan2(sinT,cosT);

    sinT2=sin(T/2);
    cosT2=cos(T/2);

    Q=[cosT2; sinT2.*normAXS(1,:); sinT2.*normAXS(2,:); sinT2.*normAXS(3,:)];

    Q(:,and(Q(1,:)==1,isnan(Q(2,:))))=repmat([1;0;0;0],1,sum(and(Q(1,:)==1,isnan(Q(2,:)))));
    Q(:,and(Q(1,:)==-1,isnan(Q(2,:))))=repmat([0;0;1;0],1,sum(and(Q(1,:)==-1,isnan(Q(2,:)))));

    R11=(1-2*(Q(3,:).^2+Q(4,:).^2));
    R12=2*(Q(2,:).*Q(3,:)-Q(1,:).*Q(4,:));
    R13=2*(Q(2,:).*Q(4,:)+Q(1,:).*Q(3,:));

    R21=2*(Q(2,:).*Q(3,:)+Q(1,:).*Q(4,:));
    R22=(1-2*(Q(2,:).^2+Q(4,:).^2));
    R23=2*(Q(3,:).*Q(4,:)-Q(1,:).*Q(2,:));

    R31=2*(Q(2,:).*Q(4,:)-Q(1,:).*Q(3,:));
    R32=2*(Q(3,:).*Q(4,:)+Q(1,:).*Q(2,:));
    R33=(1-2*(Q(2,:).^2+Q(3,:).^2));

    n_dof=size(problem.elements.bounded.Kdiag,2);
    P=1:12:n_dof;
    
    rTb=[P P P
         P+1 P+1 P+1
         P+2 P+2 P+2
         P+3 P+3 P+3
         P+4 P+4 P+4
         P+5 P+5 P+5
         P+6 P+6 P+6
         P+7 P+7 P+7
         P+8 P+8 P+8
         P+9 P+9 P+9
         P+10 P+10 P+10
         P+11 P+11 P+11];
    cTb=[P P+1 P+2
         P P+1 P+2
         P P+1 P+2
         P+3 P+4 P+5
         P+3 P+4 P+5
         P+3 P+4 P+5
         P+6 P+7 P+8
         P+6 P+7 P+8
         P+6 P+7 P+8
         P+9 P+10 P+11
         P+9 P+10 P+11
         P+9 P+10 P+11];
    vTb=[R11 R12 R13
         R21 R22 R23
         R31 R32 R33
         R11 R12 R13
         R21 R22 R23
         R31 R32 R33
         R11 R12 R13
         R21 R22 R23
         R31 R32 R33
         R11 R12 R13
         R21 R22 R23
         R31 R32 R33];


    problem.elements.bounded.T=sparse(rTb,cTb,vTb);
    problem.elements.bounded.T=problem.elements.bounded.T';

    if ~isempty(problem.elements.unbounded.parisA)
    
        x1=problem.elements.node(:,problem.elements.unbounded.parisA);
        x2=problem.elements.node(:,problem.elements.unbounded.parisB);
    
        v=(x2-x1);
        normv=v./vecnorm(v);

        problem.elements.unbounded.sij=(v-(vecnorm(v)-(problem.elements.radius(problem.elements.unbounded.parisA)+problem.elements.radius(problem.elements.unbounded.parisB))).*normv)/2;
        problem.elements.unbounded.sij1=normv.*problem.elements.radius(problem.elements.unbounded.parisA);
        problem.elements.unbounded.sij2=normv.*problem.elements.radius(problem.elements.unbounded.parisB);
      
            AXS=[zeros(1,length(normv(1,:)));-normv(3,:);normv(2,:)];

    normAXS=AXS./vecnorm(AXS);

    cosT=[normv(1,:)];
    sinT=(normv(2,:).^2+normv(3,:).^2).^(1/2);

    T=atan2(sinT,cosT);

    sinT2=sin(T/2);
    cosT2=cos(T/2);

    Q=[cosT2; sinT2.*normAXS(1,:); sinT2.*normAXS(2,:); sinT2.*normAXS(3,:)];

    Q(:,and(Q(1,:)==1,isnan(Q(2,:))))=repmat([1;0;0;0],1,sum(and(Q(1,:)==1,isnan(Q(2,:)))));
    Q(:,and(Q(1,:)==-1,isnan(Q(2,:))))=repmat([0;0;1;0],1,sum(and(Q(1,:)==-1,isnan(Q(2,:)))));

    R11=(1-2*(Q(3,:).^2+Q(4,:).^2));
    R12=2*(Q(2,:).*Q(3,:)-Q(1,:).*Q(4,:));
    R13=2*(Q(2,:).*Q(4,:)+Q(1,:).*Q(3,:));

    R21=2*(Q(2,:).*Q(3,:)+Q(1,:).*Q(4,:));
    R22=(1-2*(Q(2,:).^2+Q(4,:).^2));
    R23=2*(Q(3,:).*Q(4,:)-Q(1,:).*Q(2,:));

    R31=2*(Q(2,:).*Q(4,:)-Q(1,:).*Q(3,:));
    R32=2*(Q(3,:).*Q(4,:)+Q(1,:).*Q(2,:));
    R33=(1-2*(Q(2,:).^2+Q(3,:).^2));

    n_dof=size(problem.elements.unbounded.Kdiag,2);
    P=1:12:n_dof;
    
    rTu=[P P P
         P+1 P+1 P+1
         P+2 P+2 P+2
         P+3 P+3 P+3
         P+4 P+4 P+4
         P+5 P+5 P+5
         P+6 P+6 P+6
         P+7 P+7 P+7
         P+8 P+8 P+8
         P+9 P+9 P+9
         P+10 P+10 P+10
         P+11 P+11 P+11];
    cTu=[P P+1 P+2
         P P+1 P+2
         P P+1 P+2
         P+3 P+4 P+5
         P+3 P+4 P+5
         P+3 P+4 P+5
         P+6 P+7 P+8
         P+6 P+7 P+8
         P+6 P+7 P+8
         P+9 P+10 P+11
         P+9 P+10 P+11
         P+9 P+10 P+11];
    vTu=[R11 R12 R13
         R21 R22 R23
         R31 R32 R33
         R11 R12 R13
         R21 R22 R23
         R31 R32 R33
         R11 R12 R13
         R21 R22 R23
         R31 R32 R33
         R11 R12 R13
         R21 R22 R23
         R31 R32 R33];

        problem.elements.unbounded.T=sparse(rTu,cTu,vTu);
        problem.elements.unbounded.T=problem.elements.unbounded.T';
    else

        Tu=[];
        sij=[];

    end

    end


end
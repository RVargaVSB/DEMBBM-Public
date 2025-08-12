function [problem] = BoundedMatrixAssembler(problem)

        LambdaVec=zeros(1,size(problem.elements.bounded.pairs,1));
        L=problem.elements.radius(problem.elements.bounded.pairs(1,:))+problem.elements.radius(problem.elements.bounded.pairs(2,:));
        LambdacontPosVec=problem.elements.radius(problem.elements.bounded.pairs(1,:))<=problem.elements.radius(problem.elements.bounded.pairs(2,:));
        LambdaVec(LambdacontPosVec)=problem.elements.lambdas(problem.elements.bounded.pairs(1,LambdacontPosVec));
        LambdaVec(~LambdacontPosVec)=problem.elements.lambdas(problem.elements.bounded.pairs(2,~LambdacontPosVec));
        Rb=LambdaVec.*min(problem.elements.radius(problem.elements.bounded.pairs(1,:)),problem.elements.radius(problem.elements.bounded.pairs(2,:)));
        Ab=pi*Rb.^2;
        Ib=pi*Rb.^4/4;
        Ix=2*Ib;
        Wy=pi*Rb.^3/4;
        Wx=Wy*2;
        E = problem.material{1}.E;
        G = problem.material{1}.G;

        fik=12*E.*Ib./(G*Ab.*L.^2);
        
        K1_1=E*Ab./L;
        K2_2=12*E.*Ib./((1+fik).*L.^3);
        K4_4=G.*Ix./L;
        K5_5=(4+fik)*E.*Ib./((1+fik).*L);
        K2_6=6*E.*Ib./((1+fik).*L.^2);
        K4_11=(2-fik).*E.*Ib./((1+fik).*L); 

        Kc1_1=1./Ab;
        Kc2_2=4/3./Ab;
        Kc3_3=4/3./Ab;
        Kc4_4=1./Wy;
        Kc5_5=1./Wy;
        Kc6_6=1./Wx;

        meq=(problem.elements.mass(problem.elements.bounded.pairs(1,:)).*problem.elements.mass(problem.elements.bounded.pairs(2,:)))./(problem.elements.mass(problem.elements.bounded.pairs(1,:))+problem.elements.mass(problem.elements.bounded.pairs(2,:)));
        Jeq=(problem.elements.inertia(problem.elements.bounded.pairs(1,:)).*problem.elements.inertia(problem.elements.bounded.pairs(2,:)))./(problem.elements.inertia(problem.elements.bounded.pairs(1,:))+problem.elements.inertia(problem.elements.bounded.pairs(2,:)));

        n_dof=size(K1_1,2)*12;
        n_dofr=size(K1_1,2)*12;
        P=1:12:n_dof;
        Ps=1:12:n_dofr;

        rKc=[P,P,P+1,P+1,P+1,P+1,P+2,P+2,P+2,P+2,P+3,P+3,P+4,P+4,P+4,P+4,P+5,P+5,P+5,P+5,P+6,P+6,P+7,P+7,P+7,P+7,P+8,P+8,P+8,P+8,P+9,P+9,P+10,P+10,P+10,P+10,P+11,P+11,P+11,P+11];
        cKc=[P,P+6,P+1,P+5,P+7,P+11,P+2,P+4,P+8,P+10,P+3,P+9,P+2,P+4,P+8,P+10,P+1,P+5,P+7,P+11,P,P+6,P+1,P+5,P+7,P+11,P+2,P+4,P+8,P+10,P+3,P+9,P+2,P+4,P+8,P+10,P+1,P+5,P+7,P+11];
        vKc=[K1_1,-K1_1,K2_2,K2_6,-K2_2,K2_6,K2_2,-K2_6,-K2_2,-K2_6,K4_4,-K4_4,-K2_6,K5_5,K2_6,K4_11,K2_6,K5_5,-K2_6,K4_11,-K1_1,K1_1,-K2_2,-K2_6,K2_2,-K2_6,-K2_2,K2_6,K2_2,K2_6,-K4_4,K4_4,-K2_6,K4_11,K2_6,K5_5,K2_6,K4_11,-K2_6,K5_5];
        
        rKcr=[Ps,Ps+1,Ps+2,Ps+3,Ps+4,Ps+5,Ps+6,Ps+7,Ps+8,Ps+9,Ps+10,Ps+11];
        cKcr=[Ps,Ps+1,Ps+2,Ps+3,Ps+4,Ps+5,Ps+6,Ps+7,Ps+8,Ps+9,Ps+10,Ps+11];
        vKcr=[Kc1_1,Kc2_2,Kc3_3,Kc4_4,Kc5_5,Kc6_6,Kc1_1,Kc2_2,Kc3_3,Kc4_4,Kc5_5,Kc6_6];
        
        problem.elements.bounded.Kbcross=sparse(rKcr,cKcr,vKcr);
        
        problem.elements.bounded.Kdiag=sparse(rKc,cKc,vKc);
        
        dofs=size(problem.elements.node,2)*6;
        dpos=1:6:dofs;
        rcM=[dpos, dpos+1, dpos+2, dpos+3, dpos+4, dpos+5];
        vM=[problem.elements.mass, problem.elements.mass, problem.elements.mass, problem.elements.inertia,problem.elements.inertia, problem.elements.inertia];
        problem.mass_matrix=sparse(rcM,rcM,vM);
    
        c.as=problem.elements.bounded.pairs(1,:);
        c.bs=problem.elements.bounded.pairs(2,:);
        
        rL=[c.as*6-5; c.as*6-4; c.as*6-3; c.as*6-2; c.as*6-1; c.as*6; c.bs*6-5; c.bs*6-4; c.bs*6-3; c.bs*6-2; c.bs*6-1; c.bs*6];
        cL=[P; P+1; P+2; P+3; P+4; P+5; P+6; P+7; P+8; P+9; P+10; P+11];
        vL=true(size(cL));
        
        problem.elements.bounded.Localization=sparse(rL,cL,vL);
end


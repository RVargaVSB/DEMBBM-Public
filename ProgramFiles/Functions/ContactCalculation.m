function [problem] = ContactCalculation(problem)

        unbound=problem.elements.unbounded;
        problem.currentState.unboundedPairs = [problem.elements.unbounded.parisA; problem.elements.unbounded.parisB];

        if sum(problem.elements.bounded.tensionVectorPrev) > sum(problem.elements.bounded.tensionVector)
            n_dof=size(problem.elements.bounded.tensionVector,2);
            DG = spdiags(problem.elements.bounded.tensionVector',0,n_dof,n_dof);
            RMK = kron(DG, eye(12));   
            problem.elements.bounded.Kdiag=RMK*problem.elements.bounded.Kdiag*RMK;
        end

        problem.elements.bounded.tensionVectorPrev = problem.elements.bounded.tensionVector;

        pos=find(problem.elements.unboundParticles);

        if ~isempty(pos)
            nc=length(pos);
            xuc=problem.elements.node(1,pos);
            yuc=problem.elements.node(2,pos);
            zuc=problem.elements.node(3,pos);
            xall=problem.elements.node(1,:);
            yall=problem.elements.node(2,:);
            zall=problem.elements.node(3,:);
            np=length(xall);
            RMall=repmat(problem.elements.radius,nc,1);
            Rc=problem.elements.radius(pos);
            RMc=repmat(Rc,np,1);
            Rd2=RMall+RMc';
            AUXa=repmat(xall,nc,1);
            AUYa=repmat(yall,nc,1);
            AUZa=repmat(zall,nc,1);
            AUXc=repmat(xuc,np,1);
            AUYc=repmat(yuc,np,1);
            AUZc=repmat(zuc,np,1);
            GapMatrix=round(((AUXa-AUXc').^2+(AUYa-AUYc').^2+(AUZa-AUZc').^2).^(1/2),10);
            ContactMatrix=GapMatrix<=Rd2;
            [p1,p2]=find(ContactMatrix);
            cc=[pos(p1); p2'];
            p1=p1';
            p2=p2';
            p1(:,cc(1,:)==cc(2,:))=[];
            p2(:,cc(1,:)==cc(2,:))=[];
            cc(:,cc(1,:)==cc(2,:))=[];
            idx = sub2ind(size(GapMatrix),p1,p2);
            gv=GapMatrix(idx)-Rd2(idx);
            cc=sort(cc);
            cb=[problem.elements.bounded.pairs(1,:); problem.elements.bounded.pairs(2,:)];
            cb(:,problem.elements.bounded.tensionVector==0)=[];
            [~,~,ib] = intersect(cb',cc','rows');
            cc(:,ib)=[];
            gv(ib)=[];
            cc=[cc; gv];
            cc = unique(cc.', 'rows', 'stable').';
        
            problem.elements.unbounded.gap=cc(3,:);
            problem.elements.unbounded.parisA=cc(1,:);
            problem.elements.unbounded.parisB=cc(2,:);
        
            [problem.elements.unbounded.intposXA,problem.elements.unbounded.intposXB,problem.elements.unbounded.intposYA,problem.elements.unbounded.intposYB,problem.elements.unbounded.intposZA,problem.elements.unbounded.intposZB,problem.elements.unbounded.exist]=deal(zeros(size(problem.elements.unbounded.parisA)));
    else
            [problem.elements.unbounded.parisB,problem.elements.unbounded.parisA,problem.elements.unbounded.intposXA,problem.elements.unbounded.intposXB,problem.elements.unbounded.intposYA,problem.elements.unbounded.intposYB,problem.elements.unbounded.intposZA,problem.elements.unbounded.intposZB,problem.elements.unbounded.exist]=deal([]);   
    end


    if  ~isempty(problem.elements.unbounded.parisA)
        if ~isempty(unbound.parisA)
        ccp=[unbound.parisA; unbound.parisB];
        cc=cc(1:2,:);
        [v,ib]=ismember(ccp.',cc.','rows');
        ib(ib==0)=[];
        problem.elements.unbounded.intposXA(ib)=unbound.intposXA(v);
        problem.elements.unbounded.intposXB(ib)=unbound.intposXB(v);
        problem.elements.unbounded.intposYA(ib)=unbound.intposYA(v);
        problem.elements.unbounded.intposYB(ib)=unbound.intposYB(v);
        problem.elements.unbounded.intposZA(ib)=unbound.intposZA(v);
        problem.elements.unbounded.intposZB(ib)=unbound.intposZB(v);
        problem.elements.unbounded.exist(ib)=unbound.exist(v);
        end
    end
   

end

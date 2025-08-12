function [bounded] = BoundedContactsElementAssambler(problem)

    RAM=problem.solver.ram*1024*1024*1024/8;

    n_e=size(problem.elements.node,2);

    if n_e<sqrt(RAM)

       xe=problem.elements.node(1,:);
       ye=problem.elements.node(2,:);
       ze=problem.elements.node(3,:);
       Rad=repmat(problem.elements.radius,n_e,1);
       Rad2=Rad+Rad';
       AUX=repmat(xe,n_e,1);
       AUY=repmat(ye,n_e,1);
       AUZ=repmat(ze,n_e,1);
       GapMatrix=round(((AUX-AUX').^2+(AUY-AUY').^2+(AUZ-AUZ').^2).^(1/2),10);
       ContactMatrix=GapMatrix<=Rad2;
       ContactMatrix = ContactMatrix - eye(size(problem.elements.node,2));
       ContactMatrix=triu(ContactMatrix);
       [as,bs]=find(ContactMatrix);
              
       as=as';
       bs=bs';
       
       bounded.pairs=[as;bs];

    else
        as=[];
        bs=[];

        maxtics=floor(sqrt(RAM));
        nummat=ceil(n_e/maxtics);

        for i = 1:nummat
            for j=i:nummat
                istart=(i-1)*maxtics+1;
                if i == nummat
                    iend=n_e;
                else
                    iend=i*maxtics;
                end
                jstart=(j-1)*maxtics+1;
                if j == nummat
                    jend=n_e;
                else
                    jend=j*maxtics;
                end

                XI=problem.elements.node(:,istart:iend)';
                XJ=problem.elements.node(:,jstart:jend)';
                RI=problem.elements.radius(istart:iend);
                RJ=problem.elements.radius(jstart:jend);
                n_j=iend-istart+1;
                n_i=jend-jstart+1;
                RI=repmat(RI,n_i,1);
                RJ=repmat(RJ,n_j,1);
                RC=RI+RJ';
                Gapmatrix=pdist2(XI,XJ);
                ContactMatrix=Gapmatrix<=RC';
                if i==j
                    ContactMatrix = ContactMatrix - eye(n_i);
                    ContactMatrix=triu(ContactMatrix);
                end
                [ads,bds]=find(ContactMatrix);
                ads=istart+ads-1;
                bds=jstart+bds-1;
                as=[as ads'];
                bs=[bs bds'];  
            end
        end
        
        bounded.pairs=[as;bs];
    end

    bounded.tensionVector=ones(size(as));
    bounded.tensionVectorPrev =  bounded.tensionVector;
    
    n = numel(as);                                       
    matIds = problem.elements.material(as);  
    
    bounded.sigmaMax = arrayfun(@(i) normrnd(1, problem.material{matIds(i)}.sigmaRandom) * problem.material{matIds(i)}.fsmax, (1:n)');
    bounded.sigmaMin = arrayfun(@(i) normrnd(1, problem.material{matIds(i)}.sigmaRandom) * problem.material{matIds(i)}.fsmin, (1:n)');
    bounded.tauMax   = arrayfun(@(i) normrnd(1, problem.material{matIds(i)}.sigmaRandom) * problem.material{matIds(i)}.ftmax, (1:n)');

    bounded.sigLim = [bounded.sigmaMax'; bounded.sigmaMax'; bounded.tauMax'; bounded.tauMax'; bounded.sigmaMax'; bounded.sigmaMax'; bounded.tauMax'; bounded.tauMax'];
end
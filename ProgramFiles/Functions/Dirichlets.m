function [set] = Dirichlets(dirichlets,problem)
        
    logvec=false(size(problem.geometrySet.node,2),6);
    for i = 1:length(dirichlets)
        points=dirichlets{i}.position;
        dir=dirichlets{i}.dofs;
        
        dir=repmat(dir,size(points,2),1);
    
        logvec(points,:)=dir;
    end
    logvec=logvec';
    set=logvec(:);

end


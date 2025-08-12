function [problem,a,v,du] = SolverFunction(problem,type)


    if type == "full"
        dof=size(problem.K,1);
        I=speye(dof);
        Z=sparse(dof,dof);
        
        if problem.currentState.updateTrasform
            problem.C = problem.mass_matrix*problem.solver.alphaC+problem.K*problem.solver.betaC;
            A = problem.K*problem.currentState.deltaTime^2*problem.solver.beta+problem.C*problem.currentState.deltaTime*problem.solver.gamma+problem.mass_matrix;
            problem.A = (A+A')/2;
        end

        
        vel=problem.elements.velocity(:);
        axel=problem.elements.acceleration(:);
        Fext=problem.elements.externalforces(:);
        Fint=problem.elements.internalforces(:);
        F=Fext-Fint;
        
        as=problem.currentState.deltaTime*vel+1/2*(I-2*problem.solver.beta*I)*axel*problem.currentState.deltaTime^2;
        bs=vel+(I-problem.solver.gamma*I)*problem.currentState.deltaTime*axel;
        
        a=zeros(dof,1);
        
        B=F-problem.K*as-problem.C*bs;
    
        a(~problem.dirichlets)=problem.A(~problem.dirichlets,~problem.dirichlets)\B(~problem.dirichlets); 
    
        v=bs+problem.currentState.deltaTime*problem.solver.gamma*a;
        du=as+problem.currentState.deltaTime^2*problem.solver.beta*a;
    end

    if type == "choll"

        dof=size(problem.K,1);
        I=speye(dof);
        Z=sparse(dof,dof);
        
        vel=problem.elements.velocity(:);
        axel=problem.elements.acceleration(:);
        Fext=problem.elements.externalforces(:);
        Fint=problem.elements.internalforces(:);
        F=Fext-Fint;

        if problem.currentState.updateTrasform
            problem.C = problem.mass_matrix*problem.solver.alphaC+problem.K*problem.solver.betaC;
            A = problem.K*problem.currentState.deltaTime^2*problem.solver.beta+problem.C*problem.currentState.deltaTime*problem.solver.gamma+problem.mass_matrix;
            A = (A+A')/2;
            problem.LL = chol(A(~problem.dirichlets, ~problem.dirichlets), 'lower');
        end

        as=problem.currentState.deltaTime*vel+1/2*(I-2*problem.solver.beta*I)*axel*problem.currentState.deltaTime^2;
        bs=vel+(I-problem.solver.gamma*I)*problem.currentState.deltaTime*axel;
        
        a=zeros(dof,1);
        
        B=F-problem.K*as-problem.C*bs;

        rhs = B(~problem.dirichlets);
        L = problem.LL;
        y = L \ rhs;
        a(~problem.dirichlets) = L' \ y;

        v=bs+problem.currentState.deltaTime*problem.solver.gamma*a;
        du=as+problem.currentState.deltaTime^2*problem.solver.beta*a;
    end

    if type == "explicit-Euler"

        dof = size(problem.K,1);

        axel=problem.elements.acceleration(:);
        vel  = problem.elements.velocity(:);
        Fext = problem.elements.externalforces(:);
        Fint = problem.elements.internalforces(:);

        if problem.currentState.updateTrasform
            problem.C = problem.mass_matrix * problem.solver.alphaC + problem.K * problem.solver.betaC;
        end

        Fdump = problem.C * vel;
        F = Fext - Fint - Fdump;

        a = zeros(dof,1);
        a(~problem.dirichlets) = F(~problem.dirichlets) ./ diag(problem.mass_matrix(~problem.dirichlets, ~problem.dirichlets));

        v = vel + problem.currentState.deltaTime * a;
        du = problem.currentState.deltaTime * v;
    end

     if type == "explicit-CDM"

        dof = size(problem.K,1);

        if problem.solver.stepCount == 1 
           u_old = problem.elements.displacement(:);
        else
           u_old = problem.currentState.u_old;
        end

        u = problem.elements.displacement(:);
        axel=problem.elements.acceleration(:);
        vel  = problem.elements.velocity(:);
        Fext = problem.elements.externalforces(:);
        Fint = problem.elements.internalforces(:);

        if problem.currentState.updateTrasform
            problem.C = problem.mass_matrix * problem.solver.alphaC + problem.K * problem.solver.betaC;
        end

        Fdump = problem.C * vel;
        F = Fext - Fint - Fdump;

        a = zeros(dof,1);
        a(~problem.dirichlets) = F(~problem.dirichlets) ./ diag(problem.mass_matrix(~problem.dirichlets, ~problem.dirichlets));

        u_new = 2 * u - u_old + problem.currentState.deltaTime^2 * a;
        problem.currentState.u_old = u;
        
        du = u_new - u;
        v = (u_new - u_old) / (2*problem.currentState.deltaTime);
    end

    if type == "cg"
        dof = size(problem.K,1);
        I = speye(dof);
    
        vel  = problem.elements.velocity(:);
        axel = problem.elements.acceleration(:);
        Fext = problem.elements.externalforces(:);
        Fint = problem.elements.internalforces(:);
        F = Fext - Fint;
    
        if problem.currentState.updateTrasform
            problem.C = problem.mass_matrix * problem.solver.alphaC + problem.K * problem.solver.betaC;
            A = problem.K * problem.currentState.deltaTime^2 * problem.solver.beta ...
              + problem.C * problem.currentState.deltaTime * problem.solver.gamma ...
              + problem.mass_matrix;
            A = (A + A') / 2;
            problem.A_iter = A;
    
        end
    
        as = problem.currentState.deltaTime * vel + ...
             0.5 * (I - 2 * problem.solver.beta * I) * axel * problem.currentState.deltaTime^2;
        bs = vel + (I - problem.solver.gamma * I) * problem.currentState.deltaTime * axel;
    
        B = F - problem.K * as - problem.C * bs;
    
        a = zeros(dof,1);
        free = ~problem.dirichlets;

        A=problem.A_iter;
    
        [asol, ~] = pcg(@(x) A(free,free)*x, B(free), problem.solver.solverOpt, 300);
        a(free) = asol;
    
        v = bs + problem.currentState.deltaTime * problem.solver.gamma * a;
        du = as + problem.currentState.deltaTime^2 * problem.solver.beta * a;
    end
    
    if type == "cg+ichol"
        dof = size(problem.K,1);
        I = speye(dof);
        free = ~problem.dirichlets;
    
        vel  = problem.elements.velocity(:);
        axel = problem.elements.acceleration(:);
        Fext = problem.elements.externalforces(:);
        Fint = problem.elements.internalforces(:);
        F = Fext - Fint;
    
        if problem.currentState.updateTrasform
            problem.C = problem.mass_matrix * problem.solver.alphaC + problem.K * problem.solver.betaC;
            A = problem.K * problem.currentState.deltaTime^2 * problem.solver.beta ...
              + problem.C * problem.currentState.deltaTime * problem.solver.gamma ...
              + problem.mass_matrix;
            A = (A + A') / 2;
            problem.A_iter = A;
            problem.M1 = ichol(A(free,free));
        end
    
        as = problem.currentState.deltaTime * vel + ...
             0.5 * (I - 2 * problem.solver.beta * I) * axel * problem.currentState.deltaTime^2;
        bs = vel + (I - problem.solver.gamma * I) * problem.currentState.deltaTime * axel;
    
        B = F - problem.K * as - problem.C * bs;
    
        a = zeros(dof,1);

        A=problem.A_iter;
    
        [asol, ~] = pcg(@(x) A(free,free)*x, B(free), problem.solver.solverOpt, 300,problem.M1,problem.M1');
        a(free) = asol;
    
        v = bs + problem.currentState.deltaTime * problem.solver.gamma * a;
        du = as + problem.currentState.deltaTime^2 * problem.solver.beta * a;
    end

end


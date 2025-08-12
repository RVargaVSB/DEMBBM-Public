function [elements] = SetElements(problem)



n = size(problem.geometrySet.node,2);

elements.lambdas = problem.lambda(problem.geometrySet.prop);

elements.node=problem.geometrySet.node;
elements.nodestart=problem.geometrySet.node;
elements.radius = problem.geometrySet.radius;


elements.volume=4/3*pi*elements.radius.^3;

roTable = cellfun(@(m) m.ro, problem.material);
elements.mass = elements.volume .* roTable(problem.geometrySet.propMat);
elements.inertia=2/5*elements.mass.*elements.radius.^2;
elements.material = problem.geometrySet.propMat;

elements.velocity=zeros(6,n);
elements.acceleration=zeros(6,n);
elements.displacement=zeros(6,n);

elements.internalforces=zeros(6,n);
elements.externalforces=zeros(6,n);

elements.stress=zeros(6,n);

elements.unboundParticles=zeros(1,n);


end


function [] = DrawNow(problem)

     if problem.drawer.basicPlots(1) == 1
        figure(problem.drawer.f1);
        if problem.solver.timestep>0
           clf
        end
        subplot(2,2,1)
        scatter3(problem.elements.nodestart(1,:)+problem.elements.displacement(1,:)*problem.drawer.deforScale,problem.elements.nodestart(2,:)+problem.elements.displacement(2,:)*problem.drawer.deforScale,problem.elements.nodestart(3,:)+problem.elements.displacement(3,:)*problem.drawer.deforScale,problem.elements.radius*5000,"red","filled")
        colorbar
        axis equal
    
        subplot(2,2,2)
        hold on    
        plot(problem.drawer.maxU)
        plot(problem.drawer.minU)
        subplot(2,2,3)
        hold on
        plot(problem.drawer.maxV)
        plot(problem.drawer.minV)
        subplot(2,2,4)
        hold on
        plot(problem.drawer.maxA)
        plot(problem.drawer.minA)
        drawnow
        hold off
    end

    
    if problem.drawer.basicPlots(2)==1

% Vykreslení napětí v elementech
    figure(problem.drawer.f2);
        if problem.solver.timestep>0
           clf
        end
    subplot(2,2,1)   
    scatter3(problem.elements.node(1,:),problem.elements.node(2,:),problem.elements.node(3,:),problem.elements.radius*5000,problem.elements.stress(problem.drawer.stressPlot(1),:)',"filled");
    colorbar
    colormap turbo
    axis equal
    subplot(2,2,2)
    scatter3(problem.elements.node(1,:),problem.elements.node(2,:),problem.elements.node(3,:),problem.elements.radius*5000,problem.elements.stress(problem.drawer.stressPlot(2),:)',"filled");
    colorbar
    colormap turbo
    axis equal
    subplot(2,2,3)
    scatter3(problem.elements.node(1,:),problem.elements.node(2,:),problem.elements.node(3,:),problem.elements.radius*5000,problem.elements.stress(problem.drawer.stressPlot(3),:)',"filled");
    colorbar
    colormap turbo
    axis equal
    subplot(2,2,4)
    scatter3(problem.elements.node(1,:),problem.elements.node(2,:),problem.elements.node(3,:),problem.elements.radius*5000,problem.elements.stress(problem.drawer.stressPlot(4),:)',"filled");
    colorbar
    colormap turbo
    hold on
    axis equal
    drawnow
    end

    if problem.drawer.basicPlots(3)==1

    figure(problem.drawer.f3);
        if problem.solver.timestep>0
           clf
        end
    hold on
    subplot(2,2,1)
    scatter3(problem.result.nodeContact(1,:),problem.result.nodeContact(2,:),problem.result.nodeContact(3,:),50,max(problem.result.sigContact([2,6],:)),"filled");
    colorbar
    view(3)
    axis equal
    subplot(2,2,2)
    scatter3(problem.result.nodeContact(1,:),problem.result.nodeContact(2,:),problem.result.nodeContact(3,:),50,max(problem.result.sigContact([3,4,7,8],:)),"filled");
    colorbar
    view(3)
    axis equal
    subplot(2,2,3)
    scatter3(problem.result.nodeContact(1,:),problem.result.nodeContact(2,:),problem.result.nodeContact(3,:),50,min(problem.result.sigContact([1,2],:)),"filled");
    colorbar
    view(3)
    axis equal
    subplot(2,2,4)
    scatter3(problem.result.nodeCrack(1,:),problem.result.nodeCrack(2,:),problem.result.nodeCrack(3,:),50,"red","filled");
    colorbar
    view(3)
    axis equal

    drawnow
    end


end


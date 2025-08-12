figure
hold on

% === Všechny kontakty modře ===
scatter3(problem.result.nodeContact(1,:), ...
         problem.result.nodeContact(2,:), ...
         problem.result.nodeContact(3,:), ...
         500, 'blue', 'filled');

% === Trhliny červeně (bez nultého bodu, a jen nenulové) ===
crackNodes = problem.result.nodeCrack(:,2:end);  % vynecháme první
nonzeroIdx = any(crackNodes ~= 0, 1);            % sloupce kde je aspoň něco nenulové
crackNodes = crackNodes(:, nonzeroIdx);          % jen tyto

scatter3(crackNodes(1,:), crackNodes(2,:), crackNodes(3,:), ...
         500, 'red', 'filled');
axis off           % skryje osy a popisky
box off            % odstraní rámeček kolem grafu (pokud tam je)
view(3)
axis equal
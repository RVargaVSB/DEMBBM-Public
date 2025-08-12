function [problem] = Konzolka()

%% Parametry
r = 0.025;
d = 2*r;
r_small = r * (sqrt(2) - 1);

% Sloup
sloup_x = 0.6;
sloup_y = 0.6;
sloup_z = 3;

% Konzolka
konzolka_x = 0.6;
konzolka_y = sloup_y;
konzolka_z = 0.4;
konzolka_z_offset = 2;

%% === SLOUP ===
[xs, ys, zs] = ndgrid(0:d:sloup_x, 0:d:sloup_y, 0:d:sloup_z);
Xs_sloup = xs(:);
Ys_sloup = ys(:);
Zs_sloup = zs(:);

[xs_m, ys_m, zs_m] = ndgrid(d/2:d:sloup_x-d/2, d/2:d:sloup_y-d/2, d/2:d:sloup_z-d/2);
Xs_small_sloup = xs_m(:);
Ys_small_sloup = ys_m(:);
Zs_small_sloup = zs_m(:);

%% === KONZOLKA (rozšířená výška, ořízneme spodní klín) ===
[xk_full, yk_full, zk_full] = ndgrid( ...
    sloup_x : d : sloup_x + konzolka_x, ...
    0 : d : konzolka_y, ...
    konzolka_z_offset : d : konzolka_z_offset + 2 * konzolka_z);
Xk_all = xk_full(:);
Yk_all = yk_full(:);
Zk_all = zk_full(:);

% Přímka mezi body: (sloup_x, konzolka_z_offset) → (sloup_x + konzolka_x, konzolka_z_offset + konzolka_z)
relX = (Xk_all - sloup_x) ./ konzolka_x;
Zline = konzolka_z_offset + relX * konzolka_z;

% Zachovat body NAD nebo NA přímce
mask = Zk_all >= Zline;
Xk = Xk_all(mask);
Yk = Yk_all(mask);
Zk = Zk_all(mask);

% Vyřazení konzolových kuliček, které jsou uvnitř sloupu
idx_konzolka_valid = Xk > sloup_x;
Xk = Xk(idx_konzolka_valid);
Yk = Yk(idx_konzolka_valid);
Zk = Zk(idx_konzolka_valid);

%% === MALÉ KULIČKY do konzolky (se stejným oříznutím)
[xk_s, yk_s, zk_s] = ndgrid( ...
    sloup_x + r : d : sloup_x + konzolka_x - r, ...
    d/2 : d : konzolka_y - d/2, ...
    konzolka_z_offset + r : d : konzolka_z_offset + 2 * konzolka_z - r);
Xk_small_all = xk_s(:);
Yk_small_all = yk_s(:);
Zk_small_all = zk_s(:);

relX_s = (Xk_small_all - sloup_x) ./ konzolka_x;
Zline_s = konzolka_z_offset + relX_s * konzolka_z;
mask_s = Zk_small_all >= Zline_s + r_small; 

Xk_small = Xk_small_all(mask_s);
Yk_small = Yk_small_all(mask_s);
Zk_small = Zk_small_all(mask_s);


idx_konzolka_small_valid = Xk_small > sloup_x;
Xk_small = Xk_small(idx_konzolka_small_valid);
Yk_small = Yk_small(idx_konzolka_small_valid);
Zk_small = Zk_small(idx_konzolka_small_valid);

%% === SPOJENÍ VŠECH BODŮ
X_all = [Xs_sloup; Xs_small_sloup; Xk; Xk_small];
Y_all = [Ys_sloup; Ys_small_sloup; Yk; Yk_small];
Z_all = [Zs_sloup; Zs_small_sloup; Zk; Zk_small];

R_all = [ ...
    r       * ones(size(Xs_sloup)); ...
    r_small * ones(size(Xs_small_sloup)); ...
    r       * ones(size(Xk)); ...
    r_small * ones(size(Xk_small)) ];

problem.geometrySet.node=[X_all'; Y_all'; Z_all'];
problem.geometrySet.radius=R_all';
problem.geometrySet.prop=[ones(size(Xs_sloup))',ones(size(Xs_small_sloup))'*2,ones(size(Xk))',ones(size(Xk_small))'*2];
problem.geometrySet.propMat=ones(size(problem.geometrySet.prop));

%% Boundary conditions


dirichlets{1}.position = [find(and(problem.geometrySet.node(1,:)==0,problem.geometrySet.node(2,:)==0))];
dirichlets{1}.dofs=[1 1 1 1 1 1];

problem.dirichlets=Dirichlets(dirichlets,problem);

%% Force conditions

       problem.forceSetting.type = "time";
       problem.forceSetting.startTime = 0;
       problem.forceSetting.endTime = 0.1;
       problem.forceSetting.stopOnCrack = 0;
       problem.force{1}.value = -500000;
       problem.force{1}.direction = 3;

       problem.force{1}.position = [find( ...
            problem.geometrySet.node(1,:) >= 0.8 & problem.geometrySet.node(1,:) <= 1.0 & ...
            problem.geometrySet.node(2,:) >= 0.2 & problem.geometrySet.node(2,:) <= 0.4 & ...
            problem.geometrySet.node(3,:) == 2.8 )];


end


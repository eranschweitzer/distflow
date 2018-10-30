function lateral_plot(bus, branch, refn, refl, varargin)

if nargin < 3
    refn = [];
    refl = [];
else
    if nargin < 4
      refl = [];
    end
end
refflag = ~isempty(refn);
refflagl= ~isempty(refl);
diffplot = ismember('diffplot', varargin);
sflag = isfield(branch,'S');
%% setup
ephasing = {branch.phase}.';
nphasing = {bus.phase}.';
nb = length(nphasing);
f = [branch.f].';
t = [branch.t].';

%% find 3 phase branches
t3ph = zeros(length(t),3);
f3ph = zeros(length(f),3);

for k = 1:length(ephasing)
    t3ph(k,ephasing{k}) = t(k);
    f3ph(k,ephasing{k}) = f(k);
end

% mask is 1  in entr l, phi, if line l contains phase phi.
mask = t3ph & f3ph;
%% connection matrices per phase
T = {t(mask(:,1)), t(mask(:,2)), t(mask(:,3))};
F = {f(mask(:,1)), f(mask(:,2)), f(mask(:,3))};

% Incidence matrix
M = cellfun(@(x,y) sparse([1:length(x), 1:length(x)], [x, y], [ones(1,length(x)), -ones(1,length(x))], length(x), nb),...
    F, T, 'UniformOutput', false);

% laplacian
L = cellfun(@(x) x'*x, M, 'UniformOutput', false);

% node degree
d = cellfun(@(x) diag(x), L, 'UniformOutput', false); 

% adjacency
A = cellfun(@(x, y) diag(x) - y, d, L, 'UniformOutput', false);

%% from laterals
% laterals are defined by the number of leaf nodes.
% essentially a lateral is a path to a leaf.
laterals = cellfun(@(x) cell(sum(x==1)-1,1), d, 'UniformOutput', false);
if sflag
  laterals_s = laterals;
end

for phi = 1:3
    leaves = find(d{phi} == 1);
    leaves = leaves(leaves ~= 1); % remove sourse node
    for l = 1:length(leaves)
        p = dijkstra(A{phi}, 1, leaves(l)); % path to source
        % branch ids are equal to TO node ids - 1. Therefore remove source
        % and decrement count by 1
        pedge = p(2:end) - 1;
        vm = cellfun(@(x,y) x(y==phi), {bus(p).vm}, {bus(p).phase});
        if sflag
            S = cellfun(@(x,y) x(y==phi), {branch(pedge).S}, {branch(pedge).phase});
            if refflagl
              Sref = cellfun(@(x,y) x(y==phi), {refl(pedge).S}, {refl(pedge).phase});
              laterals_s{phi}{l} = [ensure_col_vect(pedge) ensure_col_vect(S) ensure_col_vect(Sref)];
            else
              laterals_s{phi}{l} = [ensure_col_vect(pedge) ensure_col_vect(S)];
            end
        end
        if ~refflag
            laterals{phi}{l} = [ensure_col_vect(p) ensure_col_vect(vm)];
        else
            vmref = cellfun(@(x,y) x(y==phi), {refn(p).vm}, {refn(p).phase});
            laterals{phi}{l} = [ensure_col_vect(p) ensure_col_vect(vm) ensure_col_vect(vmref)];
        end
    end
end

%%
cmap = colormap('lines');
for phi = 1:3
    subplot(3,1,phi)
    for l = 1:length(laterals{phi})
        if refflag && diffplot
            plot(laterals{phi}{l}(:,1),laterals{phi}{l}(:,3) - laterals{phi}{l}(:,2),'Color', cmap(1,:), 'Marker','.', 'MarkerSize', 10);
        else
            if refflag
                plot(laterals{phi}{l}(:,1), laterals{phi}{l}(:,3),'k--', 'Marker','.', 'MarkerSize', 10);
                hold on;
            end
            plot(laterals{phi}{l}(:,1), laterals{phi}{l}(:,2),'Color', cmap(1,:), 'Marker','.', 'MarkerSize', 10);
        end
        if l == 1
            hold on;
        end
    end
    if refflag && ~diffplot
        legend('opendss', 'distflow')
    end
    title(sprintf('Phase %d', phi))
    xlabel('Node Number')
    if diffplot
        ylabel('Vm_{OpenDss} - Vm_{DistFlow} [p.u]')
    else
        ylabel('Vm [p.u]')
    end
    set(gca, 'FontSize', 16)
    xlim([1, nb])
end

if sflag
  figure;
  for phi = 1:3
    subplot(3,1,phi)
    for l = 1:length(laterals_s{phi})
        if refflagl && diffplot
            yyaxis left;
            plot(laterals_s{phi}{l}(:,1),real(laterals_s{phi}{l}(:,3) - laterals_s{phi}{l}(:,2)),'-','Color', cmap(1,:), 'Marker','.', 'MarkerSize', 10);
            yyaxis right;
            plot(laterals_s{phi}{l}(:,1),imag(laterals_s{phi}{l}(:,3) - laterals_s{phi}{l}(:,2)),'-','Color', cmap(2,:), 'Marker','.', 'MarkerSize', 10);
        else
            if refflagl
                yyaxis left;
                plot(laterals_s{phi}{l}(:,1), real(laterals_s{phi}{l}(:,3)),'k--', 'Marker','.', 'MarkerSize', 10);
                yyaxis right;
                plot(laterals_s{phi}{l}(:,1), imag(laterals_s{phi}{l}(:,3)),'k:', 'Marker','.', 'MarkerSize', 10);
                hold on;
            end
            yyaxis left;
            plot(laterals_s{phi}{l}(:,1), real(laterals_s{phi}{l}(:,2)),'-','Color', cmap(1,:), 'Marker','.', 'MarkerSize', 10);
            yyaxis right;
            plot(laterals_s{phi}{l}(:,1), imag(laterals_s{phi}{l}(:,2)),'-','Color', cmap(2,:), 'Marker','.', 'MarkerSize', 10);
        end
        if l == 1
            hold on;
        end
    end
    if refflag && ~diffplot
        legend('opendss', 'distflow')
    end
    yyaxis left;
    title(sprintf('Phase %d', phi))
    xlabel('Branch Number')
    if diffplot
        ylabel('P_{OpenDss} - P_{DistFlow} [p.u]')
        yyaxis right;
        ylabel('Q_{OpenDss} - Q_{DistFlow} [p.u]')
    else
        ylabel('P [p.u]')
        yyaxis right;
        ylabel('Q [p.u]')
    end
    set(gca, 'FontSize', 16)
    xlim([1, nb-1])
  end
end


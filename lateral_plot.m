function lateral_plot(bus, branch, ref, varargin)

if nargin < 3
    refflag = 0;
else
    refflag = 1;
end

diffplot = ismember('diffplot', varargin);

ephasing = {branch.phase}.';
nphasing = {bus.phase}.';
nb = length(nphasing);
f = [branch.f].';
t = [branch.t].';

t3ph = zeros(length(t),3);
f3ph = zeros(length(f),3);

for k = 1:length(ephasing)
    t3ph(k,ephasing{k}) = t(k);
    f3ph(k,ephasing{k}) = f(k);
end

mask = t3ph & f3ph;

T = {t(mask(:,1)), t(mask(:,2)), t(mask(:,3))};
F = {f(mask(:,1)), f(mask(:,2)), f(mask(:,3))};

M = cellfun(@(x,y) sparse([1:length(x), 1:length(x)], [x, y], [ones(1,length(x)), -ones(1,length(x))], length(x), nb),...
    F, T, 'UniformOutput', false);

L = cellfun(@(x) x'*x, M, 'UniformOutput', false);

d = cellfun(@(x) diag(x), L, 'UniformOutput', false);

A = cellfun(@(x, y) diag(x) - y, d, L, 'UniformOutput', false);

laterals = cellfun(@(x) cell(sum(x==1)-1,1), d, 'UniformOutput', false);

for phi = 1:3
    leaves = find(d{phi} == 1);
    leaves = leaves(leaves ~= 1); % remove sourse node
    for l = 1:length(leaves)
        p = dijkstra(A{phi}, 1, leaves(l));
        vm = cellfun(@(x,y) x(y==phi), {bus(p).vm}, {bus(p).phase});
        if ~refflag
            laterals{phi}{l} = [ensure_col_vect(p) ensure_col_vect(vm)];
        else
            vmref = cellfun(@(x,y) x(y==phi), {ref(p).vm}, {ref(p).phase});
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




function singlephase_lateral_plot(mpc,labels)

ctmp = colormap('lines');
cmap = [[0 0 0]; ctmp([1,2,4,5,6,7],:)];
line_styles = {'-','--', ':', '-.'};
define_constants;
if iscell(mpc)
  ncases = length(mpc); 
else
  ncases = 1;
  mpc = {mpc};
end
if nargin < 2
  labels = cellfun(@(x) sprintf('case%d',x), num2cell(1:ncases),'UniformOutput', false);
end
laterals = singlephase_laterals(mpc{1});

%% collect values
vm = cell(length(laterals.nid),1);
pf = cell(length(laterals.bid),1);
qf = cell(length(laterals.bid),1);
for l = 1:length(vm)
  pn = laterals.nid{l};
  pb = laterals.bid{l};
  vm{l} = cellfun(@(x) x.bus(pn,VM), mpc,'UniformOutput',false);
  pf{l} = cellfun(@(x) x.branch(pb,PF), mpc,'UniformOutput',false);
  qf{l} = cellfun(@(x) x.branch(pb,QF), mpc,'UniformOutput',false);
end

%%
for l = 1:length(vm)
  %%%%% voltage magnitude
  subplot(2,3,[1,4])
  for k = 1:ncases
    plot(laterals.nid{l},vm{l}{k},'Color',cmap(mod(k-1,8)+1,:), ...
      'LineStyle', line_styles{mod(k-1,4)+1}, 'Marker', '.', 'MarkerSize', 14);
    if (l == 1) && (k == 1)
      hold on;
    end
  end
  %%%%% real power
  subplot(2,3,[2,3])
  for k = 1:ncases
    plot(laterals.bid{l},pf{l}{k},'Color',cmap(mod(k-1,8)+1,:), ...
      'LineStyle', line_styles{mod(k-1,4)+1}, 'Marker', '.', 'MarkerSize', 14);
    if (l == 1) && (k == 1)
      hold on;
    end
  end
  %%%%% reactive power
  subplot(2,3,[5,6])
  for k = 1:ncases
    plot(laterals.bid{l}, qf{l}{k},'Color',cmap(mod(k-1,8)+1,:), ...
      'LineStyle', line_styles{mod(k-1,4)+1}, 'Marker', '.', 'MarkerSize', 14);
    if (l == 1) && (k == 1)
      hold on;
    end
  end
end

subplot(2,3,[1,4])
legend(labels);
xlabel('Node Number')
ylabel('Vm [p.u.]')
set(findall(gca, 'Type', 'Line'),'LineWidth',2);
set(gca,'FontSize', 16);

subplot(2,3,[2,3])
legend(labels);
xlabel('Branch Number')
ylabel('Pf [MW]')
set(findall(gca, 'Type', 'Line'),'LineWidth',2);
set(gca,'FontSize', 16);

subplot(2,3,[5,6])
legend(labels);
xlabel('Branch Number')
ylabel('Qf [MVAr]')
set(findall(gca, 'Type', 'Line'),'LineWidth',2);
set(gca,'FontSize', 16);

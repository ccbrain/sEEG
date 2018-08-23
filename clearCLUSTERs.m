preCluster=who('CLUSTER*');
noCluster=who('NOCLUSTER*');

for i=1:length(preCluster)
    clear(preCluster{i});
end

for i=1:length(noCluster)
    clear(noCluster{i});
end
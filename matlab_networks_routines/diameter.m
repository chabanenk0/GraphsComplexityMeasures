% The longest shortest path between any two nodes nodes in the network
% INPUTS: adjacency matrix, adj
% OUTPUTS: network diameter, diam
% Other routines used: simple_dijkstra.m
% GB, Last updated: June 8, 2010

function [diam,dij] = diameter(adj)

%from ave_path_length.m
%for i=1:n; dij=[dij; simple_dijkstra(adj,i) ]; end
%from smooth_diameter.m
%for i=1:n; dij=[dij; simple_dijkstra(adj,i)]; end
%from vertex_eccentrality
%for s=1:n; ec(s)=max( simple_dijkstra(adj,s) ); end
%from closeness
%for i=1:length(adj); C(i)=1/sum( simple_dijkstra(adj,i) ); end

diam=0;
dij=[];
for i=1:size(adj,1)
    d=simple_dijkstra(adj,i);
    dij=[dij;d];
    diam = max([max(d),diam]);
end
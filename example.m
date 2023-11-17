% This example is from http://www.pnas.org/cgi/doi/10.1073/pnas.0305937101
% with 11 Nodes, which are Cln3; MBF; SBF; Cln1,2; Cdh1; Swi5; Cdc20; Clb5,6; Sic1; Clb1,2; Mcm1

m = 11;   % number of nodes, could be set to less than or equal to 11
n = 4;    % number of nodes selected for states, and others are considered as controls

if m > 11
	error('The number of nodes cannot be large than 11.');
end
if n >= m
	error('The number of selected nodes must be small than the number of all nodes.');
end

% the adjacent matrix of the Boolean network
bn_adj = [
    0	0	0	0	0	0	0	0	0	0	0;
    1	0	0	0	0	0	0	0	0	-1	0;
    1	0	0	0	0	0	0	0	0	-1	0;
    0	0	1	0	0	0	0	0	0	0	0;
    0	0	0	-1	0	0	1	-1	0	-1	0;
    0	0	0	0	0	0	1	0	0	-1	1;
    0	0	0	0	0	0	0	0	0	1	1;
    0	1	0	0	0	0	-1	0	-1	0	0;
    0	0	0	-1	0	1	1	-1	0	-1	0;
    0	0	0	0	-1	0	-1	1	-1	0	1;
    0	0	0	0	0	0	0	1	0	1	0];
self_degradation_nodes = [1 4 6 7 11]; % the node index set with self degradation
if m < 11
	self_degradation_nodes = setdiff(self_degradation_nodes,(m+1):11);
end
d = 2^m;
n2 = 2^n;
d2 = 2^(m-n);
bn_adj = bn_adj(1:m,1:m);
theta = zeros(m,1); % threshhold for each node
%T = zeros(m,d); % the transition matrix of the Boolean network
TT = zeros(n,d); % the transition matrix of the logical equations
bn_adj_new = zeros(m,m);

P = wij(n2,d2); % the swap matrix

N = nchoosek(1:m, n); % choose n state nodes from all m nodes
for k = 1:size(N,1) % loop for selected state nodes
	nodes_selected = N(k,:);
	nodes_control = setdiff(1:m,nodes_selected);
	% put all the selected state nodes before the rest control nodes
	nodes_neworder = [nodes_selected, nodes_control];
	bn_adj_new = bn_adj(nodes_neworder, nodes_neworder);
	self_degradation_nodes_new = zeros(size(self_degradation_nodes));
	for i = 1:length(self_degradation_nodes)
		self_degradation_nodes_new(i) = find(nodes_neworder == self_degradation_nodes(i));
	end
	F_OL = zeros(n2,d);
	for j = 1:d
		x = dec2bin(j-1,m) - '0';
		x1 = bn_adj_new*x';
		g_theta = find(x1>theta);
		l_theta = find(x1<theta);
		e_theta = find(x1==theta);
		x1(g_theta) = 1;
		x1(l_theta) = 0;
		x1(e_theta) = x(e_theta);
		%T(g_theta,j) = 1;
		%T(l_theta,j) = 0;
		%T(e_theta,j) = x(e_theta);
		for sdn = self_degradation_nodes_new
			if sdn<=m && x(sdn)==1 && any(find(e_theta==sdn))
				x1(sdn) = 0;
				%T(sdn,j) = 0;
			end
		end
		for i = 1:n % nodes_selected
			if x(i) == x1(i)
				TT(i,j) = 0;
			else
				TT(i,j) = 1;
			end
		end
		G = s2v(TT(1,j));
		for i = 2:n
			G = kron(G,s2v(TT(i,j)));
		end
		F_OL(:,j) = G;
	end
	F = F_OL*P;
	for s = 1:n2
		F_s = F(:,s:n2:end); % Get the block (F_OL*P)_s
		if rank(F_s) < n2 % check the rank condition
			%disp(nodes_selected);
			break
		end
		if s == n2
			fprintf(['When nodes [',strcat(repmat('%i ',[1,m-n])), '] are selected as controls, the system is feedback shapable!\n'],nodes_control);
		end
	end
end

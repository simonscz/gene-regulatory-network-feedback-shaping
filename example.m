m = 11;   % number of nodes, could be set to less than or equal to 11
n = 4;    % number of nodes selected for states, and others are considered as controls

d = 2^m; % number of states
% the adjacent matrix of the Boolean network
bn_adj = [
    0	-1	0	0	0	0	0	0	0	0	-1;
    0	0	0	0	0	-1	0	0	0	0	0;
    1	0	0	0	0	0	0	0	0	0	0;
    0	0	0	0	0	1	0	0	1	1	0;
    0	0	1	0	0	0	1	0	0	0	0;
    0	-1	0	1	0	0	-1	0	1	-1	0;
    0	0	0	0	0	0	0	1	0	0	0;
    0	1	0	1	0	1	0	0	0	0	-1;
    0	1	0	0	0	-1	0	-1	0	1	1;
    0	0	0	0	0	0	0	0	0	0	1;
    0	0	0	0	0	0	0	1	0	0	0];
bn_adj = bn_adj(1:m,1:m);
theta = zeros(m,1); % threshhold for each node
T = zeros(m,d); % initialize the transition matrix of the Boolean network
TT = zeros(n,d); % initialize the transition matrix of the logical equations
for i1 = 1:8
    for i2 = i1+1 : 9
        for i3 = i2+1: 10
            for i4 = i3+1:11 % loop for selected state nodes
                    F_OL = zeros(2^n,2^m);
                    for j = 1:d
		                x = dec2bin(j-1,m) - '0';
		                x1 = bn_adj*x';
		                g_theta = find(x1>theta);
		                l_theta = find(x1<theta);
		                e_theta = find(x1==theta);
		                x1(g_theta) = 1;
		                x1(l_theta) = 0;
		                x1(e_theta) = x(e_theta);
		                T(g_theta,j) = 1;
		                T(l_theta,j) = 0;
		                T(e_theta,j) = x(e_theta);
		                k = 1;
		                for i = [i1 i2 i3 i4]
			                if x(i) == x1(i)
				                TT(k,j) = 0;
			                else
				                TT(k,j) = 1;
			                end
			                k = k + 1;
                        end
                        G=s2v(TT(1,j));
                        for i = 2:n
                            G = kron(G,s2v(TT(i,j)));
                        end
                        for i = 1:2^n
                            F_OL(i,j) = G(i);
                        end
                    end
	                P = wij(2^n,2^(m-n));
	                F = F_OL*P;
                    for s = 1:2^n
                        FF = zeros(2^n,2^(m-n));

                        %分块，FF为论文中(F_OL*P)_s
                        for i = 1:2^n
                            for j = 1:2^(m-n)
                                FF(i,j) = F(i,s+(j-1)*2^n);
                            end
                        end
                        if rank(FF)<2^n
                            break 
                        end
                        if s == 2^n
                            ii = [i1,i2,i3,i4];
                            fprintf('状态量为%i,%i,%i,%i时，是feedback shapable的\n',ii);
                        end
                    end
            end
        end
    end
end
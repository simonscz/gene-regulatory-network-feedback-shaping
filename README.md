# Feedback Shaping of Gene Regulatory Network 

## Requirements

The codes are implemented in Matlab/Octave, and no special toolboxes are needed.

## Instructions

The example is a 11-node cell-cycle regulatory network of budding yeast from a PNAS paper ([http://www.pnas.org/cgi/doi/10.1073/pnas.0305937101]). 

### Network Dynamics

The dynamics of the network is defined over a directed graph, where each arc $(j, i)$ is associated with a weight $`a_{ij} \in\{-1,1\}`$. Then the evolution of the Boolean states is governed by

$$ 
\begin{array}{l}
x_i(t+1)=f_i(\sum_{j\in{N_i}}a_{ij}x_j(t)),i=2,3,5,8,9,10;\\
x_i(t+1)=f_i^*(\sum_{j\in{N_i}}a_{ij}x_j(t)),i=1,4,6,7,11,
\end{array}
$$

where

$$
f_i(a)=\begin{cases}
1, & \text{if}~ a > 0;\\
0, & \text{if}~ a < 0;\\
x_i, & \text{if}~ a = 0,
\end{cases}
$$

and

$$
f_i^*(a)=\begin{cases}
1, & \text{if}~ a > 0;\\
0, & \text{if}~ a \le 0.
\end{cases}
$$

### Control Nodes Selection

We select some nodes out of the eleven as control nodes. When a node $k$ is selected as a control node, we assume its state will be manipulated by external factors, and in this case, $x_k(t)$ becomes a control input that can take any binary value at any time. Then we verify whether the resulting controlled cell-cycle is feedback shapable by checking the rank condition.

### The Codes

example.m is the main procudre which will show what nodes are selected as controls will make the system feedback shapable.

s2v.m converts a Boolean scalar to its vector form, precisely, which maps 0 to [1,0]' and 1 to [0,1]'.
  
wij.m produces swap matrix, which is taken from the STP toolbox for Matlab/Octave.
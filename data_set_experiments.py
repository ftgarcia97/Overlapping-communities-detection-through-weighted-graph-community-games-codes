from exact_formulations import F_ShJK,F_ShMod,find_stable_partition
from weights_generator import weights_generator,expected_weights_generator,approx_expected_weights_generator
from LSE_heuristic import LSE_heuristic
from Overlapping_LFR_benchmark_generator import Overlapping_LFR_benchmark_generator
import omega_index
from NMI_overlapping import onmi 

comu1=[{3,4,5},{1,2}]
comu2=[{1,2},{3,4,5}]
print(onmi([{3,4,5},{1,2}],[{1,2},{3,4,5}]))
nc=11
p=2


edges=open('data/metabolic.txt').read()
edges= [item.split() for item in edges.split('\n')[:-1]]
edges=[tuple([int(i[0]),int(i[1])]) for i in edges if(i[0]!=i[1])]
nodes=set([i for (i,j) in edges]+[j for (i,j) in edges])
print(len(nodes))
# W=weights_generator(nodes, edges)
# W_e=expected_weights_generator(nodes, edges)
# W_e=approx_expected_weights_generator(nodes, edges)

# W_star={}
# for i in nodes:
#     for j in nodes:
#         if(i!=j):
#             W_star[i,j]=W[i,j]-W_e[i,j]
#         else:
#             W_star[i,j]=0
            
banned_partitions=[]

# t_max=10
# obj_val={}
# partition={}
# bridge_nodes={}
# for t in range(t_max):
#     obj_val[t],partition[t],bridge_nodes[t],initial_sol=LSE_heuristic(nodes,edges,W_star,nc,p,banned_partitions)
#     print(obj_val[t],partition[t],bridge_nodes[t])
#     banned_partitions.append(initial_sol)
# arg_max=max(obj_val, key=obj_val.get)
# print(obj_val[arg_max],partition[arg_max],bridge_nodes[arg_max])


nc=25
N = 500
gamma = 2
beta = 1
mu = 0.2
mu_puente=0.7
k_min=20
k_max=30
s_min=20
s_max=30
p=3
N_o=20
edges,nc,S,bridge_nodes=Overlapping_LFR_benchmark_generator(N,N_o,k_min,k_max,s_min,s_max,gamma,beta,mu,mu_puente,p,nc)
print(omega_index.Omega(S,S).omega_score)
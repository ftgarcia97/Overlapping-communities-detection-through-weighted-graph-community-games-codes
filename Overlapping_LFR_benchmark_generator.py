# -*- coding: utf-8 -*-
"""
Created on Wed Feb  2 09:49:38 2022

@author: 34625
"""
import random
import networkx as nx
import math

##  Overlapping LFR benchmark generator algorithm
def Overlapping_LFR_benchmark_generator(N,N_o,k_min,k_max,s_min,s_max,gamma,beta,mu,mu_o,p,nc):
    degrees={}
    nodes=range(1,N+1)
    bridge_nodes=[]
    while(len(bridge_nodes)<N_o):
        bridge=random.choice(nodes)
        if(bridge not in bridge_nodes):
            bridge_nodes.append(bridge)
    infeasible=True
    while(infeasible):
        for i in nodes:
            degrees[i]=random.choices(range(k_min,k_max+1),weights=[k**(-gamma) for k in range(k_min,k_max+1)])[0]
        if(sum([degrees[i] for i in nodes]) % 2 ==0):
            infeasible=False
        for i in bridge_nodes:
            if(math.ceil(degrees[i]*(1-mu_o))*p>degrees[i]):
                infeasible=True
    community_sizes=[]
    sum_sizes=0
    
    while(sum_sizes<N+N_o*(p-1)):
        if(N+N_o*(p-1)-sum_sizes<s_min):
            sum_sizes=0
            community_sizes=[]
        size=random.choices(range(s_min,s_max+1),weights=[k**(-beta) for k in range(s_min,s_max+1)])[0]
        if(size>N+N_o*(p-1)-sum_sizes):
            size=N+N_o*(p-1)-sum_sizes
        community_sizes.append(size)
        sum_sizes=sum_sizes+size
        if(len(community_sizes)>nc):
            sum_sizes=0
            community_sizes=[]
    nc=len(community_sizes)
    
    S={}
    for s in range(1,nc+1):
        S[s]=set()
    C={}
    for i in bridge_nodes:
        C[i]=set()
    candidate_nodes=set(nodes)
    while(candidate_nodes!=set()):
        node=random.choice(list(candidate_nodes))
        community=random.choice(range(1,nc+1))
        if(node in bridge_nodes):
            if(node not in S[community]):
                S[community].add(node)
                C[node].add(community)
                if(len(C[node])==p):
                    candidate_nodes=candidate_nodes-{node}
                    
        else:
            S[community].add(node)
            candidate_nodes=candidate_nodes-{node}
        if(len(S[community])>community_sizes[community-1]):
            out_node=random.choice(list(S[community]))
            S[community]=S[community]-{out_node}
            candidate_nodes.add(out_node)
            if(out_node in bridge_nodes):
                C[out_node]=C[out_node]-{community}
    internal_degrees={}
    external_degrees={}
    for i in nodes:
        internal_degrees[i]=0
    for s in range(1,nc+1):
        for i in S[s]:
            if(i in bridge_nodes):
                internal_degrees[i,s]=math.ceil(degrees[i]*(1-mu_o))
            else:
                internal_degrees[i,s]=math.ceil(degrees[i]*(1-mu))
            internal_degrees[i]=internal_degrees[i]+internal_degrees[i,s]
            
    for s in range(1,nc+1):
        degrees_sum=0
        for i in S[s]:
            degrees_sum=degrees_sum+internal_degrees[i,s]
            
        if(degrees_sum % 2 !=0):
            nodes_with_internal_space=[i for i in S[s] if(internal_degrees[i]<degrees[i])]
            if(len(nodes_with_internal_space)>0):
                node=random.choice(nodes_with_internal_space)
                internal_degrees[node,s]=internal_degrees[node,s]+1
                internal_degrees[node]=internal_degrees[node]+1
                degrees_sum=degrees_sum+1
            else:
                node=random.choice(list(S[s]))
                internal_degrees[node,s]=internal_degrees[node,s]-1
                internal_degrees[node]=internal_degrees[node]-1
                degrees_sum=degrees_sum-1
            
            
    edges=[]
    for i in nodes:
        external_degrees[i]=degrees[i]-internal_degrees[i]
    for s in range(1,nc+1):
        community=list(S[s])
        random_graph=nx.configuration_model([internal_degrees[i,s] for i in community])
        edges=edges+[(community[i],community[j]) for (i,j) in random_graph.edges()]
    random_graph=nx.configuration_model([external_degrees[i] for i in nodes])
    edges=edges+[(nodes[i],nodes[j]) for (i,j) in random_graph.edges()]
    for i in nodes:
        if(degrees[i]!=len([1 for (j1,j2) in edges if(i==j1)])+len([1 for (j1,j2) in edges if(i==j2)])):
            print("ERROR",i)
    for s in range(1,nc+1):
        for i in list(S[s]):
            if(internal_degrees[i,s]>len([1 for (j1,j2) in edges if(i==j1 and j2 in list(S[s]))])+len([1 for (j1,j2) in edges if(i==j2 and j1 in list(S[s]))])):
                print("ERROR",i,s)
    return edges,nc,S,bridge_nodes
    

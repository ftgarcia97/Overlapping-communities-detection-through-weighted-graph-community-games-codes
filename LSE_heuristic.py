import random
from exact_formulations import find_stable_partition,F_ShMod,find_stable_partition3

## Local stability exploration algorithm:
    ## W: Weights - Expected Weights
    ## nodes: Set of nodes
    ## nc: Maximum number of communities
    ## p: Maximum number of communities to which a node can belong to
def LSE_heuristic(nodes,edges,W,nc,p,banned_partitions):
    adjacent_degree={}
    for i in nodes:
        adjacent_degree[i]=len([(j,k) for (j,k) in edges if(j==i or k==i)]) 
    stability=False
    while(stability==False):
        ## Random initial partition
        # initial_partition=random_partition(nodes,nc)
        initial_sol,initial_partition=find_stable_partition(W,nodes,nc,p,banned_partitions)
        # initial_partition=find_stable_partition3(W,nodes,nc,p,banned_partitions)
        # obj_val1,initial_partition=F_ShMod(W,nodes,nc,p)
        ## S[k]: Set of nodes that belong to community k
        S={}
        
        ## C[i]: Index Set of the communities to which i belongs
        C={}
        for i in nodes:
            C[i]=[]
        for k in range(1,nc+1):
            S[k]=initial_partition[k-1]
            for i in S[k]:
                C[i].append(k)
        
        ## obj_val: Objective function value
        obj_val=obj_func(nodes, W, S, nc)
        ## Initialize the maximum delta
        delta_max=1
        ## Initialize FeasibleMoves set
        moves=[]
        stability=is_stable(nodes,W,S,nc)
        # print(stability,initial_partition)
        while(delta_max>1.e-10 and stability==False):
            
            ## The initial partition is not stable so we have to find first an stable one
            stable_moves=[]
            delta={}
            
            no_stable_communities=[]
            for s in range(1,nc+1):
                for i in S[s]:
                    if(sum([W[i,j] for j in S[s] if(j!=i)])-sum([W[i,j] for j in nodes if(j not in S[s])])<0):
                        no_stable_communities.append(s)
                        break
            
            ## Add movement
            d=1
            for i in nodes:
                for s in (set(no_stable_communities)-set(C[i])):
                    if(p>=(len(C[i])+1)):
                        S[s]=S[s]+[i]
                        if(no_stable_communities==[i]):
                            if(is_stable(nodes,W,S,nc)):
                                stable_moves.append((i,s,d))
                        S[s]=list(set(S[s])-{i})
                        delta[i,s,d]=sum([W[i,j] for j in S[s] if(len(set(C[i]) & set(C[j]))==0)])
                    else:
                        delta[i,s,d]=0
            
            ## Remove movement
            d=2
            for i in nodes:
                for s in (set(C[i]) & set(no_stable_communities)):
                    if(len(C[i])>1):
                        S[s]=list(set(S[s])-{i})
                        if(no_stable_communities==[i]):
                            if(is_stable(nodes,W,S,nc)):
                                stable_moves.append((i,s,d))
                        S[s]=S[s]+[i]
                        delta[i,s,d]=-sum([W[i,j] for j in S[s] if(len(set(C[i]) & set(C[j]))==1 and j!=i)])
                    else:
                        delta[i,s,d]=0
            
            ## Swap movement
            # d=3
            # for i in nodes:
            #     for j in nodes:
            #         if(i<j):
            #             for s in (set(C[i])-set(C[j])) & set(no_stable_communities):
            #                 for r in (set(C[j])-set(C[i])) & set(no_stable_communities):
            #                     S[s]=list(set(S[s]+[j])-{i})
            #                     S[r]=list(set(S[r]+[i])-{j})
            #                     if(no_stable_communities==[i] or no_stable_communities==[j] or no_stable_communities==[i,j]):
            #                         if(is_stable(nodes,W,S,nc)):
            #                             stable_moves.append((i,s,j,r,d))
            #                     S[s]=list(set(S[s]+[i])-{j})
            #                     S[r]=list(set(S[r]+[j])-{i})
            #                     delta[i,s,j,r,d]=sum([W[i,k] for k in S[r] if(len(set(C[i]) & set(C[k]))==0 and k!=j)])+sum([W[j,k] for k in S[s] if(len(set(C[j]) & set(C[k]))==0 and k!=i)])-sum([W[i,k] for k in S[s] if(len(set(C[i]) & set(C[k]))==1 and k!=i and k not in S[r])])-sum([W[j,k]for k in S[r] if(len(set(C[k]) & set(C[j]))==1 and k!=j and k not in S[s])])   
            if(len(stable_moves)>0):
                delta_stable={}
                for t in stable_moves:
                    if(len(t)==3):
                        i,s,d=t
                        delta_stable[i,s,d]=delta[i,s,d]
                    if(len(t)==5):
                        i,s,j,r,d=t
                        delta_stable[i,s,j,r,d]=delta[i,s,j,r,d]
                arg_max=max(delta_stable, key=delta_stable.get)
                if(len(arg_max)==3):
                    node_max,com_max,d_max=arg_max
                    delta_max=delta_stable[node_max,com_max,d_max]
                if(len(arg_max)==5):
                    node1_max,com1_max,node2_max,com2_max,d_max=arg_max
                    delta_max=delta_stable[node1_max,com1_max,node2_max,com2_max,d_max]
                stability=True
            else:
                arg_max=max(delta, key=delta.get)
                if(len(arg_max)==3):
                    node_max,com_max,d_max=arg_max
                    delta_max=delta[node_max,com_max,d_max]
                if(len(arg_max)==5):
                    node1_max,com1_max,node2_max,com2_max,d_max=arg_max
                    delta_max=delta[node1_max,com1_max,node2_max,com2_max,d_max]
                          
            if(delta_max>1.e-10 or stability):
                if(d_max==1):
                    S[com_max]=S[com_max]+[node_max]
                    moves.append((node_max,com_max,d_max,delta_max))
                    obj_val=obj_val+delta_max
                    C[node_max]=C[node_max]+[com_max]                        
                if(d_max==2):
                    S[com_max]=list(set(S[com_max])-{node_max})
                    moves.append((node_max,com_max,d_max,delta_max))
                    obj_val=obj_val+delta_max
                    C[node_max]=list(set(C[node_max])-{com_max})
                if(d_max==3):
                    S[com1_max]=list(set(S[com1_max]+[node2_max])-{node1_max})
                    S[com2_max]=list(set(S[com2_max]+[node1_max])-{node2_max})
                    moves.append((node1_max,com1_max,node2_max,com2_max,d_max,delta_max))
                    obj_val=obj_val+delta_max
                    C[node1_max]=list(set(C[node1_max]+[com2_max])-{com1_max})
                    C[node2_max]=list(set(C[node2_max]+[com1_max])-{com2_max})
            stability=is_stable(nodes,W,S,nc)
                    
    
    ## Now we have an initial stable partition
    
    ## Initialize the maximum delta again
    delta_max=1
    
    while(delta_max>1.e-10):
        delta={}
        print(moves)
        ## Add movement
        d=1
        for i in nodes:
            for s in (set(range(1,nc+1))-set(C[i])):
                if(p>=(len(C[i])+1)):
                    S[s]=S[s]+[i]
                    move_stability=is_stable_partial(nodes,W,S,nc,s)
                    S[s]=list(set(S[s])-{i})
                    if(move_stability):
                        delta[i,s,d]=sum([W[i,j] for j in S[s] if(len(set(C[i]) & set(C[j]))==0)])
                    else:
                        delta[i,s,d]=0
                else:
                    delta[i,s,d]=0
        
        ## Remove movement
        d=2
        for i in nodes:
            for s in C[i]:
                if(len(C[i])>1):
                    S[s]=list(set(S[s])-{i})
                    move_stability=is_stable_partial(nodes,W,S,nc,s)
                    S[s]=S[s]+[i]
                    if(move_stability):
                        delta[i,s,d]=-sum([W[i,j] for j in S[s] if(len(set(C[i]) & set(C[j]))==1 and j!=i)])
                    else:
                        delta[i,s,d]=0
                else:
                    delta[i,s,d]=0
        
        ## Swap movement
        # d=3
        # for i in nodes:
        #     for j in nodes:
        #         if(i<j):
        #             for s in (set(C[i])-set(C[j])):
        #                 for r in (set(C[j])-set(C[i])):
        #                     S[s]=list(set(S[s]+[j])-{i})
        #                     S[r]=list(set(S[r]+[i])-{j})
        #                     move_stability=is_stable_partial(nodes,W,S,nc,s) and is_stable_partial(nodes,W,S,nc,r)
        #                     S[s]=list(set(S[s]+[i])-{j})
        #                     S[r]=list(set(S[r]+[j])-{i})
        #                     if(move_stability):
        #                         delta[i,s,j,r,d]=sum([W[i,k] for k in S[r] if(len(set(C[i]) & set(C[k]))==0 and k!=j)])+sum([W[j,k] for k in S[s] if(len(set(C[j]) & set(C[k]))==0 and k!=i)])-sum([W[i,k] for k in S[s] if(len(set(C[i]) & set(C[k]))==1 and k!=i and k not in S[r])])-sum([W[j,k]for k in S[r] if(len(set(C[k]) & set(C[j]))==1 and k!=j and k not in S[s])])   
        #                     else:
        #                         delta[i,s,j,r,d]=0


        arg_max=max(delta, key=delta.get)
        if(len(arg_max)==3):
            node_max,com_max,d_max=arg_max
            delta_max=delta[node_max,com_max,d_max]
        if(len(arg_max)==5):
            node1_max,com1_max,node2_max,com2_max,d_max=arg_max
            delta_max=delta[node1_max,com1_max,node2_max,com2_max,d_max]
                      
        if(delta_max>1.e-10):
            if(d_max==1):
                S[com_max]=S[com_max]+[node_max]
                moves.append((node_max,com_max,d_max,delta_max))
                obj_val=obj_val+delta_max
                C[node_max]=C[node_max]+[com_max]                        
            if(d_max==2):
                S[com_max]=list(set(S[com_max])-{node_max})
                moves.append((node_max,com_max,d_max,delta_max))
                obj_val=obj_val+delta_max
                C[node_max]=list(set(C[node_max])-{com_max})
            if(d_max==3):
                S[com1_max]=list(set(S[com1_max]+[node2_max])-{node1_max})
                S[com2_max]=list(set(S[com2_max]+[node1_max])-{node2_max})
                moves.append((node1_max,com1_max,node2_max,com2_max,d_max,delta_max))
                obj_val=obj_val+delta_max
                C[node1_max]=list(set(C[node1_max]+[com2_max])-{com1_max})
                C[node2_max]=list(set(C[node2_max]+[com1_max])-{com2_max})

    
    #plt.savefig('heuristica1_football_lambda01.jpg')
    partition=[]
    for s in range(1,nc+1):
        if(len(S[s])>0):
            partition.append(set(S[s]))
            
    if(abs(obj_val-obj_func(nodes,W,S,nc))>1.e-5):
        print("ERROR Obj Val")
        obj_val=0
    if(is_stable(nodes,W,S,nc)==False):
        print("ERROR Stability")
        obj_val=0
        
    bridge_nodes=[]
    for i in nodes:
        if(len(C[i])>=2):
            bridge_nodes.append(i)
    
    return obj_val,partition,bridge_nodes,initial_sol

def random_partition(iterable,k):
  results = [[] for i in range(k)]
  for value in iterable:
    x = random.randrange(k)
    results[x].append(value)
  return results

def obj_func(nodes,W,community,nc):
    obj_val=0
    for i in nodes:
        for j in nodes:
            if(i<j):
                k=1
                while(k<=nc):
                    if(i in community[k] and j in community[k]):
                        obj_val=obj_val+W[i,j]
                        break
                    k=k+1
    return obj_val

def is_stable(nodes,W,community,nc):
    for s in range(1,nc+1):
        for i in community[s]:
            if(sum([W[i,j] for j in community[s] if(j!=i)])-sum([W[i,j] for j in nodes if(j not in community[s])])<0):
                return False
    return True

def is_stable_partial(nodes,W,community,nc,s):
    for i in community[s]:
        if(sum([W[i,j] for j in community[s] if(j!=i)])-sum([W[i,j] for j in nodes if(j not in community[s])])<0):
            return False
    return True


from probabilities import prob_adj,prob_CN,prob_triangle

## Weights built by Jonnalagadda and Kuppusamy:
    ## nodes: Set of nodes
    ## edges: Set of edges
def weights_generator(nodes,edges):
    adjacent_degree={}
    for i in nodes:
        adjacent_degree[i]=len([(j,k) for (j,k) in edges if(j==i or k==i)]) 
          
    neighbors={}
    for i in nodes:
        neighbors[i]=[j for j in nodes if(j!=i and ((i,j) in edges or (j,i) in edges))]
    
    # Parameters P_ij y CN_ij
    P,CN=P_CN_generator(nodes,edges)
        
    # Weights W_ij
    W={}
    for i in nodes:
        for j in nodes:
            if(i==j or adjacent_degree[i]==0 or adjacent_degree[j]==0):
                W[i,j]=0
            else:
                if(i in neighbors[j]):
                    if(adjacent_degree[i]==1 or adjacent_degree[j]==1):
                        W[i,j]=P[i,j]
                    else:
                        W[i,j]=2*CN[i,j]+P[i,j]             
                else:
                    W[i,j]=(CN[i,j]-P[i,j])/4
                    
    # Return the new weights
    return W


## Exact Expected Weights:
    ## nodes: Set of nodes
    ## edges: Set of edges
def expected_weights_generator(nodes,edges):
    
    n_edges=len(edges)
    adjacent_degree={}
    for i in nodes:
        adjacent_degree[i]=len([(j,k) for (j,k) in edges if(j==i or k==i)]) 
          
    neighbors={}
    for i in nodes:
        neighbors[i]=[j for j in nodes if(j!=i and ((i,j) in edges or (j,i) in edges))]
    
    P,CN=P_CN_generator(nodes,edges)
    
    ## Calculate all the different probabilities  
    non_repeated_degrees={adjacent_degree[i] for i in nodes}
    non_repeated_prob_adj={}
    non_repeated_prob_CN={}
    non_repeated_prob_triangle={}
    for i in non_repeated_degrees:
        for j in non_repeated_degrees:
            if(i<=j):
                non_repeated_prob_adj[i,j]=prob_adj(i,j,n_edges)
                if(i<j):
                    non_repeated_prob_adj[j,i]=non_repeated_prob_adj[i,j]
    for i in non_repeated_degrees:
        for j in non_repeated_degrees:
            for k in non_repeated_degrees:
                if(i<=j):
                    non_repeated_prob_CN[i,j,k]=prob_CN(i,j,k,n_edges)
                    if(i<j):
                        non_repeated_prob_CN[j,i,k]=non_repeated_prob_CN[i,j,k]
                if(i<=j<=k):
                    non_repeated_prob_triangle[i,j,k]=prob_triangle(i,j,k,n_edges)
                    non_repeated_prob_triangle[j,i,k]=non_repeated_prob_triangle[i,j,k]
                    non_repeated_prob_triangle[i,k,j]=non_repeated_prob_triangle[i,j,k]
                    non_repeated_prob_triangle[j,k,i]=non_repeated_prob_triangle[i,j,k]
                    non_repeated_prob_triangle[k,i,j]=non_repeated_prob_triangle[i,j,k]
                    non_repeated_prob_triangle[k,j,i]=non_repeated_prob_triangle[i,j,k]
    
    # Exact Expected Weights W^e_ij
    W_e={}
    for i in nodes:
        for j in nodes:
            if(i==j):
                W_e[i,j]=0
            if(i<j):
                if(adjacent_degree[i]==0 or adjacent_degree[j]==0):
                    W_e[i,j]=0
                else:
                    if(adjacent_degree[i]==1 or adjacent_degree[j]==1):
                        W_e[i,j]=non_repeated_prob_adj[adjacent_degree[i],adjacent_degree[j]]*P[i,j]+(sum([non_repeated_prob_CN[adjacent_degree[i],adjacent_degree[j],adjacent_degree[k]] for k in nodes if(k!=i and k!=j)]))*P[i,j]/4
                    else:
                        W_e[i,j]=3*non_repeated_prob_adj[adjacent_degree[i],adjacent_degree[j]]*P[i,j]+(sum([(non_repeated_prob_CN[adjacent_degree[i],adjacent_degree[j],adjacent_degree[k]]+7*non_repeated_prob_triangle[adjacent_degree[i],adjacent_degree[j],adjacent_degree[k]]) for k in nodes if(k!=i and k!=j)]))*P[i,j]/4
                W_e[j,i]=W_e[i,j]
    
    # Return exact expected weights
    return W_e


## Expected Weights approximation:
    ## nodes: Set of nodes
    ## edges: Set of edges
def approx_expected_weights_generator(nodes,edges):
    
    # W: Weights built by Jonnalagadda and Kuppusamy
    W=weights_generator(nodes,edges)
    weights_degree={}
    for i in nodes:
        weights_degree[i]=sum([W[i,k] for k in nodes if(k!=i)])
    total_weights=sum([weights_degree[i] for i in nodes])
    
    # Expected Weights approximation W^e_ij
    W_e={}
    for i in nodes:
        for j in nodes:
            W_e[i,j]=weights_degree[i]*weights_degree[j]/total_weights
    
    # Return expected weights approximation
    return W_e
                

## P and CN parameters generator:
    ## nodes: Set of nodes
    ## edges: Set of edges
def P_CN_generator(nodes,edges):
    adjacent_degree={}
    for i in nodes:
        adjacent_degree[i]=len([(j,k) for (j,k) in edges if(j==i or k==i)]) 
          
    neighbors={}
    for i in nodes:
        neighbors[i]=[j for j in nodes if(j!=i and ((i,j) in edges or (j,i) in edges))]
    
    # Parameters P_ij y CN_ij
    P={}
    CN={}
    for i in nodes:
        for j in nodes:
            if(i==j or adjacent_degree[i]==0 or adjacent_degree[j]==0):
                P[i,j]=0
                CN[i,j]=0
            else:
                P[i,j]=1/adjacent_degree[i]+1/adjacent_degree[j]
                CN[i,j]=(len([k for k in nodes if(k in neighbors[i] and k in neighbors[j] and k!=i and k!=j)])+1)*P[i,j]
                
    return P,CN

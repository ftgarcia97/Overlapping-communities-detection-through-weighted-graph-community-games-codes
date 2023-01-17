import gurobipy as gp
from gurobipy import GRB


## Formulation F_{Sh-JK}:
    ## W: Weights
    ## nodes: Set of nodes
    ## nc: Maximum number of communities
def F_ShJK(W,nodes,nc):
    m_ShJK = gp.Model("model F_ShJK")
    
    # Variables
    x={}
    z={}
    h={}

    for s in range(1,nc+1):
        for i in nodes:
            x[i,s]=m_ShJK.addVar(vtype=GRB.BINARY,lb=0,ub=1,name="x_"+str(i)+"_"+str(s))
            for r in range(s+1,nc+1):
                h[i,s,r]=m_ShJK.addVar(vtype=GRB.CONTINUOUS,lb=0,ub=1,name="h_"+str(i)+"_"+str(s)+"_"+str(r))
            for j in nodes:
                if(i<j):
                    z[i,j,s]=m_ShJK.addVar(vtype=GRB.CONTINUOUS,lb=0,ub=1,name="z_"+str(i)+"_"+str(j))

        
    
    # Objective function
    m_ShJK.setObjective(sum([sum([sum([W[i,j]*z[i,j,s] for i in nodes if(i<j)]) for j in nodes]) for s in range(1,nc+1)]), GRB.MAXIMIZE)
    
    # Constraints
    for i in nodes:
        m_ShJK.addConstr(sum([x[i,s] for s in range(1,nc+1)])>=1)
    for s in range(1,nc):
        m_ShJK.addConstr(sum([x[i,s] for i in nodes])>= sum([x[i,s+1] for i in nodes]))
    for s in range(1,nc+1):
        for i in nodes:
            m_ShJK.addConstr(sum([x[j,s]*W[i,j] for j in nodes if(j!=i)])>=x[i,s]*sum([W[i,j] for j in nodes if(j!=i)])/2)
            for j in nodes:
                if(i<j):
                    m_ShJK.addConstr(z[i,j,s]<=x[i,s])
                    m_ShJK.addConstr(z[i,j,s]<=x[j,s])
                    m_ShJK.addConstr(z[i,j,s]>=x[i,s]+x[j,s]-1)                   
        for r in range(s+1,nc+1):
            for i in nodes:
                m_ShJK.addConstr(h[i,s,r]<=x[i,r])
                m_ShJK.addConstr(h[i,s,r]<=1-x[i,s])
                m_ShJK.addConstr(x[i,r]-x[i,s]-h[i,s,r] <= 0)
                m_ShJK.addConstr(sum([h[j,s,r] for j in nodes])>=x[i,r])
        
    
    # Solve the model
    m_ShJK.setParam("LogToConsole",0)
    m_ShJK.optimize()
    
    partition=[]
    for s in range(1,nc+1): 
        community=[]
        for i in  nodes:
            if(x[i,s].x>=0.5):
                community.append(i)
        if(len(community)==0):
            break
        else:
            partition.append(community)
    # Return objective value and final partition
    return m_ShJK.objVal,partition



## Formulation F_{Sh-Mod}:
    ## W: Weights - Expected Weights
    ## nodes: Set of nodes
    ## nc: Maximum number of communities
    ## p: Maximum number of communities to which a node can belong to
def F_ShMod(W,nodes,nc,p):
    m_ShMod = gp.Model("model F_ShMod")
    
    # Variables
    x={}
    z={}
    y={}

    for s in range(1,nc+1):
        for i in nodes:
            x[i,s]=m_ShMod.addVar(vtype=GRB.BINARY,lb=0,ub=1,name="x_"+str(i)+"_"+str(s))
            for j in nodes:
                if(i<j):
                    z[i,j,s]=m_ShMod.addVar(vtype=GRB.CONTINUOUS,lb=0,ub=1,name="z_"+str(i)+"_"+str(j))
    
    for i in nodes:
        for j in nodes:
            if(i<j):
                y[i,j]=m_ShMod.addVar(vtype=GRB.CONTINUOUS,lb=0,ub=1, name="y"+str(i)+"_"+str(j))

        
    
    # Objective function
    m_ShMod.setObjective(sum([sum([W[i,j]*y[i,j] for i in nodes if(i<j)]) for j in nodes]), GRB.MAXIMIZE)
    
    # Constraints
    for i in nodes:
        m_ShMod.addConstr(sum([x[i,s] for s in range(1,nc+1)])>=1)
        m_ShMod.addConstr(sum([x[i,s] for s in range(1,nc+1)])<=p)
    for s in range(1,nc):
        m_ShMod.addConstr(sum([x[i,s] for i in nodes])>= sum([x[i,s+1] for i in nodes]))
    for s in range(1,nc+1):
        for i in nodes:
            m_ShMod.addConstr(sum([x[j,s]*W[i,j] for j in nodes if(j!=i)])>=x[i,s]*sum([W[i,j] for j in nodes if(j!=i)])/2+(1-x[i,s])*sum([W[i,j] for j in nodes if(W[i,j]<0 and j!=i)]))
            for j in nodes:
                if(i<j):
                    m_ShMod.addConstr(z[i,j,s]<=x[i,s])
                    m_ShMod.addConstr(z[i,j,s]<=x[j,s])
                    m_ShMod.addConstr(z[i,j,s]>=x[i,s]+x[j,s]-1)                   
    for i in nodes:
        for j in nodes:
            if(i<j):
                for s in range(1,nc+1):
                    m_ShMod.addConstr(z[i,j,s]<=y[i,j])
                m_ShMod.addConstr(y[i,j]<=sum([z[i,j,s] for s in range(1,nc+1)]))
        
    
    # Solve the model
    #m_ShMod.setParam("LogToConsole",0)
    m_ShMod.optimize()
    
    partition=[]
    for s in range(1,nc+1): 
        community=[]
        for i in  nodes:
            if(x[i,s].x>=0.5):
                community.append(i)
        if(len(community)==0):
            break
        else:
            partition.append(community)
    # Return objective value and final partition
    return m_ShMod.objVal,partition

def find_stable_partition(W,nodes,nc,p,banned_partitions):
    
    m_stable = gp.Model("model stable")
    n=len(nodes)
    # Variables
    x={}
    z={}
    y={}
    
    for s in nodes:
        for i in nodes:
            if(i<=s):
                x[i,s]=m_stable.addVar(vtype=GRB.BINARY,lb=0,ub=1,name="x_"+str(i)+"_"+str(s))

        
    
    # Objective function
    m_stable.setObjective(1, GRB.MAXIMIZE)
    
    # Constraints
    m_stable.addConstr(sum([x[s,s] for s in nodes])<=nc)
    for i in nodes:
        m_stable.addConstr(sum([x[i,s] for s in nodes if(i<=s)])>=1)
        m_stable.addConstr(sum([x[i,s] for s in nodes if(i<=s)])<=p)
        for s in nodes:
            if(i<=s):
               m_stable.addConstr(x[i,s]<=x[s,s]) 
               m_stable.addConstr(sum([x[j,s]*W[i,j] for j in nodes if(j!=i and j<=s)])>=x[i,s]*sum([W[i,j] for j in nodes if(j!=i)])/2+(1-x[i,s])*sum([W[i,j] for j in nodes if(W[i,j]<0 and j!=i and j<=s)]))
                
        
    for sol in banned_partitions:
        m_stable.addConstr(sum([sum([x[i,s] for i in nodes if(i<=s and sol[i,s]>=0.5)])+sum([1-x[i,s] for i in nodes if(i<=s and sol[i,s]<0.5)]) for s in nodes])<=(n*n+n)/2-1)
    # Solve the model
    # m_stable.setParam("LogToConsole",0)
    m_stable.optimize()
    
    partition=[]
    for s in nodes: 
        if(x[s,s].x>=0.5):
            community=[s]
            for i in  nodes:
                if(i<s and x[i,s].x>=0.5):
                    community.append(i)
            partition.append(community)
    while(len(partition)<nc):
        partition.append([])
    sol={}
    for i in nodes:
        for s in nodes:
            if(i<=s):
                sol[i,s]=x[i,s].x
        
    # Return objective value and final partition
    return sol,partition

def find_stable_partition2(W,nodes,nc,p):
    
    m_stable = gp.Model("model stable")
    
    # Variables
    x={}
    
    for s in range(1,nc+1):
        for i in nodes:
            x[i,s]=m_stable.addVar(vtype=GRB.BINARY,lb=0,ub=1,name="x_"+str(i)+"_"+str(s))

        
    
    # Objective function
    m_stable.setObjective(1, GRB.MAXIMIZE)
    
    # Constraints
    for i in nodes:
        m_stable.addConstr(sum([x[i,s] for s in range(1,nc+1)])>=1)
        m_stable.addConstr(sum([x[i,s] for s in range(1,nc+1)])<=p)
    for s in range(1,nc+1):
        for i in nodes:
            m_stable.addConstr(sum([x[j,s]*W[i,j] for j in nodes if(j!=i)])>=x[i,s]*sum([W[i,j] for j in nodes if(j!=i)])/2+(1-x[i,s])*sum([W[i,j] for j in nodes if(W[i,j]<0 and j!=i)]))
    
    
    # Solve the model
    m_stable.setParam("LogToConsole",0)
    m_stable.optimize()
    
    partition=[]
    for s in range(1,nc+1): 
        community=[]
        for i in nodes:
            if(x[i,s].x>=0.5):
                community.append(i)
        partition.append(community)
    while(len(partition)<nc):
        partition.append([])
    
        
    # Return objective value and final partition
    return partition

def find_stable_partition3(W,nodes,nc,p,banned_communities):
    
    m_stable = gp.Model("model stable")
    n=len(nodes)
    # Variables
    x={}
    
    for s in nodes:
        for i in nodes:
            if(i<=s):
                x[i,s]=m_stable.addVar(vtype=GRB.BINARY,lb=0,ub=1,name="x_"+str(i)+"_"+str(s))

        
    
    # Objective function
    m_stable.setObjective(1, GRB.MAXIMIZE)
    
    # Constraints
    m_stable.addConstr(sum([x[s,s] for s in nodes])<=nc)
    for i in nodes:
        m_stable.addConstr(sum([x[i,s] for s in nodes if(i<=s)])>=1)
        m_stable.addConstr(sum([x[i,s] for s in nodes if(i<=s)])<=p)
        for s in nodes:
            if(i<s):
               m_stable.addConstr(x[i,s]<=x[s,s]) 
    for s in nodes:
        for i in nodes:
            if(i<=s):
                m_stable.addConstr(sum([x[j,s]*W[i,j] for j in nodes if(j!=i and j<=s)])>=x[i,s]*sum([W[i,j] for j in nodes if(j!=i)])/2+(1-x[i,s])*sum([W[i,j] for j in nodes if(W[i,j]<0 and j!=i and j<=s)]))
        
    for com in banned_communities:
        m_stable.addConstr(sum([x[i,max(com)] for i in com])<=len(com)-1)
    # Solve the model
    m_stable.setParam("LogToConsole",0)
    m_stable.optimize()
    
    partition=[]
    for s in nodes: 
        if(x[s,s].x>=0.5):
            community=[s]
            for i in  nodes:
                if(i<s and x[i,s].x>=0.5):
                    community.append(i)
            partition.append(community)
    while(len(partition)<nc):
        partition.append([])
    sol={}
    for i in nodes:
        for s in nodes:
            if(i<=s):
                sol[i,s]=x[i,s].x
        
    # Return objective value and final partition
    return partition
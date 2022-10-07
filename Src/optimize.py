# Import packages
import gurobipy as gp
from gurobipy import GRB
import time

def cb(model, where):
    if where == GRB.Callback.MIPNODE:
        # Get model objective
        obj = model.cbGet(GRB.Callback.MIPNODE_OBJBST)

        # Has objective changed?
        if abs(obj - model._cur_obj) > 1e-8:
            # If so, update incumbent and time
            model._cur_obj = obj
            model._time = time.time()

    # Terminate if objective has not improved in 1s
    if time.time() - model._time > 1:
        model.terminate()

def model(G_clean, transcripts:list, norm, sparsity_constr, factor:int):

    # Extract edge type and counts from graph file into dictionary
    edges_dict = {}
    for edgeKey, edgeValue in G_clean.edges.items():
        count = edgeValue["counts"]["c"]
        if edgeValue["type"] == "SpliceJunction" or edgeValue["type"] == "Exon":
            edges_dict[edgeKey] = count
    edges = list(edges_dict.keys())

    # Create Adjazenzdictionary (path,edge): 0/1
    adj_dict = {}
    for i in range(0, len(transcripts)):
        for j in range(1, len(transcripts[i]) - 1):
            startnode = transcripts[i][j]
            endnode = transcripts[i][j + 1]
            current_edge = (startnode, endnode)
            if current_edge in edges:
                adj_dict[i, current_edge] = 1
        for edge in edges:
            if (i, edge) not in adj_dict.keys():
                adj_dict[i, edge] = 0

    # Create model for gurobi
    model = gp.Model("Transcript Expression")
    model.Params.LogToConsole = 0

    # Add variables
    no_trans = len(transcripts)
    vars = model.addVars(no_trans, vtype=GRB.INTEGER, name="expression_levels", lb=0) # Expression levels of transcripts
    helper1 = model.addVars(edges, lb=-GRB.INFINITY, vtype=GRB.CONTINUOUS, name="X") # "X" values

    # Optimization
    if sparsity_constr == "L0":
        sparsity_norm0 = model.addVar(name="Sparsity_Constraint_L0") # Sparsity constraint L0 norm
        model.addGenConstrNorm(sparsity_norm0, vars, 0)

    if norm == "L1":
        # Variable
        norm1 = model.addVars(edges, vtype=GRB.CONTINUOUS, name="L1_norm") # L1 norm (Betrag)
        # Objective function
        for j in edges:
            model.addConstr(helper1[j] == (edges_dict[j] - (gp.quicksum(adj_dict[i,j] * vars[i] for i in range(len(transcripts))))))
            model.addConstr(norm1[j] >= helper1[j])
            model.addConstr(norm1[j] >= -helper1[j])
        # Sparsity constraints
        if sparsity_constr == "L0":
            model.setObjective((gp.quicksum(norm1[j] for j in edges) + (factor*sparsity_norm0)), GRB.MINIMIZE)
        elif sparsity_constr == "L1":
            model.setObjective((gp.quicksum(norm1[j] for j in edges) + (factor*vars.sum())), GRB.MINIMIZE)
        else:
            model.setObjective(gp.quicksum(norm1[j] for j in edges), GRB.MINIMIZE)

    elif norm == "L0":
        # Variable
        norm0 = model.addVar(name="L0_norm") # L0 norm
        # Objective function
        for j in edges:
            model.addConstr(helper1[j] == (edges_dict[j] - (gp.quicksum(adj_dict[i,j] * vars[i] for i in range(len(transcripts))))))
        model.addGenConstrNorm(norm0, helper1, 0)
        # Sparsity constraints
        if sparsity_constr == "L0":
            model.setObjective(norm0 + (factor*sparsity_norm0), GRB.MINIMIZE)
        elif sparsity_constr == "L1":
            model.setObjective(norm0 + (factor*vars.sum()), GRB.MINIMIZE)
        else:
            model.setObjective(norm0, GRB.MINIMIZE)

    elif norm == "L2":
        # Variables
        norm2_1 = model.addVars(edges, vtype=GRB.CONTINUOUS, name="L2_norm_1") # Absolute value of helper1
        norm2_2 = model.addVars(edges, vtype=GRB.CONTINUOUS, name="L2_norm_2") # Square of norm2_1
        #norm2_3 = model.addVar(name="L2_norm_3") # Sum of norm2_2
        #norm2_4 = model.addVar(name="L2_norm_4") # Square root of norm2_3 = L2 norm
        # Objective function
        for j in edges:
            model.addConstr(helper1[j] == (edges_dict[j] - (gp.quicksum(adj_dict[i,j] * vars[i] for i in range(len(transcripts))))))
            model.addConstr(norm2_1[j] >= helper1[j])
            model.addConstr(norm2_1[j] >= -helper1[j])
            model.addGenConstrPow(norm2_1[j], norm2_2[j], 2)
        #model.addConstr(norm2_3 == gp.quicksum(norm2_2[j] for j in edges))
        #model.addGenConstrPow(norm2_3, norm2_4, 0.5)
        # Sparsity constraints
        if sparsity_constr == "L0":
            model.setObjective((gp.quicksum(norm2_2[j] for j in edges) + (factor*sparsity_norm0)), GRB.MINIMIZE)
        elif sparsity_constr == "L1":
            model.setObjective((gp.quicksum(norm2_2[j] for j in edges) + (factor*vars.sum())), GRB.MINIMIZE)
        else:
            model.setObjective(gp.quicksum(norm2_2[j] for j in edges), GRB.MINIMIZE)

    model._cur_obj = float('inf')
    model._time = time.time()
    model.optimize(callback=cb)
    
    # Return results
    var_dict = {}
    try:
        for var in model.getVars():
            if "expression_levels" in var.varName:
                var_dict[var.varName[18:len(var.varName)-1]] = var.X
        return var_dict
    except AttributeError:
        return None




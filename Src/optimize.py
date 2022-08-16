# Import packages
import gurobipy as gp
from gurobipy import GRB


def model(G_clean, transcripts, norm, sparsity_constr, factor):
    # Extract edge type and counts from graph file into dictionary
    edges_dict = {}
    for edgeKey, edgeValue in G_clean.edges.items():
        count = edgeValue["counts"]["c"]
        if edgeValue["type"] == "SpliceJunction" or edgeValue["type"] == "Exon":
            edges_dict[edgeKey] = count
    # Print(edges_dict)
    edges = list(edges_dict.keys())

    # Create Adjazenzmatrix (path,edge): 0/1
    adj_matrix = {}
    for i in range(0, len(transcripts)):
        for j in range(1, len(transcripts[i]) - 1):
            startnode = transcripts[i][j]
            endnode = transcripts[i][j + 1]
            current_edge = (startnode, endnode)
            if current_edge in edges:
                adj_matrix[i, current_edge] = 1
        for edge in edges:
            if (i, edge) not in adj_matrix.keys():
                adj_matrix[i, edge] = 0

    # Create model for gurobi
    model = gp.Model("Transcript Expression")
    model.Params.LogToConsole = 0

    # Add variables
    no_trans = len(transcripts)
    vars = model.addVars(no_trans, vtype=GRB.CONTINUOUS, name="expression_levels", lb=0.0) # Expression levels of transcripts
    helper1 = model.addVars(edges, lb=-GRB.INFINITY, vtype=GRB.CONTINUOUS, name="X") # "X" values

    # Optimization
    if sparsity_constr == "L0":
        sparsity_norm0 = model.addVar(name="Sparsity_Constraint_L0") # Sparsity constraint L0 norm
        model.addGenConstrNorm(sparsity_norm0, vars, 0) # Hinzufügen dieser Zeile verändert sofort die Werte der Variablen var1, warum?
        # model.setObjective((gp.quicksum(helper2[j] for j in edges) + (factor*sparsity_norm0)), GRB.MINIMIZE)
        model.addConstr(sparsity_norm0 <= factor)
    elif sparsity_constr == "L1":
        model.addConstr(vars.sum() <= factor) # Sparsity constraint L1 norm

    if norm == "L1":
        norm1 = model.addVars(edges, vtype=GRB.CONTINUOUS, name="L1_norm") # L1 norm (Betrag)
        for j in edges:
            model.addConstr(helper1[j] == (edges_dict[j] - (gp.quicksum(adj_matrix[i,j] * vars[i] for i in range(len(transcripts))))))
            model.addConstr(norm1[j] >= helper1[j])
            model.addConstr(norm1[j] >= -helper1[j])
        model.setObjective(gp.quicksum(norm1[j] for j in edges), GRB.MINIMIZE)
    elif norm == "L0":
        norm0 = model.addVar(name="L0_norm") # L0 norm
        for j in edges:
            model.addConstr(helper1[j] == (edges_dict[j] - (gp.quicksum(adj_matrix[i,j] * vars[i] for i in range(len(transcripts))))))
        model.addGenConstrNorm(norm0, helper1, 0)
        model.setObjective(norm0, GRB.MINIMIZE)
    elif norm == "L2":
        norm2_1 = model.addVars(edges, vtype=GRB.CONTINUOUS, name="L2_norm_1") # Absolute value of helper1
        norm2_2 = model.addVars(edges, vtype=GRB.CONTINUOUS, name="L2_norm_2") # Square of norm2_1
        norm2_3 = model.addVar(name="L2_norm_3") # Sum of norm2_2
        norm2_4 = model.addVar(name="L2_norm_4") # Square root of norm2_3 = L2 norm
        for j in edges:
            model.addConstr(helper1[j] == (edges_dict[j] - (gp.quicksum(adj_matrix[i,j] * vars[i] for i in range(len(transcripts))))))
            model.addConstr(norm2_1[j] >= helper1[j])
            model.addConstr(norm2_1[j] >= -helper1[j])
            model.addGenConstrPow(norm2_1[j], norm2_2[j], 2)
        model.addConstr(norm2_3 == gp.quicksum(norm2_2[j] for j in edges))
        model.addGenConstrPow(norm2_3, norm2_4, 0.5)
        model.setObjective(norm2_4, GRB.MINIMIZE)

    model.optimize()

    # Return results
    var_dict = {}
    for var in model.getVars():
        #print(var.varName)
        #print(var.X)
        if "expression_levels" in var.varName:
            var_dict[var.varName] = var.X
    return var_dict




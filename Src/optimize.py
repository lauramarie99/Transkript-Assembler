# Import packages
import gurobipy as gp
from gurobipy import GRB


def model(G_clean, transcripts):
    #extract edge type and counts from graph file into dictionary
    edges_dict = {}
    for edgeKey, edgeValue in G_clean.edges.items():
        count = edgeValue["counts"]["c"]
        if edgeValue["type"] == "SpliceJunction" or edgeValue["type"] == "Exon":
            edges_dict[edgeKey] = count
    #print(edges_dict)
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

    #create model for gurobi
    model = gp.Model("Transcript Expression")
    model.Params.LogToConsole = 0
    # Add variables
    no_trans = range(len(transcripts))
    vars = model.addVars(no_trans, vtype=GRB.CONTINUOUS, name="expression_levels", lb=0.0)
    #"X" values
    helper1 = model.addVars(edges, lb=-GRB.INFINITY, vtype=GRB.CONTINUOUS, name="X")
    #"Y" values (cannot be negative)
    helper2 = model.addVars(edges, vtype=GRB.CONTINUOUS, name="Y")


    # Define optimization problem
    for j in edges:
        model.addConstr(
            helper1[j] == edges_dict[j] - gp.quicksum(adj_matrix[i, j] * vars[i] for i in range(len(transcripts))))
        # model.addConstr(helper2[j] == gp.abs_(helper1[j]))
        #solution from gatter
        model.addConstr(helper2[j] >= helper1[j])
        model.addConstr(helper2[j] >= -helper1[j])

    model.setObjective(gp.quicksum(helper2[j] for j in edges), GRB.MINIMIZE)

    model.optimize()

    # Print results
    var_dict = {}
    for var in model.getVars():
        #print(var.varName)
        #print(var.X)
        if "expression_levels" in var.varName:
            var_dict[var.varName] = var.X
    return var_dict




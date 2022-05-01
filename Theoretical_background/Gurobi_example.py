from gurobipy import *

# Create empty model
m = Model()

# Add variables
x = m.addVar(vtype=GRB.BINARY, name="x")
y = m.addVar(vtype=GRB.BINARY, name="y")
z = m.addVar(vtype=GRB.BINARY, name="z")

# Set objective function
m.setObjective(x + y + 2*z,GRB.MAXIMIZE)

# Add contraints
c1 = m.addConstr(x + 2*y + 4*z <= 4)
c2 = m.addConstr(x + y >= 1)

# Solve model
m.optimize()

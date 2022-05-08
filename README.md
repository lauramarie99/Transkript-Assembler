# Transcript assembly

The goal of this project is to create a tool for transcript assembly. 
This tool should deliver an expressed set of isoforms and their respective expression levels.

We want to compare different methods and cost functions.

The following tasks have to be completed:

## STEP 1: Preparation

• You are provided with a set of splice graphs and their respective bins and read pairs in a custom format.

• A Python-based parsing script is also provided. You may modify and extend this script for this project.

• The script includes also functionality to write valid GTF-Files for benchmarking.

• Familiarize yourself with the custom graph format and the GTF-Format required as output.

• The project requires the use of Python 3, Gurobi with the Python 3 Interface gurobipy, and NetworkX.
Install all tools and familiarize yourself with them.

• A further script allows to compare the GTF file of an assembly to the also provided truth set. We use this
for benchmarking the assembly quality. In order to execute this script please install the tool cuffcompare.

## STEP 2: Path Enumeration

• Implement the strategies for path enumeration described above, both with and without constraints from
bins and read pairs.

• Include effective "emergency break" mechanisms to avoid excessive runtimes on too complex loci.

• Explore the properties of the computed path-sets. Compute statistics for key properties, including but
not limited to path lengths, set sizes and common edges between path.

• Consider alternative enumeration schemes, such as starting at a central node, which may be more efficient
in combination with bin constraints.

## STEP 3: Optimization with Linear and Quadratic Programming

• Implement the assembly as an application of Linear or Quadratic Programming on the edges of a splice
graph as described above.

• Compare assembly results varying the
  - path enumeration strategies as implemented in WP1,
  - norms,
  - and sparsity constraints.

• Compare results for the complete and pre-filtered splice graphs.

## STEP 4: Flow-based Optimization

• Implement the assembly as an application of the network flow based graph transformation as described above.

• Extend the path enumeration strategies of WP1 for use in flow decomposition.

• Compare assembly results varying the

  - flow decomposition strategy
  - and the used cost function, ideally relating to previously tested norms.

• Compare results for complete and pre-filtered splice graphs.

## STEP 5: Interpretation and Further Research

• Analyze and summarize the results of all workpackages.

• Bonus: Extend any workpackage by your own ideas. You are welcome to change the specs of this project
after consultation with us if you are interested in a "deeper diver" of a particular task.



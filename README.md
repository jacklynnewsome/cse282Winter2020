# cse282Winter2020
CSE assignments



current problem
A cycle that traverses each edge of a graph exactly once is called an Eulerian cycle, and we say that a graph containing such a cycle is Eulerian. The following algorithm constructs an Eulerian cycle in an arbitrary directed graph.

    EULERIANCYCLE(Graph)
        form a cycle Cycle by randomly walking in Graph (don't visit the same edge twice!)
        while there are unexplored edges in Graph
            select a node newStart in Cycle with still unexplored edges
            form Cycle’ by traversing Cycle (starting at newStart) and then randomly walking
            Cycle ← Cycle’
        return Cycle
Eulerian Cycle Problem
Find an Eulerian cycle in a graph.

Given: An Eulerian directed graph, in the form of an adjacency list.

Return: An Eulerian cycle in this graph.

Sample Dataset
0 -> 3

1 -> 0

2 -> 1,6

3 -> 2

4 -> 2

5 -> 4

6 -> 5,8

7 -> 9

8 -> 7

9 -> 6




Sample Output
6->8->7->9->6->5->4->2->1->0->3->2->6

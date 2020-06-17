package uk.ac.ebi.beam;

import java.util.ArrayDeque;
import java.util.Arrays;
import java.util.Deque;

final class KececiogluMatching {

	private final static int EvenLabel = 2;
	private final static int OddLabel = 3;
	private final static int UnreachedLabel = 4;
	private final static  int nil = -1;

	private final Edge[]    Tree;
	private final Edge[]    Bridge;
	private final int[]     Shore;
	private final int[]     Base;
	private final int[]     Label;
	private final int[]     Age;

	private final ArrayDeque<Edge> Path;
	private final ArrayDeque<Integer> alternatingTree;
	
	private final int nMatched;
	private final Matching matching;
	private final IntSet subset;
	private int Time;

	private boolean IsReached(int V) {return Label[V] != UnreachedLabel;}
	private boolean IsEven(int V)    {return Label[V] == EvenLabel;}
	private boolean IsOdd(int V)     {return Label[V] == OddLabel;}
	private int Other(Edge E, int V) {return ((E) == null) ? nil : E.other(V);}

	/**
	 *
	 * Maximum matching in general graphs using Edmond's Blossom Algorithm. This
	 * implementation was adapted from John Kececioglu and Justin Pecqueur.
	 * "Computing maximum-cardinality matchings in sparse general graphs."
	 * Proceedings of the 2nd Workshop on Algorithm Engineering, 121-132, 1998.
	 * and from their C implementation at 
	 * https://www2.cs.arizona.edu/~kece/Research/code/matching1.sh.Z
	 * 
	 * @author Ben Wolfe
	 */
	KececiogluMatching(Graph G, Matching M, int nMatched, IntSet subset)
	{
		this.Time = 1;
		this.Label = new int[G.order()];
		this.Age = new int[G.order()];
		this.Shore = new int[G.order()];
		this.Tree = new Edge[G.order()];
		this.Bridge = new Edge[G.order()];
		this.Base = new int[G.order()];
		this.alternatingTree = new ArrayDeque<Integer>(G.order());
		Path = new ArrayDeque<Edge>(G.order());
		
		
		this.matching = M;
		this.subset =subset;
		
		Arrays.fill(Base, nil);
		Arrays.fill(Label, UnreachedLabel);
		Arrays.fill(Age, 0);

		for (int V = 0; V < G.order(); V++ ) {
			if (subset.contains(V) && matching.unmatched(V) && Search(V, G)) {
				Augment(Path, alternatingTree);
				nMatched += 2;	
			}
		}
		this.nMatched = nMatched;
	}

    /**
     * Find an augmenting path. If an augmenting path
     * was found then the search must be restarted. If a blossom was detected
     * the blossom is contracted and the search continues.
     *
     * @return an augmenting path was found
     */
	boolean Search (int V, Graph G)
	{
		// label current vertex as even and record its age
		Label[V] = EvenLabel;
		Age[V] = Time++;
		boolean Found = false;

		// start alt tree with current vertex
		alternatingTree.addFirst(V);

		// create list of edges connected to current vertex 
		ArrayDeque<Edge> S = new ArrayDeque<Edge>(G.order());
		int W;
		
		// add incident edges to our stack S 
		for (Edge E : G.edges(V)) 
		{
			if (!subset.contains(E.other(V)) || E.bond() == Bond.SINGLE) 
				continue;
			
			S.addLast(E);

			// peek ahead for an augmenting path and bail early if so
			W = E.other(V);
			if (!IsReached(W) && matching.unmatched(W))
				break;
		}
		
		while (!S.isEmpty() && !Found)
		{
			Edge E = S.removeLast();

			int X = Base(E.either());
			int Y = Base(E.other(E.either()));
			if (X == Y)
				continue;
			if (!IsEven(X))
			{
				int Z = X;
				X = Y;
				Y = Z;
			}

			// found an augmenting path
			if (!IsReached(Y) && matching.unmatched(Y))
			{
				Label[Y] = OddLabel;
				Tree[Y] = E;
				Age[Y] = Time++;
				alternatingTree.addFirst(Y);
				Recover(Y);
				Found = true;
				break;

			// found a matched edge, need to add nbrs of mate of edge to DFS	
			} else if (!IsReached(Y) && matching.matched(Y))
			{
				Label[Y] = OddLabel;
				Tree[Y] = E;
				Age[Y] = Time++;
				alternatingTree.addFirst(Y);

				Edge F = G.edge(Y, matching.other(Y));
				int Z = F.other(Y);
				Label[Z] = EvenLabel;
				Age[Z] = Time++;
				alternatingTree.addFirst(Z);

				for (Edge E2: G.edges(Z)) 
				{
					if (E2 != F && subset.contains(E2.other(Z)) && E2.bond() != Bond.SINGLE)
					{
						S.addLast(E2);
						
						// peek ahead for an augmenting path and bail early if so
						W = Other(E2, Z);
						if (!IsReached(W) && matching.unmatched(W))
							break;
					}
				}
			// found a blossom, need to shrink
			} else if (IsEven(Y)) {
				Shrink(E, S, G);
			}
		}

		if (!Found)
		{
			alternatingTree.clear();
		}

		return Found;
	}

	private int Base(int i) {
		return Base[i] == nil ? i : (Base[i] = Base(Base[i]));
	}


	/**
	 *  Recover an augmenting path ending at vertex V by walking
	 *  up the tree back to the root.
	 *
	 * Records a list of the unmatched edges on Path.
	 *
	 */
	private void Recover (int V)
	{
		do
		{
			Path.addFirst(Tree[V]);
			int W = Other(Tree[V], V);
			int B = Base(W);
			Path(W, B, Path);

			V = matching.matched(B) ? matching.other(B) : nil;
		}
		while (V != nil);
	}


	/**
	 * Recursively recover the even-length piece of an alternating path
	 * that begins at vertex V with a matched edge and ends at base B
	 * of its blossom
	 *
	 * The unmatched edges on the path are added to list P, and are in arbitrary
	 * order.
	 *
	 */
	private void Path ( int V, int B, Deque<Edge> P) {

		if (V != B)
			if (IsOdd(V))
			{
				Path(Shore[V], matching.other(V), P);
				P.addFirst(Bridge[V]);
				Path(Other(Bridge[V], Shore[V]), B, P);
			}
			else if (IsEven(V))
			{
				int W = matching.other(V);
				P.addFirst(Tree[W]);
				Path(Other(Tree[W], W), B, P);
			}
	}

	/**
	 * Given an edge E between two even blossoms, shrink the implied
	 * cycle in the alternating tree into a superblossom
	 *
	 * Edges incident to odd vertices on the blossom are added to the stack S
	 * of search edges.
	 *
	 */
	private void Shrink(Edge E, ArrayDeque<Edge> S, Graph G)
	{
		boolean    Found;
		int V = E.either();
		int W = E.other(V);
		int baseV = Base(V);
		int baseW = Base(W);

		if (Age[baseW] > Age[baseV])
		{
			int temp = baseW;
			baseW = baseV;
			baseV = temp;

			temp = V;
			V = W;
			W = temp;
		}

		/*
		 * Walk up the alternating tree from vertex V to vertex A, shrinking
		 * the blossoms into a superblossom.  Edges incident to the odd vertices
		 * on the path from V to A are pushed onto stack S, to later search from.
		 */
		Found = false;
		while (baseV != baseW)
		{
			Base[baseV] = baseW;
			
			// matched edge of B, 1 step back in DFS
			Edge M = G.edge(baseV, matching.other(baseV));
			// node previous to B in DFS, the odd node
			W = Other(M, baseV);

			// bridge of w is the edge after its match in the DFS
			Bridge[W] = E;

			// shore of w is node one step forward in DFS, 
			// not necessarily the same as B which is the base of the blossom of V
			Shore[W] = V;

			// tree of w is the edge leading to w from previous node in DFS
			//Edge T = Tree[W];
			E = Tree[W];

			// look for an unmatched nbr of vertex that is being shrunken into vertex and bail early if we find one
			if (!Found) {
				for (Edge F: G.edges(W)) {
					if (F != M && F != E && subset.contains(F.other(W)) && F.bond() != Bond.SINGLE)
					{
						S.addLast(F);

						int Z = F.other(W);
						// peek ahead and bail early if we know we're going to find an augmenting path
						if (!IsReached(Z) && matching.unmatched(Z))
						{
							Found = true;
							break;
						}
					}
				}
			}

			Base[Base(W)] = baseW;

			// V is now the node before W in the DFS
			V = E.other(W);
			baseV = Base(V);
		}
	}


	/**
	 * Augment the matching along augmenting path P, and expand
	 * into singleton sets all original vertices in T
	 *
	 * This assumes list P contains only the unmatched edges on the path,
	 * and that list T contains all vertices in all blossoms in the alternating
	 * tree containing the augmenting path.
	 *
	 */
	private void Augment(ArrayDeque<Edge> path, ArrayDeque<Integer> alternatingTree) {
		Edge E = path.poll();
		while (E != null)
		{
			matching.match(E.either(), E.other(E.either()));
			E = path.poll();
		}

		Integer V = alternatingTree.poll();
		while (V != null)
		{
		    Label[V] = UnreachedLabel;
			Base[V] = nil;
			V = alternatingTree.poll();
		}
	}
    
	/**
     * Utility to maximise an existing matching of the provided graph.
     *
     * @param g a graph
     * @param m matching on the graph, will me modified
     * @param n current matching cardinality         
     * @param s subset of vertices to match
     * @return the maximal matching on the graph
     */
    static int maximise(Graph g, Matching m, int n, IntSet s) {
    	KececiogluMatching mm = new KececiogluMatching(g, m, n, s);
        return mm.nMatched;
    }
}

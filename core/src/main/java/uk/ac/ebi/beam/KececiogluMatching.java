package uk.ac.ebi.beam;

import java.util.ArrayDeque;
import java.util.Arrays;
import java.util.Deque;

final class KececiogluMatching {

	private final static int EvenLabel = 2;
	private final static int OddLabel = 3;
	private final static int UnreachedLabel = 4;
	private final static  int nil = -1;

	private final Edge[]    tree;
	private final Edge[]    bridge;
	private final int[]     shore;
	private final int[]     base;
	private final int[]     label;
	private final int[]     age;

	private final ArrayDeque<Edge> path;
	private final ArrayDeque<Integer> alternatingTree;
	private final ArrayDeque<Edge> searchStack;
	
	private final int nMatched;
	private final Matching matching;
	private final IntSet subset;
	private int time;

	private boolean isReached(int v) {return label[v] != UnreachedLabel;}
	private boolean isEven(int v)    {return label[v] == EvenLabel;}
	private boolean isOdd(int v)     {return label[v] == OddLabel;}
	private int other(Edge e, int v) {return ((e) == null) ? nil : e.other(v);}

	/**
	 *
	 * Maximum matching in general graphs using Edmond's Blossom Algorithm. This
	 * implementation was adapted from John Kececioglu and Justin Pecqueur.
	 * "Computing maximum-cardinality matchings in sparse general graphs."
	 * Proceedings of the 2nd Workshop on Algorithm Engineering, 121-132, 1998.
	 * and from their C implementation at 
	 * https://www2.cs.arizona.edu/~kece/Research/code/matching1.sh.z
	 * 
	 * @author Ben Wolfe
	 */
	KececiogluMatching(Graph g, Matching M, int nMatched, IntSet subset)
	{
		this.time = 1;
		this.label = new int[g.order()];
		this.age = new int[g.order()];
		this.shore = new int[g.order()];
		this.tree = new Edge[g.order()];
		this.bridge = new Edge[g.order()];
		this.base = new int[g.order()];
		this.alternatingTree = new ArrayDeque<Integer>(g.order());
		this.path = new ArrayDeque<Edge>(g.order());
		this.searchStack = new ArrayDeque<Edge>(g.order());
		
		this.matching = M;
		this.subset =subset;
		
		Arrays.fill(base, nil);
		Arrays.fill(label, UnreachedLabel);
		Arrays.fill(age, 0);

		for (int v = 0; v < g.order(); v++ ) {
			if (subset.contains(v) && matching.unmatched(v) && search(v, g)) {
				augment(path, alternatingTree);
				nMatched += 2;	
			}
		}
		this.nMatched = nMatched;
	}

    /**
     * Find an augmenting path. If an augmenting path
     * was found then the search must be restarted. If a blossom was detected
     * the blossom is contracted and the search continues.
     * @param v a free node in the graph from which to begin depth-first search
     * @param g the graph in which to begin depth-first search
     * @return an augmenting path was found
     */
	private boolean search (int v, Graph g)
	{
		// label current vertex as even and record its age
		label[v] = EvenLabel;
		age[v] = time++;
		boolean found = false;

		// start alt tree with current vertex
		alternatingTree.addFirst(v);

		// create list of edges connected to current vertex 
		searchStack.clear();
		int w;
		
		// add incident edges to our stack S 
		for (Edge e : g.edges(v)) 
		{
			if (!subset.contains(e.other(v)) || e.bond() == Bond.SINGLE) 
				continue;
			
			searchStack.addLast(e);

			// peek ahead for an augmenting path and bail early if so
			w = e.other(v);
			if (!isReached(w) && matching.unmatched(w))
				break;
		}
		
		while (!searchStack.isEmpty() && !found)
		{
			Edge e = searchStack.removeLast();

			int X = find(e.either());
			int y = find(e.other(e.either()));
			if (X == y)
				continue;
			if (!isEven(X))
			{
				int z = X;
				X = y;
				y = z;
			}

			// found an augmenting path
			if (!isReached(y) && matching.unmatched(y))
			{
				label[y] = OddLabel;
				tree[y] = e;
				age[y] = time++;
				alternatingTree.addFirst(y);
				recover(y);
				found = true;
				break;

			// found a matched edge, need to add nbrs of mate of edge to DFS	
			} else if (!isReached(y) && matching.matched(y))
			{
				label[y] = OddLabel;
				tree[y] = e;
				age[y] = time++;
				alternatingTree.addFirst(y);

				Edge f = g.edge(y, matching.other(y));
				int z = f.other(y);
				label[z] = EvenLabel;
				age[z] = time++;
				alternatingTree.addFirst(z);

				for (Edge e2: g.edges(z)) 
				{
					if (e2 != f && subset.contains(e2.other(z)) && e2.bond() != Bond.SINGLE)
					{
						searchStack.addLast(e2);
						
						// peek ahead for an augmenting path and bail early if so
						w = other(e2, z);
						if (!isReached(w) && matching.unmatched(w))
							break;
					}
				}
			// found a blossom, need to shrink
			} else if (isEven(y)) {
				shrink(e, g);
			}
		}

		if (!found)
		{
			alternatingTree.clear();
		}

		return found;
	}

	private int find(int i) {
		return base[i] == nil ? i : (base[i] = find(base[i]));
	}


	/**
	 *  Recover an augmenting path ending at vertex v by walking
	 *  up the tree back to the root.
	 *
	 * Records a list of the unmatched edges on Path.
	 * @param v the free node found indicating the discovery of an augmenting path
	 */
	private void recover (int v)
	{
		do
		{
			path.addFirst(tree[v]);
			int w = other(tree[v], v);
			int b = find(w);
			path(w, b, path);

			v = matching.matched(b) ? matching.other(b) : nil;
		}
		while (v != nil);
	}


	/**
	 * Recursively recover the even-length piece of an alternating path
	 * that begins at vertex v with a matched edge and ends at base b
	 * of its blossom
	 *
	 * The unmatched edges on the path are added to list P, and are in arbitrary
	 * order.
	 *@param v a node in the graph from which to walk to the root
	 *@param b the root to walk back to
	 *@param P a stack on which to place the unmatched edges along the walk
	 */
	private void path ( int v, int b, Deque<Edge> P) {

		if (v != b)
			if (isOdd(v))
			{
				path(shore[v], matching.other(v), P);
				P.addFirst(bridge[v]);
				path(other(bridge[v], shore[v]), b, P);
			}
			else if (isEven(v))
			{
				int w = matching.other(v);
				P.addFirst(tree[w]);
				path(other(tree[w], w), b, P);
			}
	}

	/**
	 * Given an edge e between two even blossoms, shrink the implied
	 * cycle in the alternating tree into a superblossom
	 *
	 * Edges incident to odd vertices on the blossom are added to the stack S
	 * of search edges.
	 *	@param e a blossom-closing edge
	 *  @param g the graph containing the discovered blossom
	 */
	private void shrink(Edge e, Graph g)
	{
		boolean    found;
		int v = e.either();
		int w = e.other(v);
		int baseV = find(v);
		int baseW = find(w);

		if (age[baseW] > age[baseV])
		{
			int temp = baseW;
			baseW = baseV;
			baseV = temp;

			temp = v;
			v = w;
			w = temp;
		}

		/*
		 * Walk up the alternating tree from vertex v to vertex a, shrinking
		 * the blossoms into a superblossom.  Edges incident to the odd vertices
		 * on the path from v to a are pushed onto stack S, to later search from.
		 */
		found = false;
		while (baseV != baseW)
		{
			base[baseV] = baseW;
			
			// matched edge of b, 1 step back in DFS
			Edge m = g.edge(baseV, matching.other(baseV));
			// node previous to b in DFS, the odd node
			w = other(m, baseV);

			// bridge of w is the edge after its match in the DFS
			bridge[w] = e;

			// shore of w is node one step forward in DFS, 
			// not necessarily the same as b which is the base of the blossom of v
			shore[w] = v;

			// tree of w is the edge leading to w from previous node in DFS
			//Edge T = Tree[w];
			e = tree[w];

			// look for an unmatched nbr of vertex that is being shrunken into vertex and bail early if we find one
			if (!found) {
				for (Edge f: g.edges(w)) {
					if (f != m && f != e && subset.contains(f.other(w)) && f.bond() != Bond.SINGLE)
					{
						searchStack.addLast(f);

						int z = f.other(w);
						// peek ahead and bail early if we know we're going to find an augmenting path
						if (!isReached(z) && matching.unmatched(z))
						{
							found = true;
							break;
						}
					}
				}
			}

			base[find(w)] = baseW;

			// v is now the node before w in the DFS
			v = e.other(w);
			baseV = find(v);
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
	 * @param path a stack of unmatched edges along the augmenting path
	 * @param alternatingTree a stack of nodes in all the blossoms in the alternating tree
	 * containing the augmenting path
	 */
	private void augment(ArrayDeque<Edge> path, ArrayDeque<Integer> alternatingTree) {
		Edge e = path.poll();
		while (e != null)
		{
			matching.match(e.either(), e.other(e.either()));
			e = path.poll();
		}

		Integer v = alternatingTree.poll();
		while (v != null)
		{
		    label[v] = UnreachedLabel;
			base[v] = nil;
			v = alternatingTree.poll();
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
    
    /**
     * Utility to maximise an existing matching of the provided graph.
     *
     * @param g a graph
     * @param m matching on the graph
     * @return the maximal matching on the graph
     */
    static int maximise(Graph g, Matching m, int n) {
        return maximise(g, m, n, IntSet.universe());
    }

    /**
     * Utility to get the maximal matching of the specified graph.
     *
     * @param g a graph
     * @return the maximal matching on the graph
     */
    static Matching maximal(Graph g) {
        Matching m = Matching.empty(g);
        maximise(g, m, 0);
        return m;
    }
}

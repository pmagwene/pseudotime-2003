"""Algorithms for PQ-tree construction

Requires Python 2.3+

-----------------------------------------------------------------------------
(c) Copyright by Paul M. Magwene, 2002-2004  (mailto:pmagwene@sas.upenn.edu)

    Permission to use, copy, modify, and distribute this software and its
    documentation for any purpose and without fee or royalty is hereby granted,
    provided that the above copyright notice appear in all copies and that
    both that copyright notice and this permission notice appear in
    supporting documentation or portions thereof, including modifications,
    that you make.

    THE AUTHOR PAUL M. MAGWENE DISCLAIMS ALL WARRANTIES WITH REGARD TO
    THIS SOFTWARE, INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND
    FITNESS, IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY SPECIAL,
    INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING
    FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT,
    NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION
    WITH THE USE OR PERFORMANCE OF THIS SOFTWARE !
-----------------------------------------------------------------------------
"""


# standard python modules
import copy, random

# third party modules
import Numeric as Num

# local modules
import graph



#------------------------------------------------------------------------------
#   ____                               _     ____             _   _
#  / ___| _   _ _ __  _ __   ___  _ __| |_  |  _ \ ___  _   _| |_(_)_ __   ___  ___
#  \___ \| | | | '_ \| '_ \ / _ \| '__| __| | |_) / _ \| | | | __| | '_ \ / _ \/ __|
#   ___) | |_| | |_) | |_) | (_) | |  | |_  |  _ < (_) | |_| | |_| | | | |  __/\__ \
#  |____/ \__,_| .__/| .__/ \___/|_|   \__| |_| \_\___/ \__,_|\__|_|_| |_|\___||___/
#              |_|   |_|
#------------------------------------------------------------------------------

def distance_matrix(data):
    """data points -> distance matrix
    
    Given an n x p data matrix (variables in columns, objects in rows),
    calculcate the euclidean distance between each object. Returns
    a n x n matrix.
    """
    Q = Num.innerproduct(data,data)
    diag = Num.diagonal(Q)
    R = Num.ones(Q.shape)*diag
    return Num.sqrt(R+Num.transpose(R) - 2*Q)
    

# Priority dictionary using binary heaps
# David Eppstein, UC Irvine, 8 Mar 2002

class PriorityDict(dict):
    def __init__(self):
        """Initialize PriorityDict by creating binary heap of pairs (value,key). 
        Note that changing or removing a dict entry will not remove the old pair 
        from the heap until it is found by smallest() or until the heap is 
        rebuilt.
        """        
        self.__heap = []
        dict.__init__(self)

    def smallest(self):
        '''Find smallest item after removing deleted items from heap.'''
        if len(self) == 0:
            raise IndexError, "smallest of empty PriorityDictionary"
        heap = self.__heap
        while heap[0][1] not in self or self[heap[0][1]] != heap[0][0]:
            lastItem = heap.pop()
            insertionPoint = 0
            while 1:
                smallChild = 2*insertionPoint+1
                if smallChild+1 < len(heap) and \
                        heap[smallChild] > heap[smallChild+1]:
                    smallChild += 1
                if smallChild >= len(heap) or lastItem <= heap[smallChild]:
                    heap[insertionPoint] = lastItem
                    break
                heap[insertionPoint] = heap[smallChild]
                insertionPoint = smallChild
        return heap[0][1]
    
    def __iter__(self):
        '''Create destructive sorted iterator of PriorityDictionary.'''
        def iterfn():
            while len(self) > 0:
                x = self.smallest()
                yield x
                del self[x]
        return iterfn()
    
    def __setitem__(self,key,val):        
        """Change value stored in dictionary and add corresponding pair to heap. 
        Rebuilds the heap if the number of deleted items grows too large, to 
        avoid memory leakage.
        """
        dict.__setitem__(self,key,val)
        heap = self.__heap
        if len(heap) > 2 * len(self):
            self.__heap = [(v,k) for k,v in self.iteritems()]
            self.__heap.sort()  # builtin sort likely faster than O(n) heapify
        else:
            newPair = (val,key)
            insertionPoint = len(heap)
            heap.append(None)
            while insertionPoint > 0 and \
                    newPair < heap[(insertionPoint-1)//2]:
                heap[insertionPoint] = heap[(insertionPoint-1)//2]
                insertionPoint = (insertionPoint-1)//2
            heap[insertionPoint] = newPair
    
    def setdefault(self,key,val):
        '''Reimplement setdefault to call our customized __setitem__.'''
        if key not in self:
            self[key] = val
        return self[key]

    def update(self, other):
        for key in other.keys():
            self[key] = other[key]




#----------------------------------------------------------------------------
#      __  __ ____ _____
#     |  \/  / ___|_   _|
#     | |\/| \___ \ | |
#     | |  | |___) || |
#     |_|  |_|____/ |_|
#
# Prim's algorith for finding MST
#       -- implemented in pure python
#----------------------------------------------------------------------------


def adjacent(W, u):
    return list(Num.nonzero(Num.greater(W[u],0)))

def mst_prim(V, W):
    """Prim's algorithm for finding the minimum spanning tree (MST) of a graph.

    Arguments:
        * V is an integer list of vertices of the graph
        * W is a weight matrix with non-zero edge weights indicating an edge.

    Returns:
        * a list of edges expressed as vertex pairs

    - Based on an original pure python implementation by Tim Wilson.
    - Sped up considerably by using PriorityDict class of D. Eppstein
    - For full discussion of algorithm see Cormen et al. 2001.
        Introduction to algorithms. MIT Press.
    """
    W = Num.array(W, Num.Float)
    INF = 10000 * max(max(W))

    # initialize and set each value of the array P (pi) to none
    # pi holds the parent of u, so P(v)=u means u is the parent of v
    P = {}

    Q = PriorityDict()
    Q[V[0]] =  0       # make the first vertex the root

    i = 1 
    for v in V[1:]:
        Q[i] = INF
        i += 1

    # loop while the min queue is not empty
    for u in Q:
        Adj = adjacent(W, u)
        for v in Adj:
            w = W[u,v]      # get the weight of the edge uv
            # proceed if v is in Q and the weight of uv is less than v's key
            try:
                if w < Q[v]:
                    # set v's parent to u
                    P[v] = u
                    # v's key to the weight of uv
                    Q[v] = w
            except KeyError:
                continue
    
    g = {}
    for key in P:
        i, j = key, P[key]
        graph.add(g, i, j, W[i][j])
    return g


#------------------------------------------------------------------------------    
#   ____  _                _            _     ____       _   _
#  / ___|| |__   ___  _ __| |_ ___  ___| |_  |  _ \ __ _| |_| |__
#  \___ \| '_ \ / _ \| '__| __/ _ \/ __| __| | |_) / _` | __| '_ \
#   ___) | | | | (_) | |  | ||  __/\__ \ |_  |  __/ (_| | |_| | | |
#  |____/|_| |_|\___/|_|   \__\___||___/\__| |_|   \__,_|\__|_| |_|
#
#------------------------------------------------------------------------------



def Dijkstra(graph, start, end=None):
    """Find shortest paths from the start vertex to all vertices nearer than or 
    equal to the end.
    
    Dijkstra's algorithm for shortest paths
    David Eppstein, UC Irvine, 4 April 2002
    
    http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/117228    

    The input graph G is assumed to have the following
    representation: A vertex can be any object that can
    be used as an index into a dictionary.  G is a
    dictionary, indexed by vertices.  For any vertex v,
    G[v] is itself a dictionary, indexed by the neighbors
    of v.  For any edge v->w, G[v][w] is the length of
    the edge.  This is related to the representation in
    <http://www.python.org/doc/essays/graphs.html>
    where Guido van Rossum suggests representing graphs
    as dictionaries mapping vertices to lists of neighbors,
    however dictionaries of edges have many advantages
    over lists: they can store extra information (here,
    the lengths), they support fast existence tests,
    and they allow easy modification of the graph by edge
    insertion and removal.  Such modifications are not
    needed here but are important in other graph algorithms.
    Since dictionaries obey iterator protocol, a graph
    represented as described here could be handed without
    modification to an algorithm using Guido's representation.

    Of course, G and G[v] need not be Python dict objects;
    they can be any other object that obeys dict protocol,
    for instance a wrapper in which vertices are URLs
    and a call to G[v] loads the web page and finds its links.
    
    The output is a pair (D,P) where D[v] is the distance
    from start to v and P[v] is the predecessor of v along
    the shortest path from s to v.
    
    Dijkstra's algorithm is only guaranteed to work correctly
    when all edge lengths are positive. This code does not
    verify this property for all edges (only the edges seen
     before the end vertex is reached), but will correctly
    compute shortest paths even for some graphs with negative
    edges, and will raise an exception if it discovers that
    a negative edge has caused it to make a mistake.
    """

    D = {}    # dictionary of final distances
    P = {}    # dictionary of predecessors
    Q = PriorityDict()   # est.dist. of non-final vert.
    Q[start] = 0
    
    for v in Q:
        D[v] = Q[v]
        if v == end: break
        
        for w in graph[v]:
            vwLength = D[v] + graph[v][w]
            if w in D:
                if vwLength < D[w]:
                    raise ValueError, \
  "Dijkstra: found better path to already-final vertex"
            elif w not in Q or vwLength < Q[w]:
                Q[w] = vwLength
                P[w] = v
    
    return (D,P)
            

def shortest_path_dijkstra(graph,start,end):
    """Find a single shortest path from the given start vertex
    to the given end vertex.

    The input has the same conventions as Dijkstra().
    The output is a list of the vertices in order along
    the shortest path.
    
    Returns None on failure.
    """

    D,P = Dijkstra(graph,start,end)
    Path = []
    try:
        while 1:
            Path.append(end)
            if end == start: break
            end = P[end]
    except KeyError:
        return None
    Path.reverse()
    return Path


def allpairs_dijkstra(graph, withprogress=False):
    spdist = {}
    keys = graph.keys()
    
    for key in keys:
        spdist[key] = {}

    for i,key in enumerate(keys):
        D,P = Dijkstra(graph,key)
        for other in keys:
            try:
                spdist[key][other] = D[other]
            except KeyError:
                spdist[key][other] = None 

        if withprogress:
            if not (i % 10):
                print "%d..." % (i),
                    
    return spdist                    



#------------------------------------------------------------------------------
#  ____  _                      _              ____       _   _
# |  _ \(_) __ _ _ __ ___   ___| |_ ___ _ __  |  _ \ __ _| |_| |__
# | | | | |/ _` | '_ ` _ \ / _ \ __/ _ \ '__| | |_) / _` | __| '_ \
# | |_| | | (_| | | | | | |  __/ ||  __/ |    |  __/ (_| | |_| | | |
# |____/|_|\__,_|_| |_| |_|\___|\__\___|_|    |_|   \__,_|\__|_| |_|



def tree_diameter(G, first = None):
    """Given a weighted tree calculates the diameter path.  Returns: (diampath, length).
    
    The diameter path is the longest of all shortest-paths in the tree.
    
    Calculation is via a double hanging path on the tree.

    The tree should be represented in the dictionary of dictionaries format.
    """
    if first is None:
        first = random.choice(G.keys())
    d, p  = Dijkstra(G, first)
    td = zip(d.values(), d.keys())
    dmax = max(td)
    second = dmax[1]

    d,p = Dijkstra(G, second)
    td = zip(d.values(), d.keys())
    dmax = max(td)
    furthest = dmax[1]

    return shortest_path_dijkstra(G, second, furthest), dmax[0]            




#------------------------------------------------------------------------------
#  ____   ___    _____
# |  _ \ / _ \  |_   _| __ ___  ___
# | |_) | | | |   | || '__/ _ \/ _ \
# |  __/| |_| |   | || | |  __/  __/
# |_|    \__\_\   |_||_|  \___|\___|


def indecisive_backbone(G, diameterpath):
    """Given a tree and the diameter path of that tree,
    returns the "indecisive backbone" of the diameter path.

    The indecisive backbone is defined as the largest continuous set of vertices on 
    the diameter path for which the first and last vertices are either both of
    degree >= 3, both of degree == 1, or if there is only one vertex of degree >= 3, then
    that vertex alone.

    Returns None if none of the diameter path nodes are indecisive.

    """

    #find first indecisive vertex
    firstidx = None
    ct = 0
    for vertex in diameterpath:
        if graph.degree(G, vertex) >= 3:
            firstidx = ct
            break
        ct += 1

    if firstidx is None:
        return None

    #find last indecisive vertex
    opp = diameterpath[:]
    opp.reverse()
    lastidx = None

    ct = 0
    for vertex in opp:
        if graph.degree(G, vertex) >= 3:
            lastidx = ct
            break
        ct += 1

    lastidx = diameterpath.index(opp[ct])
    
    return diameterpath[firstidx:lastidx+1]



def diameterpath_branches(tree, diampath):
    """Given a tree and a diameter path, find the branches of the tree that lie
    off of the diameter path.
    """
    tgraph = copy.deepcopy(tree)
    branches = {}
    for node in diampath:
        branches[node] = []
        for n in graph.neighbors(tree, node):
            if n not in diampath:
                graph.deledge(tgraph, n, node)
                reachable, pred = Dijkstra(tgraph, n)
                branches[node].append(reachable.keys())                
    
    subgraphs = {}
    edges = graph.edges(tree)
    weights = [graph.weight(tree, e[0], e[1]) for e in graph.edges(tree)]
    for root in branches.keys():
        subgraphs[root] = []
        for branch in branches[root]:
            subg = graph.subgraph(tree, branch)
            subgraphs[root].append(subg)
    
    class Results:
        pass

    results = Results()
    results.branches = branches
    results.subgraphs = subgraphs
    return results    


class Qnode(list):
    def __repr__(self):
        return str(tuple(self))

class Pnode(list):
    pass


def make_pqtree(mstree):
    """Builds a PQ-tree by walking along the diameter path indecisive backbone
    of a MST.
    """
    # find diameter path
    diampath, dpl = tree_diameter(mstree)

    # find indecisive backbone
    backbone = indecisive_backbone(mstree, diampath)

    if backbone is None:
        return Qnode(diampath)

    # find diampath branches
    dpbranches = diameterpath_branches(mstree, backbone)

    pqtree = []
    #for root in dpbranches.subgraphs.keys():
    for root in backbone:
        pqtree.append(Pnode( [Qnode((root,))] + \
                 Pnode([make_pqtree(s) for s in dpbranches.subgraphs[root]]) ))

    return Qnode(pqtree)
        

def find_pqsets(pqtree):
    """Find independent sets of PQ tree.
    
    Sets of nodes separated by decisive nodes.
    """
    indsets = []
    previous = 0
    for i in range(len(pqtree)):
        if len(pqtree[i]) == 1:
            s = pqtree[previous:i+1]
            if s:
                indsets.append(s)
            indsets.append(pqtree[i])
            previous = i
    # include the last set
    indsets.append(pqtree[previous:])
    return indsets

    
    
def snip_indecisive(mstree, pqtree):
    """Snip's all the indecisive diameter path nodes out of a pqtree"""
    snippedtree = copy.deepcopy(pqtree)
    for i in range(len(pqtree)):
        if graph.degree(mstree, pqtree[i][0][0]) > 2:
            del snippedtree[i][0]
    return snippedtree


#----------------------------------------------------------------------------
# For finding PQ-tree permutations


def permute_iter(seq):
    """Constructs all permutations of the given sequence.
    """
    if len(seq) == 1:
        yield (seq)
        raise StopIteration

    for i in range(len(seq)):
        eslice = seq[i:i+1]
        rest_iter = permute_iter(seq[:i] + seq[i+1:])
        for rest in rest_iter:
            yield (eslice + rest)

    raise StopIteration

    
def qflip_shallow(qnode):
    rev = list(qnode)
    rev.reverse()
    return rev

def pperm_shallow(pnode):
    return list(permute_iter(pnode))

def pqtree_perms(seq):
    stack = [seq[:]]
    ischanging = True
    while ischanging:
        for j in range(len(stack)):
            s = stack[j]
            for i in range(len(s)):
                if type(s[i]) == Pnode:
                    ischanging = True
                    perms = pperm_shallow(s[i])
                    if len(perms) > 1:
                        for perm in perms:
                            if perm != list(s[i]):
                                stack.append(s[:i] + perm + s[i+1:])
                                stack[j] = s[:i] + s[i] + s[i+1:] 
                    else:
                        stack[j] = s[:i] + s[i] + s[i+1:]
                    break
                elif type(s[i]) == Qnode:
                    ischanging = True
                    perm = qflip_shallow(s[i])
                    if perm != list(s[i]):
                        stack.append(s[:i] + perm + s[i+1:])
                        stack[j] = s[:i] + s[i] + s[i+1:]
                    else:
                        stack[j] = s[:i] + s[i] + s[i+1:]
                    break
                else:
                    ischanging = False
                    continue
    return stack




#------------------------------------------------------------------------------
#  ____  _                 _        ___       _             __
# / ___|(_)_ __ ___  _ __ | | ___  |_ _|_ __ | |_ ___ _ __ / _| __ _  ___ ___  ___
# \___ \| | '_ ` _ \| '_ \| |/ _ \  | || '_ \| __/ _ \ '__| |_ / _` |/ __/ _ \/ __|
#  ___) | | | | | | | |_) | |  __/  | || | | | ||  __/ |  |  _| (_| | (_|  __/\__ \
# |____/|_|_| |_| |_| .__/|_|\___| |___|_| |_|\__\___|_|  |_|  \__,_|\___\___||___/
#                    |_|

def mst_and_diameterpath(data):    
    """ data points -> Euclidean MST (as graph), diameter path (list of 
    vertices), and diameter path length.    
    """
    dist = distance_matrix(data)
    mstgraph = mst_prim(range(len(data)), dist)
    diampath = tree_diameter(mstgraph)
    return mstgraph, diampath[0], diampath[1]
    



#----------------------------------------------------------------------------
# path length finding and sorting

def path_length(path, distmtx, seglengths=None):
    """path, distance matrix -> length of path (optional: segment lengths)
    
    Given distance matrix and path (sequence of integers where values correspond
    to rows of matrix; e.g. [1,2,4,5,3,9,8,6]) calculates total path length.
    """ 
    pairs = zip(path,path[1:])
    segments = []
    tot = 0.0
    for each in pairs:
        segments.append(distmtx[each])
        tot += distmtx[each]
    
    if seglengths is None:
        return tot
    else:
        return tot, segments
        

def rank_paths(paths, dist):
    """Given a set of paths, and a distance matrix, rank those paths from shortest
    to longest.
    """
    # rank and order paths
    lengths = [path_length(i,dist) for i in paths]
    spaths = Num.take(paths, Num.argsort(lengths))
    sortedpaths = []
    for path in spaths:
        sortedpaths.append(path.tolist())
    sortedlengths = Num.sort(lengths)
    
    class Results:
        pass
    r = Results()        
    r.paths = sortedpaths
    r.lengths = sortedlengths
    return r    
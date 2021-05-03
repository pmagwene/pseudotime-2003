"""Graph algorithms

Requires Python 2.3+

-----------------------------------------------------------------------------
(c) Copyright by Paul M. Magwene, 2003-2004  (mailto:pmagwene@sas.upenn.edu)

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

Let G be a graph, with vertex set V and edge set E.  These algorithms use a 
dictionary of dictionaries representation of G. This means that G[i][j] 
exists if (i,j) is in the edge set of G.

Some common, non-modifying operations on graphs:

Operation                   Method              Function
---------                   ------              ---------

vertices of graph           graph.keys()        vertices(graph)
(out)neighbors of vertex i  graph[i].keys()     neighbors(graph, i)
in-neighbors of vertex i     none               inneighbors(graph, i)
degree of vertex i          len(graph[i])       degree(graph, i)
edges of graph               none               edges(graph)*
weight of edge (i,j)        graph[i][j]         weight(graph, i, j)
Number of edges in graph     none               numedges(graph)*
Subgraph                     none               subgraph(graph, subverts)
Subgraph on edges            none               edgesubgraph(graph, subedges)     
Connected components         none               connected_components(graph)
Connected component graphs   none               connected_subgraphs(graph)

adj matrix representation    none               asadj(graph)
as unweighted graph          none               asunweighted(graph)
graph from adjmatrix         none               fromadj(adjmtx)
graph from edgelist          none               fromedges(adjmtx)*
adjmatrix from graph         none               asadj(graph)
edgelist from graph          none               asedges(graph)

join:G1,G2 -> G{V1|V2,E1|E2} none               join(graph1,graph2)
line graph                   none               linegraph(graph)


Modifying operations on graphs:

Operation                   Method              Function
---------                   ------              ---------

add vertex i                graph[i]={}         add(graph, i)
add edge (i,j)              graph[i][j]=weight  add(graph, i, j, weight)*
delete edge (i,j)            none               deledge(graph, i, j)*
delete vertex i              none               delvertex(graph, i)





I/O for graphs:

Operation                   Method              Function
---------                   ------              ---------

save in graphviz format      none               tographviz(graph, filename)*
save in pajek format         none               topajek(graph, filename)
save partition in pajek fmt  none               topajekpart(graph, partdict, filename)
save as a text file          none               writegraph(graph, filename)
read from text file          none               readgraph(filename)


Misc support routines:

invertpartition



* has keyword argument for directed/undirected edges, defaults to undirected

Explanations of data structures/concepts

partdict: "Partition Dictionary" - a dictionary representing a partition
of a graph.  The keys of this dictionary are arbitrary partition names
(e.g. 1, "A", "protein", etc).  The values of this dictionary are Sets of 
vertices.

Line graph: Given G, the line graph G' where the vertices of G' are the edges
of G (i.e. V1'=E1}, and (V1',V2') is an edge if E1 and E2 are incident to the 
the same vertex in G.


""" 
__version__ = "0.1"


#------------------------------------------------------------------------------#
import sys, copy, string, sets, cPickle, os
from sets import Set, ImmutableSet


class GraphError(Exception):
    pass

#------------------------------------------------------------------------------#

def add (graph, i, j=None, weight=1, undirected=True, noselfloops=True):
    """Add a vertex or edge to a graph."""
    if j is None:
        if i not in graph.keys():
            graph[i] = {}
    else:
        if (i==j) and noselfloops:
            return
        if i not in graph:
            graph[i] = {}
        if j not in graph:
            graph[j] = {}
        graph[i][j] = weight
        if undirected:
            graph[j][i] = weight

     
def deledge (graph, i, j, undirected=True):
    """Delete edge (i,j) from the graph."""
    del graph[i][j]
    if undirected: del graph[j][i]
    

def delvertex (graph, i):
    """Delete the vertex i from the graph."""
    del graph[i]
    for key in graph.keys():
        try:
            del graph[key][i]
        except KeyError:
            pass                

def vertices (graph):
    """Returns list of vertices of the graph."""
    return graph.keys()
    

def neighbors (graph, i, k=1):
    """Returns (k)neighbors of vertex i.
    
    The k-neighbors of i are all vertices no more than k steps away from i.
    """
    if k == 1:
        return Set(graph[i].keys())
    else:
        neigh = Set(neighbors(graph, i))
        for other in list(neigh):
            neigh |= neighbors(graph, other, k-1)
    try:
        neigh.remove(i)
    except KeyError:
        pass
    return neigh


def neighborset(graph, queryset, k=1):
    """Finds all neighbors of vertices in the query set. 
    
    The neighbor set does not include any elements of the query set.
    """
    neighset = Set()
    for v in queryset:
        neighset |= Set(neighbors(graph,v,k))
    return neighset - Set(queryset)
      
    
def inneighbors (graph, i):
    """Gives the in-neighbors of v in G.
    """
    nin = Set()
    for vert in graph.keys():
        if i in graph[vert]:
            nin.add(vert)
    return nin    

    

def degree (graph, i):
    """Degree of vertex i.
    """
    return len(graph[i])    


def edges (graph, undirected=True):
    """Returns Set of edges."""
    edges = sets.Set()
    seen = sets.Set()
    for key in graph.keys():
        for other in graph[key]:
            if undirected:
                keyother = Set([key,other])
                if keyother not in seen:
                    edges.add((key,other))
                    seen.add(keyother)                
            else:
                edges.add((key,other))
    return edges
    

def sortededges (graph):
    """Returns edges of G, sorted from from strongest to weakest.
    """
    gedges = list(edges(graph))
    weights = []
    idx = 0
    for e in gedges:
        weights.append((graph[e[0]][e[1]], idx))
        idx += 1
    weights.sort()
    sedges = [gedges[w[1]] for w in weights]
    sedges.reverse()
    return sedges



def weight (graph, i, j):
    """Return weight of edge (i,j)."""
    return graph[i][j]


def numedges (graph, undirected=True):
    """Returns number of edges in the graph."""
    nonzero = sum( [len(graph[key]) for key in graph.keys()])
    if undirected: return nonzero/2
    else: nonzero
    

def subgraph (graph, vertices):
    """Get the subgraph determined by the given vertices.
    """
    sub = {}
    for v in vertices:
        sub[v] = {}
        for other in Set(graph[v]) & Set(vertices):
            sub[v][other] = graph[v][other]
    return sub
    
def edgesubgraph (graph, edges, undirected=True):
    """Get the subgraph determined by the given edges.
    """
    sub = {}
    for edge in edges:
        x,y = edge
        add(sub, x, y, graph[x][y], undirected=undirected)
    return sub
        

def join (g1, g2):
    """g1,g2 -> G = {V1 | V2, E1 | E2}.
    """
    g3 = copy.deepcopy(g1)
    for k in g2:
        g3[k] = copy.deepcopy(g2[k])
    return g3


def linegraph(graph):
    """Invert the vertices/edges of the graph.
    
    Edges become vertices; vertices become edges.
    """
    G = {}
    nbrdict = {}
    gedges = edges(graph)
    for edge in gedges:
        x,y = edge
        ekey = (edge, graph[x][y])
        G[ekey] = {}
        nbrdict[ekey] = Set()
        # neighbors of x
        for z in graph[x]:
            if (x,z) in gedges:
                nbrdict[ekey].add(((x,z),graph[x][z]))
            elif (z,x) in gedges:
                nbrdict[ekey].add(((z,x),graph[z][x]))
            else:
                raise GraphError("Problem constructing line graph.")
        # neighbors of y
        for z in graph[y]:
            if (y,z) in gedges:
                nbrdict[ekey].add(((y,z),graph[y][z]))
            elif (z,y) in gedges:
                nbrdict[ekey].add(((z,y),graph[z][y])) 
            else:
                raise GraphError("Problem constructing line graph.")
    
    for ekey in G.keys():
        for nbr in nbrdict[ekey]:
            G[ekey][nbr] = 1

    return G
    
def reversegraph(graph):
    """ If (i,j) is an edge in G, (j,i) is an edge in G'.
    """
    rG = {}
    for i in graph:
        if i not in rG:
            rG[i] = {}
        for j in graph[i]:
            add(rG, j, i, graph[i][j],undirected=False)
    return rG                



def asunweighted (graph):
    """ graph with weighted edges -> graph with all edge weights = 1
    """
    ugraph = copy.deepcopy(graph)
    for v1 in ugraph:
        for v2 in ugraph[v1]:
            ugraph[v1][v2] = 1
    return ugraph



def applyedgefunc (graph, func):
    """ G -> G' where G'[i][j] = func(G[i][j]).
    
    Returns new graph.
    """
    gcopy = copy.deepcopy(graph)
    for v1 in gcopy:
        for v2 in gcopy[v1]:
            gcopy[v1][v2] = func(gcopy[v1][v2])
    return gcopy                


def swapvertices (graph, v1, v2, undirected=True, noselfloops=True):
    """graph, v1, v2 -> edges(v1) = edges(v2), edges(v2) = edges(v1).
    
    Modifies the graph in place.
    """
    v1neigh = copy.deepcopy(graph[v1])
    v2neigh = copy.deepcopy(graph[v2])
    delvertex(graph, v1)
    delvertex(graph, v2)
    add(graph, v1)
    add(graph, v2)
    for other in v1neigh:
        add(graph, v2, other, v1neigh[other], undirected=undirected, noselfloops=True)
    for other in v2neigh:
        add(graph, v1, other, v2neigh[other], undirected=undirected, noselfloops=True)   


def len_cmp (x,y):
    if len(x) > len(y):
        return 1
    elif len(x) < len(y):
        return -1
    else:
        return 0


def connected_components (graph):
    """Find the connected components of a graph.
    
    Each connected component is represented as a list of vertices.
    """
    djset = DisjointSets(graph.keys())
    c = djset.connected_components(edges(graph)).values()
    
    # sort list so biggest connected component is at front
    c.sort(len_cmp)
    c.reverse()
    return c  

def connected_subgraphs (graph):
    """Finds the connected components of a graph.
    
    Return a list of connected components represented as subgraphs.    
    """
    cc = connected_components(graph)
    subs = []
    for sub in cc:    
        subs.append(subgraph(graph, sub))
    return subs
    

class DisjointSets (dict):
    def __init__(self, items):
        self.parent = {}
        self.rank = {}
        
        # make_set
        for item in items:
            self.parent[item] = item
            self.rank[item] = 0

    def link(self, x, y):
        if self.rank[x] > self.rank[y]:
            self.parent[y] = x
        else:
            self.parent[x] = y
            if self.rank[x] == self.rank[y]:
                self.rank[y] += 1

    def union(self, x, y):
        self.link(self.find_set(x), self.find_set(y))

    def find_set(self, x):
        if x is not self.parent[x]:
            self.parent[x] = self.find_set(self.parent[x])
        return self.parent[x]

    def connected_components(self, connections):
        for c in connections:
            if self.find_set(c[0]) != self.find_set(c[1]):
                self.union(c[0], c[1])
        [self.find_set(x) for x in self.parent.keys()]
        by_parent = {}
        for x,y in zip(self.parent.values(), self.parent.keys()):
            if by_parent.has_key(x):
                by_parent[x].append(y)
            else:
                by_parent[x] = [y]
        return by_parent                

    def same_component(self, u, v):
        return self.find_set(u) == self.find_set(v)    
        



#------------------------------------------------------------------------------#
# Matrix based graph representation functions

def zeroarray (dim1, dim2, atype=type(0)):
    """Create array of zeros
    """
    zero = atype(0)
    row = list([zero]*dim2)
    return [row[:] for i in range(dim1)]
    

def asadj (graph):            
    """Return adjacencency matrix representation of graph.
        
    Returns a adjacency mtx representation of the graph plus a list which
    gives the graph labels in the order in which they appear in the adj matrix.
    """
    adj = []   
    a2g = []
    verts = graph.keys()
    verts.sort()
    
    adj = zeroarray(len(verts),len(verts))    
    for i,v in enumerate(verts):        
        a2g.append(v)
        for other in graph[v].keys(): 
            j = verts.index(other)           
            adj[i][j] = graph[v][other]
        
    return adj, a2g        


def fromadj (adjmtx, labels=None, undirected=True, noselfloops=True):
    """Create dict-of-dict graph representation from adjacency matrix.
    
    Adjacency matrix can be weighted.
    """
    nverts = len(adjmtx)
    graph = {}
    for i, row in enumerate(adjmtx):
        ilabel = i
        if labels is not None:
            ilabel = labels[i]
        if ilabel not in graph:
            graph[ilabel] = {}        
        for other in [j for j in xrange(nverts) if row[j]]:
            if noselfloops:
                if i == other: continue
            jlabel = other
            if labels is not None:            
                jlabel = labels[other]
            if jlabel not in graph:
                graph[jlabel] = {}
            graph[ilabel][jlabel] = adjmtx[i][other]
            if undirected:
                graph[jlabel][ilabel] = adjmtx[i][other]
    return graph        


def asincidence (graph):
    """ undirected graph -> incidence matrix.
    
    """
    n = len(graph)
    v = graph.keys()
    v.sort()
    e = list(edges(graph, undirected=True))
    e.sort()
    ne = len(e)
    I = zeroarray(n, ne)
    for k, edge in enumerate(e):
        i,j= edge
        i = v.index(i)
        j = v.index(j)     
        I[i][k] = 1
        I[j][k] = -1
    return I            
        

def asedges (graph, undirected=True):            
    """graph -> edge list, weights list.
    """
    edgelist = []    
    weights = []
    keys = graph.keys()
    keys.sort()
    
    for key in keys:
        if not len(graph[key]):
            edgelist.append( (key,) )   # singleton
            weights.append(None)
            continue
        for other in graph[key].keys():
            if undirected:
                if other > key:
                    edgelist.append( (key,other) )
                    weights.append(graph[key][other])
            else:
                edgelist.append( (key,other) )
                weights.append(graph[key][other])
                                
    return edgelist, weights       


def fromedges (edges, weights=None, undirected=True):
    """Create dict-of-dict graph representation from list of edges."""
    graph = {}    
    for i, edge in enumerate(edges):
        if len(edge) == 0:
            pass
        if len(edge) == 1:
            graph[edge[0]] = {}
        else:
            if weights is None:
                add(graph, edge[0], edge[1], weight=1, undirected=undirected)
            else:
                add(graph, edge[0], edge[1], weight=weights[i], undirected=undirected)
    return graph            


def adj2edgelist (adjmtx, undirected=True):
    """Create an edge list representation of a graph from an adj matrix.
    """
    edgelist, weights = [],[]
    if undirected:
        for i in range(len(adjmtx)-1):  
            nonzero = [x != 0 for x in adjmtx[i]]  
            if True not in nonzero: # singleton
                edgelist.append((i,))
                weights.append(None)
                continue            
            for j, val in enumerate(nonzero[i+1:]):
                if val is True:
                    edgelist.append((i, i+1+j))
                    weights.append(adjmtx[i][i+1+j])               
    else:
        for i in range(len(adjmtx)):  
            nonzero = [x != 0 for x in adjmtx[i]]  
            if True not in nonzero: # no out edges for this vertex
                edgelist.append((i,))
                weights.append(None)
                continue
            for j, val in enumerate(nonzero):
                if val is True:
                    edgelist.append((i,j))
                    weights.append(adjmtx[i][j])
    return edgelist, weights                            

#------------------------------------------------------------------------------#
# I/O functions

def tographviz (graph, filename, labels=None, partdict=None, undirected=True):
    """Saves the FOCI network in a file which can be read by Graphviz.
    
    labels is an optional list of vertex labels.
    
    Graphviz is a program for visualizing graphs.  It can be found on the web 
    at: 
    
        http://www.research.att.com/sw/tools/graphviz/
    
    """
    f = open(filename, 'w')
    
    if undirected:
        begin = 'Graph G{\n'
        sep = " -- "
    else:
        begin = 'digraph G{\n'
        sep = " -> "
    
    f.write(begin)

    keys = graph.keys()
    keys.sort()
    for key in keys:
        nbrs = neighbors(graph, key)
        
        if len(nbrs) == 0:
            f.write(str(key) + ";\n")
            continue 
                       
        for other in nbrs:
            attrstr = " [ weight=%s ] " % str(graph[key][other])
            if type(key) == type(""): keystr = '"' + key + '"'
            else: keystr = str(key)
            if type(other) == type(""): otherstr = '"' + other + '"'
            else: otherstr = str(other)
            if undirected:
                if other > key:
                    f.write(keystr + sep + otherstr + attrstr + ";\n")    
            else:
                f.write(keystr + sep + otherstr + attrstr + ";\n")

    if partdict is not None:
        partkeys = partdict.keys()
        partkeys.sort()
        for key in partkeys:
            f.write("subgraph cluster%s{\n" % str(key))
            for each in partdict[key]:
                f.write('\t"%s";\n' % str(each))
            f.write("}\n")            

    
    if labels is not None:
        for key in graph.keys():
            f.write(str(key) + ' [label="' + str(labels[key]) + '"];\n')
                        
    f.write('}')
    f.close()



def topajek (graph, filename, labeldict=None, undirected=True):     
    """Saves graph in file format which can be read by Pajek.
    
    Pajek is a program for visualizing large graphs.  It can be found on the web 
    at: 
    
        http://vlado.fmf.uni-lj.si/pub/networks/pajek/
    
    """                
    f = open(filename, 'w')
    verts = graph.keys()
    verts.sort()
    
    f.write('*Vertices  ' + str(len(graph)) + "\n" )
    if labeldict is not None:
        for i, vert in enumerate(verts):
            try:
                lbl = labeldict[vert]
            except KeyError:
                lbl = vert
            if len(lbl):
                f.write(str(i+1) + (' "' + str(lbl) + '"') + "\n")
            else:
                f.write(str(i+1) + (' "' + str(vert) + '"') + "\n")
    else:
        for i, vert in enumerate(verts):
            f.write(str(i+1) + (' "' + str(vert) + '"') + "\n")

    if undirected:
        f.write('*Edges\n')
    else:
        f.write('*Arcs\n')
    for edge in edges(graph,undirected=undirected):
        weight = graph[edge[0]][edge[1]]
        i = verts.index(edge[0])
        j = verts.index(edge[1])
        estring = "%d    %d    %02f" % (i+1, j+1, float(weight))
        f.write(estring + "\n")
    f.close()      


def topajekpart (graph, partdict, filename, nogroup = -99):
    """Saves a partition of a graph in format readable by Pajek.
    """
    verts = graph.keys()
    verts.sort()

    parts = partdict.keys()
    parts.sort()
    
    allint = True
    for key in parts:
        if type(key) != type(0):
            allint = False
            break
            
    v2p = {}
    for i,key in enumerate(parts):
        for v in partdict[key]:
            if allint:
                v2p[v] = key
            else:
                v2p[v] = i          
    
    f = open(filename, 'w')
    f.write('*Vertices ' + str(len(graph)) + "\n")
    for v in verts:
        try:
            f.write(str(v2p[v]) + "\n")
        except KeyError:        # the vertex is not part of a partition
            f.write("%d\n" % nogroup)
    f.close()


def readgraph (filename, isweighted=False, undirected=True, offset=0, commentchar="#"):
    """Parse a graph file, returns a dict-of-dict representation of a grpah.
    
    A graph file gives the edges of the graph. It's of the form:
    1   2   1
    1   3   2 ... etc
    
    Where the first two columns give vertex indices, and the third 
    column (optional) gives edge weights.  The file can have additional 
    columns (e.g. annotation), but this is currently ignored.    
    
    The optional argument 'offset' specifies an integer value used to offset
    the vertex indices.  This if useful for shifting index values.  For
    example, if the graph in the file is indexed from "1" but you'd rather
    have it index from "0" then set offset=-1.
    """
    f = open(filename, 'r')
    graph = {}
    for line in f.readlines():
        if line[0] == commentchar:
            continue
        words = line.split()
        if not len(words):
            continue
        v1 = int(words[0]) + offset
        try:
            v2 = int(words[1]) + offset
            graph.setdefault(v1, {})
            graph.setdefault(v2, {})
        except IndexError:
            graph.setdefault(v1, {})
            continue
        weight = 1
        if isweighted:
            try:
                weight = int(words[3])
            except IndexError:
                pass
        graph[v1][v2] = weight
        if undirected:
            graph[v2][v1] = weight
    return graph


def writegraph (graph, filename, undirected=True, withprogress=False):
    """Write adj. dict to a text file.      
    
    A graph file gives the edges of the graph. It's of the form:
    1   2   1
    1   3   1 ... etc
    
    The first two entries give vertex indices.
    The third entry of each line is the edge weight.
    """
    f = open(filename, 'w')
    ct = 0
    keys = graph.keys()
    keys.sort()
    for key in keys:
        s = ""
        if not len(graph[key]):
            s += "%s\n" % str(key)   # lone vertex
        okeys = graph[key].keys()
        okeys.sort()
        
        if not undirected:
            for other in okeys:
                s += "%s\t%s\t%s\n" % (str(key), str(other), str(graph[other][key]))  
        else:
            for other in okeys:
                if other > key:
                    s += "%s\t%s\t%s\n" % (str(key), str(other), str(graph[other][key]))  
                                    
        f.write(s)                
        ct += 1
        if withprogress:
            print "\b"*10 + str(ct),
            sys.stdout.flush()                  
    f.close()            



        

#------------------------------------------------------------------------------#
# Misc. routines

def aspartition(groups):
    """A sequence of groups -> partition dict.
    """
    partdict = {}
    for i,grp in enumerate(list(groups)):
        partdict[i] = Set(grp)
    return partdict                


def invertpartition(partdict):
    m2g = {}
    for grp in partdict:
        for mmbr in partdict[grp]:
            m2g[mmbr] = grp
    return m2g            

def partitionitems(partdict):
    """partdict -> Set containing all items that appear in a partition.
    """
    all = Set()
    for key in partdict:
        all |= Set(partdict[key])
    return all        
    
    
#------------------------------------------------------------------------------#

# Test graph

testgraph = {0: {1: 1, 2: 1, 3: 1, 5: 1}, 1: {0: 1, 2: 1, 3: 1}, 
2: {0: 1, 1: 1, 3: 1, 4: 1}, 3: {0: 1, 1: 1, 2: 1, 8: 1}, 4: {2: 1}, 
5: {0: 1, 10: 1, 12: 1, 6: 1}, 6: {8:1, 9: 1, 5: 1, 7: 1}, 7: {8: 1, 9: 1, 6: 1, 
15: 1}, 8: {9: 1, 3: 1, 13: 1, 6: 1, 7: 1}, 9: {8: 1, 14: 1, 13: 1, 6: 
1, 7: 1}, 10: {11: 1, 5: 1}, 11: {10: 1, 12:1}, 12: {11: 1, 5: 1}, 13: 
{8: 1, 9: 1}, 14: {9: 1, 15: 1}, 15: {14: 1, 7: 1}}


def test():
    pass            
        
#------------------------------------------------------------------------------#        


if __name__ == "__main__":
    test()

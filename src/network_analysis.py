import sys
sys.path.append("/mnt/chaelab/rachelle/src")

import csv
import itertools
import networkx as nx
import pandas as pd

# import basics

dir_root = "/mnt/chaelab/rachelle"
dir_proj = dir_root + "/zzOtherzz/XiaoMei/HDstress"


#####
## copied from basics
#####

def get_count_dict(iterable):
    output = {}
    for e in iterable:
        output[e] = output.get(e, 0) + 1
    return output

## returns list of sets
def merge_overlapping(*iters):
    last_len = len(iters)
    iters = [set(e) for e in iters]
    new_iters = []
    while True:
        while len(iters) > 0:
            s1 = iters.pop()
            i = 0
            while i < len(iters):
                if s1 & iters[i]:
                    s2 = iters.pop(i)
                    s1.update(s2)
                else:
                    i += 1
            new_iters.append(s1)
        if len(new_iters) == last_len:
            break
        else:
            last_len = len(new_iters)
            iters = new_iters
            new_iters = []
    return new_iters


############
##  network-analysis
############

## add nodes
def nodes_iter(fname):
    with open(fname, 'r') as f:
        header = f.readline() ## skip header
        header = header.rstrip().split('\t')
        for entry in f:
            entry_split = entry.rstrip().split('\t')
            metadata = {header[i]: entry_split[i] for i in range(len(header))}
            group, module = metadata["group_module"].split('_')
            metadata["group"] = group
            metadata["module"] = module
            ## entry_split[0] is gid; everything else is metadata
            yield (entry_split[0], metadata)
    return

## add edges
## for data from 2024
def edges_iter_v1(fname):
    with open(fname, 'r') as f:
        header = f.readline() ## skip header
        header = header.rstrip().split('\t')
        for entry in f:
            entry_split = entry.rstrip().split('\t')
            ## entry_split[2] is regulatoryGene, entry_split[3] is targetGene
            yield entry_split[2:4] + [{header[i]: entry_split[i] for i in range(len(header))}]
    return

## for 20250210 data
def edges_iter_v2(fname):
    with open(fname, 'r') as f:
        header = f.readline() ## skip header
        header = header.rstrip().split('\t')
        for entry in f:
            entry_split = entry.rstrip().split('\t')
            ## entry_split[2] is regulatoryGene, entry_split[3] is targetGene
            yield entry_split[:2] + [{header[i]: entry_split[i] for i in range(len(header))}]
    return

## function to apply nx.union to multiple graphs simultaneously
def multigraph_union(seed_graph, *graphs):
    if not graphs:
        return seed_graph
    for graph in graphs:
        seed_graph = nx.union(seed_graph, graph)
    return seed_graph

## stats class
class GraphStats():
    def __init__(self, G, parent_G, iteration = 0):
        self.iteration = iteration
        self.parent_G = parent_G
        self.G = G
        self.children = []
        self.group_counts = {}
        self._make_group_counts()
    def _make_group_counts(self):
        # valid_nodes = [node for node in self.G.nodes.values() if "group" in node]
        # groups = [node["group"] for node in valid_nodes]
        groups = [node["group"] for node in self.G.nodes.values()]
        # self.group_counts = basics.get_count_dict(groups)
        self.group_counts = get_count_dict(groups) # function copied to this script
    @property
    def size(self):
        return len(self.G)

## class to track GraphStats of successively pruned networks
class SuccessivelyPrunedGraphStats():
    def __init__(self, G, graphstats_class = GraphStats, filter_function = lambda G, nid:True):
        self.G = G
        self.graphstats_class = graphstats_class
        # [self.make_graphstats(sg, parent_G = sg, iteration = curr_iteration) for i, sg
        #                                 in enumerate(subgraphs) if i != i_subgraph]
        self.newest_subgraphs = [G.subgraph(c).copy() for c in nx.weakly_connected_components(G)]
        self.graphstats = {0: [self.make_graphstats(sg, G, 0) for sg in self.newest_subgraphs]} ## indexed by iteration
        self.filter_function = filter_function
        self.last_pruned_in_nodes = [] ## nodes with edges pointing FROM them to the last pruned node (i.e. in to pruned node)
        self.last_pruned_out_nodes = [] ## nodes with edges pointing TO them from the last pruned node (i.e. out of pruned node)
    def valid_prune_nid(self, nid):
        return nid in self.G.nodes and self.filter_function(self.G, nid)
    ## PRUNE FUNCTION!!! VERY IMPORTANT!!! Should be defined in chlild classes before use
    def next_prune_nid(self):
        pass
    def make_graphstats(self, G, *args, **kwargs):
        return self.graphstats_class(G, *args, **kwargs)
    @property
    def last_iteration(self):
        return max(self.graphstats)
    @property
    def curr_iteration(self):
        return self.last_iteration + 1
    def graphstats_at_iteration(self, i):
        return self.graphstats.get(i, [])
    def subgraphs_at_iteration(self, i):
        return [graphstats.G for graphstats in self.graphstats_at_iteration(i)]
    def graph_at_iteration(self, i):
        subgraphs = self.subgraphs_at_iteration(i)
        return multigraph_union(*subgraphs)
    def graphstats_list(self):
        return list(itertools.chain(*self.graphstats.values()))
    @property
    def last_subgraphs(self):
        return self.subgraphs_at_iteration(self.last_iteration)
    @property
    def last_graph(self):
        return self.graph_at_iteration(self.last_iteration)
    def set_graphstats_at_iteration(self, i, graphstats):
        self.graphstats[i] = list(graphstats)
    def prune(self):
        ## get next valid node id
        nid = self.next_prune_nid()
        ## exit if no valid nid
        if nid is None: return
        ## prep variables
        curr_iteration = self.curr_iteration
        curr_graphstats = []
        subgraphs = self.last_subgraphs
        ## iterate through subgraphs to find the one that needs to be split (i.e. the one containing nid)
        for i_subgraph, subgraph in enumerate(subgraphs):
            if nid in subgraph.nodes:
                ## update in and out nodes
                self.last_pruned_out_nodes = [edge[1] for edge in subgraph.out_edges(nid)]
                self.last_pruned_in_nodes = [edge[0] for edge in subgraph.in_edges(nid)]
                ## create a copy of the subgraph and split that
                subgraph_copy = subgraph.copy()
                subgraph_copy.remove_node(nid)
                ## get resultant subgraphs
                subsubgraphs = [subgraph_copy.subgraph(c).copy()
                                for c in nx.weakly_connected_components(subgraph_copy)]
                ## get stats for unmodified subgraphs
                curr_graphstats.extend([self.make_graphstats(sg, parent_G = sg, iteration = curr_iteration)
                                        for i, sg in enumerate(subgraphs) if i != i_subgraph])
                ## get stats for newly split subgraphs
                curr_graphstats.extend([self.make_graphstats(ssg, parent_G = subgraph,
                                                             iteration = curr_iteration)
                                        for i, ssg in enumerate(subsubgraphs)])
                ## update with stats for current iteration and break
                self.set_graphstats_at_iteration(curr_iteration, curr_graphstats)
                self.newest_subgraphs = subsubgraphs
                break
        return


## prune nodes in order of decreasing out-degree (determined by initial graph)
## sort_function takes G and outputs nodes in sorted list in order of pruning
class SuccessivelyPrunedGraphStats_SortAtInit(SuccessivelyPrunedGraphStats):
    def __init__(self, G, sort_function, *args, **kwargs):
        super().__init__(G, *args, **kwargs)
        self.sort_function = sort_function ## not really used beyond __init__, but we'll keep it just in case
        self.prune_order = [nid for nid in sort_function(G) if self.valid_prune_nid(nid)]
        self.next_prune_id = 0
    def next_prune_nid(self):
        while self.next_prune_id < len(self.prune_order):
            self.next_prune_id += 1
            return self.prune_order[self.next_prune_id-1]
        return None

## prune nodes with greatest out-degree (re-calculated at each iteration)
## okay maybe this class does do something ## this class doesn't actually do anything lol it's just here so we have a corresponding class for dynamically determining next node to prune that is same level as SuccessivelyPrunedGraphStats_SortAtInit
class SuccessivelyPrunedGraphStats_SortDynamic(SuccessivelyPrunedGraphStats):
    def __init__(self, G, *args, **kwargs):
        super().__init__(G, *args, **kwargs)
        self.valid_prune_nids = set(nid for nid in self.G.nodes if self.valid_prune_nid(nid))
    def remove_valid_prune_nid(self, nid):
        self.valid_prune_nids.remove(nid)

## prune nodes in order of decreasing out-degree (determined by initial graph)
class SuccessivelyPrunedGraphStats_MaxOutAtInit(SuccessivelyPrunedGraphStats_SortAtInit):
    def __init__(self, G, *args, **kwargs):
        super().__init__(G, (lambda G:sorted(G.nodes, key = lambda nid:G.out_degree(nid), reverse = True)), ## sort_function
                         *args, **kwargs)

class SuccessivelyPrunedGraphStats_MaxOutAtInit2TF(SuccessivelyPrunedGraphStats_SortAtInit):
    def __init__(self, G, *args, **kwargs):
        super().__init__(
            G,
            (lambda G:sorted(G.nodes,
                             key = lambda nid:len([e for e in G.out_edges(nid) if G.nodes[e[1]]["isTF"] == '1']),
                             reverse = True)), ## sort_function
            *args, **kwargs)

## prune nodes with greatest out-degree (re-calculated at each iteration)
class SuccessivelyPrunedGraphStats_MaxOutDynamic(SuccessivelyPrunedGraphStats_SortDynamic):
    def __init__(self, G, *args, **kwargs):
        super().__init__(G, *args, **kwargs)
        self.last_pruned_in_nodes = set(G.nodes)
        self.out_degrees = {nid: G.out_degree(nid) for nid in self.valid_prune_nids}
    def update_out_degrees(self):
        for nid in set(self.last_pruned_in_nodes).intersection(self.valid_prune_nids):
            self.out_degrees[nid] -= 1
        return
    def next_prune_nid(self):
        self.update_out_degrees()
        if not self.valid_prune_nids: return None
        nid_to_prune = max(self.valid_prune_nids, key = lambda nid:self.out_degrees.get(nid,0))
        self.remove_valid_prune_nid(nid_to_prune)
        return nid_to_prune

## prune nodes with greatest out-degree (re-calculated at each iteration)
class SuccessivelyPrunedGraphStats_MaxOutDynamic2TF(SuccessivelyPrunedGraphStats_MaxOutDynamic):
    def __init__(self, G, *args, **kwargs):
        super().__init__(G, *args, **kwargs)
        self.out_degrees = {nid: len([e for e in G.out_edges(nid) if G.nodes[e[1]]["isTF"] == '1'])
                            for nid in self.valid_prune_nids}

## idk why I have to do this but just plugging the lambda into SuccessivelyPruneGraphStats_MaxOutAtInit's init statement causes them all to share the last group :/
def make_group_restricted_filter(group):
    grp = group ## this line is crucial for stopping lambda from referencing
    return lambda G,x:G.nodes[x]["group"] == grp

## generates separate SPGS objects for pruning nodes of different groups
def make_pruneG_multigroup(G, SPGS_class, groups, include_all = True,
                           and_filter = lambda G,x:True, misc_args = [], misc_kwargs = {}):
    pruneG_multi = {}
    for group in groups:
        group_filter = make_group_restricted_filter(group)
        def new_filter_function(G,x):
            return group_filter(G,x) and and_filter(G,x)
        pruneG_multi[group] = SPGS_class(G, *misc_args, filter_function = new_filter_function, **misc_kwargs)
    pruneG_multi["all"] = SPGS_class(G, *misc_args, filter_function = lambda G, x: True, **misc_kwargs)
    return pruneG_multi

def execute_pruneG_multigroup(pruneG_multigroup, maxdepth = 50, print_iteration = False):
    if print_iteration: print_i = lambda i:print(i)
    else: print_i = lambda i:None
    ## start pruning
    for i in range(maxdepth):
        print_i(i)
        for group, SPGS in pruneG_multigroup.items():
            SPGS.prune()
    return

def write_pruneG_multigroup_v1(pruneG_multigroup, fout):
    with open(fout, "w+") as f:
        _ = f.write('\t'.join(["prune_grp", "iteration", "graph_id", "parent_graph_id", "size"] +
                              [f"size_{grp}" for grp in groups]) + '\n')
        for group, group_stats in pruneG_multigroup.items():
            for entry in group_stats.graphstats_list():
                _ = f.write('\t'.join(map(str,[group, entry.iteration, id(entry.G),
                                               id(entry.parent_G), entry.size] +
                                          [entry.group_counts.get(grp, 0) for grp in groups])) + '\n')
    return

## d_pruneG_multigroups: should be a dictionary of {"<colname for prune algo>": <pruneG_multigroup dict>}
##  where <pruneG_multigroup dict> is a dictionary of {"<group pruned>": <SuccessivelyPrunedGraphStats obj>}
## iterations: iterable of prune iterations to write subnetwork membership information for
## encoding:
##  0 for nodes present in largest subnetwork at given iteration.
##  1 for nodes present in next largest subnetwork, so on and so forth.
##  -1 for nodes not in network at all (i.e. pruned)
def write_pruneG_multigroup_membershipatiter_v1(d_pruneG_multigroups, fout, iterations):
    nid_all = set()
    for pruneG_multigroup in d_pruneG_multigroups.values():
        for pruneG in pruneG_multigroup.values():
            nid_all.update(pruneG.G.nodes)
    nid_all = list(nid_all)
    df = pd.DataFrame({"nid": nid_all})
    for prune_algo, pruneG_multigroup in d_pruneG_multigroups.items():
        for group, group_stats in pruneG_multigroup.items():
            for iteration in iterations:
                if iteration not in group_stats.graphstats: continue
                iter_graphstats = group_stats.graphstats[iteration]
                sorted_graphstats = sorted(iter_graphstats,
                                           key = lambda graphstats:len(graphstats.G),
                                           reverse = True)
                membership = {nid: grp
                              for grp, graphstats in enumerate(sorted_graphstats)
                              for nid in graphstats.G.nodes}
                df_tmp = pd.DataFrame(
                    {"nid": nid_all,
                     f"{prune_algo}_{group}_{iteration}": [membership.get(nid, -1) for nid in nid_all]}
                )
                df = df.join(df_tmp.set_index("nid"), on = "nid", how = "left")
    df.to_csv(fout, index = False, header = True, sep = '\t', quoting = csv.QUOTE_NONE)
    return

# edges_iter = edges_iter_v1
edges_iter = edges_iter_v2

def make_graph_from_edges_and_nodes_file(f_edges, f_nodes):
    G = nx.MultiDiGraph()
    G.add_nodes_from(nodes_iter(f_nodes))
    keys = G.add_edges_from(edges_iter(f_edges))
    return G


## shared variables
groups = ("ECR", "C", "HD", "CHD")

## random node id in GRN that we can use to troubleshoot 'BraA03g055090.3.5C'

## parse & prune networks
network_names = ("GRN_simplify", "keyGRN_simplify", "TFTF_simplify")
# network_name = "GRN_simplify"
print_iteration = True
include_all_group = True
write_pruneG_multigroup = write_pruneG_multigroup_v1
write_pruneG_multigroup_membershipatiter = write_pruneG_multigroup_membershipatiter_v1
maxdepth = 100
write_membership_iters = [25, 50]
for network_name in network_names:
    print(network_name)
    ## input
    f_edges = dir_proj + f"/data/network/{network_name}.edges"
    f_nodes = dir_proj + f"/data/network/{network_name}.nodes"
    ## make graph
    G = make_graph_from_edges_and_nodes_file(f_edges, f_nodes)
    ########### PRUNE ORDER DETERMINED BY INITIAL GRAPH ##################
    ## prune nodes in order of decreasing out-degree (determined by initial graph)
    pruneG_serial_maxinit = make_pruneG_multigroup(G, SuccessivelyPrunedGraphStats_MaxOutAtInit,
                                                   groups, include_all = include_all_group)
    execute_pruneG_multigroup(pruneG_serial_maxinit, maxdepth = maxdepth, print_iteration = print_iteration)
    write_pruneG_multigroup(pruneG_serial_maxinit,
                            fout = dir_proj + f"/results/network/{network_name}.prune.maxOutAtInit.tsv")
    ## prune nodes in order of decreasing out-degree to TFs (determined by initial graph)
    pruneG_serial_maxinit2TF = make_pruneG_multigroup(G, SuccessivelyPrunedGraphStats_MaxOutAtInit2TF,
                                                      groups, include_all = include_all_group)
    execute_pruneG_multigroup(pruneG_serial_maxinit2TF, maxdepth = maxdepth, print_iteration = print_iteration)
    write_pruneG_multigroup(pruneG_serial_maxinit2TF,
                            fout = dir_proj + f"/results/network/{network_name}.prune.maxOutAtInit2TF.tsv")
    ########### PRUNE ORDER RE-DETERMINED AT EACH ITERATION ##############
    ## prune nodes with max out-degree at each iteration
    pruneG_serial_maxdynamic = make_pruneG_multigroup(G, SuccessivelyPrunedGraphStats_MaxOutDynamic,
                                                      groups, include_all = include_all_group)
    execute_pruneG_multigroup(pruneG_serial_maxdynamic, maxdepth = maxdepth, print_iteration = print_iteration)
    write_pruneG_multigroup(pruneG_serial_maxdynamic,
                            fout = dir_proj + f"/results/network/{network_name}.prune.maxOutDynamic.tsv")
    ## prune nodes with max out-degree to TFs at each iteration
    pruneG_serial_maxdynamic2TF = make_pruneG_multigroup(G, SuccessivelyPrunedGraphStats_MaxOutDynamic2TF,
                                                         groups, include_all = include_all_group)
    execute_pruneG_multigroup(pruneG_serial_maxdynamic2TF, maxdepth = maxdepth,
                              print_iteration = print_iteration)
    write_pruneG_multigroup(pruneG_serial_maxdynamic2TF,
                            fout = dir_proj + f"/results/network/{network_name}.prune.maxOutDynamic2TF.tsv")
    ## write subnetwork membership at specific iterations
    write_pruneG_multigroup_membershipatiter(
        {"maxOutAtInit": pruneG_serial_maxinit,
         "maxOutAtInit2TF": pruneG_serial_maxinit2TF,
         "maxOutDynamic": pruneG_serial_maxdynamic,
         "maxOutDynamic2TF": pruneG_serial_maxdynamic2TF},
        fout = dir_proj + f"/results/network/{network_name}.prune.membershipAtIter.tsv",
        iterations = write_membership_iters)


## find master controller of cHD (and maybe +HD+C too?) nodes

class Node():
    ## GNT: GraphNodeTracker object (wraps MultiDiGraph)
    def __init__(self, GNT, nid):
        self.CN = None ## CyclicNode object; else None if node is not part of a cylic graph
        self.GNT = GNT
        self.G = GNT.G
        self.id = nid
        self.children = [] ## direct children
        self.parents = [] ## direct parents
        self.cache = {}
    def set_cyclic_node(self, CN):
        self.CN = CN
        for node in self.children:
            if self in node.parents:
                node.parents.remove(self)
            if node is not CN and CN not in node.parents:
                node.parents.append(CN)
        for node in self.parents:
            if self in node.children:
                node.children.remove(self)
            if node is not CN and CN not in node.children:
                node.children.append(CN)
        return
    def node(self):
        return self if self.CN is None else self.CN
    def raw_node(self):
        return self
    def is_terminal(self):
        return len(self.children) == 0
    def all_descendants(self):
        return self.children + list(itertools.chain(*[node.all_descendants() for node in self.children]))
    def internal_children(self):
        return [node for node in self.children if not node.is_terminal()]
    def internal_descendants(self):
        return self.internal_children() + list(itertools.chain(*[node.internal_descendants() for node in self.children]))
    def terminal_children(self):
        return [node for node in self.children if node.is_terminal()]
    def terminal_descendants(self):
        return self.terminal_children() + list(itertools.chain(*[node.terminal_descendants() for node in self.children]))
    def get_property(self, property_name, unlist = False): ## returns list by default; this is so it is compatible with CyclicNode
        result = self.G.nodes[self.id].get(property_name, None)
        if unlist: return result
        return [result]
    def _link_nodes(self):
        self.parents = [self.GNT.get_node(in_edge[0]) for in_edge in self.G.in_edges(self.id)]
        self.children = [self.GNT.get_node(out_edge[1]) for out_edge in self.G.out_edges(self.id)]
        return
    def set_cache(self, varname, value):
        self.cache[varname] = value
    def set_terminals_bottom_up(self, varname, f):
        for node in self.internal_children():
            node.set_terminals_bottom_up(varname, f)
        for node in self.terminal_children():
            node.set_cache(varname, f(node))
        return
    def set_terminals_top_down(self, varname, f):
        for node in self.terminal_children():
            node.set_cache(varname, f(node))
        for node in self.internal_children():
            node.set_terminals_bottom_up(varname, f)
        return
    def set_internals_bottom_up(self, varname, f):
        if self.is_terminal(): return
        for node in self.internal_children():
            node.set_internals_bottom_up(varname, f)
        self.set_cache(varname, f(self))
        return
    def set_terminals_top_down(self, varname, f):
        if self.is_terminal(): return
        self.set_cache(varname, f(self))
        for node in self.internal_children():
            node.set_internals_bottom_up(varname, f)
        return

## this is a pseudo-node, which will replace Node for all traversals
class CyclicNode(Node):
    def __init__(self, GNT, cid, nodes):
        super().__init__(GNT, cid)
        self.nodes = nodes
        self.nids = [node.id for node in self.nodes]
        self.children = []
        self.parents = []
        self._set_children()
        self._set_parents()
    def node(self): return self
    def get_property(self, property_name, **kwargs): ## returns list (kwargs is to sink unlist)
        return list(itertools.chain(*[node.get_property(property_name) for node in self.nodes]))
    def _set_children(self):
        self.children = [node for node in set(itertools.chain(*[node.children for node in self.nodes]))
                         if node not in self.nodes and node is not self] ## children should not contain nodes in cycle
        ## propagate to child nodes
        for node in self.children:
            node.parents = [node for node in node.parents if node not in self.nodes]
            node.parents.append(self)
    def _set_parents(self):
        self.parents = [node for node in set(itertools.chain(*[node.parents for node in self.nodes]))
                        if node not in self.nodes and node is not self] ## parents should not contain nodes in cycle
        ## propagate to parent nodes
        for node in self.parents:
            node.children = [node for node in node.children if node not in self.nodes]
            node.children.append(self)

## this object basically tracks nodes in a directed graph like they're in a tree
## merges nodes in cyclic graph into CyclicNode pseudo-node)
class GraphNodeTracker():
    def __init__(self, G):
        self.G = G ## DO NOT MUTATE THIS GRAPH; properties of GraphNodeTracker & its Nodes will not be updated
        self._raw_nodes = {nid: Node(self, nid) for nid in G.nodes} ## does not contain CyclicNode objects
        self._cyclic_nodes = {}
        self._nodes = {} ## replaces relevant Node objects with CyclicNode objects
        self.root_nodes = []
        self.terminal_nodes = []
        self._set_nodes()
    @property
    def nodes(self):
        return self._nodes.values()
    @property
    def raw_nodes(self):
        return self._raw_nodes.values()
    def get_raw_node(self, nid): ## does not map to CyclicObject
        return self._raw_nodes[nid]
    def get_node(self, nid, raw = False): ## will map nodes in cyclic graphs to relevant CyclicNode object & return that instead
        raw_node = self._nodes.get(nid, None)
        if raw_node is None: raw_node = self.get_raw_node(nid)
        if raw: return raw_node
        return raw_node.node() ## returns CyclicNode if current node is part of a cycle
    def _set_nodes(self):
        for node in self.raw_nodes:
            node._link_nodes()
        self._set_cyclic_nodes()
        self._set_root_nodes()
        self._set_terminal_nodes()
        return
    def _set_root_nodes(self):
        self.root_nodes = [node for node in self.nodes if len(node.parents) == 0]
    def _set_terminal_nodes(self):
        self.terminal_nodes = [node for node in self.nodes if len(node.children) == 0]
    def _set_cyclic_nodes(self):
        cyclic_subgraphs = nx.simple_cycles(self.G)
        ## merge overlapping cyclic graphs
        # merged_subgraphs = basics.merge_overlapping(*cyclic_subgraphs)
        merged_subgraphs = merge_overlapping(*cyclic_subgraphs) ## function copied to this script
        ## duplicate _raw_nodes into _nodes
        self._nodes = {nid: node for nid, node in self._raw_nodes.items()}
        ## make CyclicNode objects & insert into _nodes
        for i, subgraph in enumerate(merged_subgraphs):
            cid = f"PseudoNode_{i}"
            new_CN = CyclicNode(self, cid, [self.get_raw_node(nid) for nid in subgraph])
            self._cyclic_nodes[cid] = new_CN
            ## remove Node from subgraphs, parents, and children and insert CyclicNode
            for nid in subgraph:
                node = self.get_raw_node(nid)
                node.set_cyclic_node(new_CN) ## removes Node from children and parent & replaces with CyclicNode (also, sets Node.CN)
                del self._nodes[nid]
            self._nodes[cid] = new_CN
        return


## function to set how many cHD terminal nodes are pointed to directly by each internal node
def make_get_children_property_value_counts(property_name, value):
    def get_children_property_value_counts(node):
        return sum(sum(v == value for v in child_node.get_property(property_name))
                   for child_node in node.children)
    return get_children_property_value_counts

def make_tally_descendant_cache(property_name, tallied_property_name):
    def tally_descendant_cache(node):
        return node.cache.get(property_name, 0) + sum(child_node.cache.get(tallied_property_name, 0)
                                                      for child_node in node.children)
    return tally_descendant_cache

def write_GNT_cache_v1(gnt, fout, names):
    with open(fout, "w+") as f:
        header = ["id", "nid", "is_pseudonode"] + list(names)
        _ = f.write('\t'.join(header) + '\n')
        for node in gnt.nodes:
            _ = f.write('\t'.join(map(
                str,
                ([node.id, node.id if not isinstance(node, CyclicNode) else ','.join(node.nids), isinstance(node, CyclicNode)] +
                 [node.cache.get(name, None) for name in names]))) + '\n')
    return

## parse & prune networks
network_names = ("GRN_simplify", "keyGRN_simplify", "TFTF_simplify")
# network_name = "GRN_simplify"
write_gnt = write_GNT_cache_v1
for network_name in network_names:
    print(network_name)
    ## input
    f_edges = dir_proj + f"/data/network/{network_name}.edges"
    f_nodes = dir_proj + f"/data/network/{network_name}.nodes"
    ## make graph
    G = make_graph_from_edges_and_nodes_file(f_edges, f_nodes)
    ## make object to track nodes
    gnt = GraphNodeTracker(G)
    for node in gnt.nodes:
        node.set_cache("group", ','.join(node.get_property("group", unlist = False)))
    for node in gnt.root_nodes:
        node.set_internals_bottom_up("descendants", lambda node:len(tuple(node.all_descendants())))
        for group in groups:
            node.set_internals_bottom_up(f"child_{group}", make_get_children_property_value_counts("group", group))
            node.set_internals_bottom_up(f"descendant_{group}", make_tally_descendant_cache(f"child_{group}", f"descendant_{group}"))
            node.set_internals_bottom_up(f"descendant_fraction_{group}",
                                         lambda node:node.cache.get(f"descendant_{group}", 0)/node.cache["descendants"])
    ## write
    write_gnt(gnt, dir_proj + f"/results/network/{network_name}.descendantStats.tsv",
              (["group", "descendants"] +
               list(itertools.chain(*[[f"child_{group}", f"descendant_{group}", f"descendant_fraction_{group}"] for group in groups]))))

# nodes = {Node(G, nid) for nid in G.nodes}

# nodes_with_outdegree = set(nid for nid in G.nodes if G.out_degree(nid) > 0)
# terminal_out_nodes = [nid for nid in nodes_with_outdegree if sum(edge[1] in nodes_with_outdegree for edge in G.out_edges(nid)) == 0]

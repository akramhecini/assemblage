import networkx as nx
import random
import statistics


def read_fastq(fasta_q):
  with open(fasta_q, 'r') as f:
    for count, line in enumerate(f, start=1):
      if (count/2)%2 == 1: #selectionner seulement les lignes du style 2,6,10 et pas 4 et 8
        yield (line.strip("\n"))


def cut_kmer(seq,k_mer):
    for i in range(len(seq)-k_mer+1):
        yield seq[i:i+k_mer]


def build_kmer_dict(fasta_q,k_mer): # finalement il renvoie ça que pour la derniere seq ! faut t'il conctner les seq?
        dic = {}

        it = read_fastq(fasta_q)

        for sequence in it:
            it_tmp = cut_kmer(sequence,k_mer) #renvoie un itér des seq de taille k_mer
            list_tmp_kmer = list(it_tmp)
            for sous_k_mer in list_tmp_kmer:

                if sous_k_mer in dic:
                    dic[sous_k_mer] += 1
                else:
                      dic[sous_k_mer] = 1


        return dic



def build_graph(dic_kmer) :

    Grph = nx.DiGraph()

    for kmer, valeur in dic_kmer.items():

        nd1 = kmer[:-1]

        nd2 = kmer[1:]

        Grph.add_edge(nd1 , nd2 , weight = valeur)

    return Grph


def get_starting_nodes(graph):

    node_entry = []

    for node in graph.nodes:

        if len(list(graph.predecessors(node))) == 0:

            node_entry.append(node)

    return node_entry


def get_sink_nodes(graph):

    node_out = []

    for node in graph.nodes:

        if len(list(graph.successors(node))) == 0:

            node_out.append(node)

    return node_out



def get_contigs(graph, entree, sortie):

    paths = []

    for nd in entree:

        for nd2 in sortie:

            path = list(nx.all_simple_paths(graph, nd, nd2))

            longu = len(path)

            if longu > 0:

                paths.append(path)

    reslt = ()

    for i in range(len(paths)):
        for ii in range(len(paths[i])):

            contig = str(paths[i][ii][0])

            for j in range(1, len(paths[i][ii])):

                tmp = str(paths[i][ii][j][-1])
                contig += tmp[-1]

            reslt.append([contig, len(contig)])

    return reslt




def fill(text, width=80):
    """Split text with a line returnto respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))


def std(liste):
    """ a function that takes a list of values and return its stdev"""

    res = statistics.stdev(liste)
    return res



def path_average_weight(graph, path):

    """ a function that take a graph and a path
    and it returns the average weigth of a sub graohe"""

    sous_graphe = graph.subgraph(path)

    liste_w = []

    for edge in sous_graphe.edges(data=True):
        liste_w.append(edge[2]['weight'])

    return statistics.mean(liste_w)


def remove_paths(graph, liste_path, delete_entry_node, delete_sink_node):

    """ Une fonction qui prend un graphe et une liste de chemin, delete_entry_node pour
indiquer si les noeuds d’entrée seront supprimés et delete_sink_node pour indiquer si
les noeuds de sortie seront supprimés et retourne un graphe nettoyé des chemins
indésirables."""

    graphe_tmp = graph

    for i in range(len(liste_path)):

        graphe_tmp.remove_nodes_from(liste_path[i][1:-1])

        if delete_entry_node == True:
            graphe_tmp.remove_node(liste_path[i][0]) #supprimer le noeud d entree
        if delete_sink_node == True:
            graphe_tmp.remove_node(liste_path[i][-1])#supprimer le noeud de sortie

    return graphe_tmp

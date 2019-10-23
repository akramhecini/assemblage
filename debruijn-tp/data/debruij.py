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

    for kmer, valeur in enumerate(dic_kmer.items()):

        nd1 = kmer[:-1]

        nd2 = kmer[1:]


        Grph.add_edge(nd1 , nd2 , weight = valeur)

    return Grph


def get_starting_nodes(graph):
    
    node_entry = []

    for node in graph.nodes:

        if len(list(graph.predecessors(node))) == 0:

            nody_entry.append(node)

    return node_entry


def get_sink_nodes(graph):
  
    node_out = []

    for node in graph.nodes:

        if len(list(graph.successors(node))) == 0:

            node_out.append(node)

    return sortis



def get_contigs(graph, entres, sortis):
    """Retourne liste de tuple(contig, taille du contig)"""
    
    paths = []

    for nd in entres:

        for nd2 in sortis:

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
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))


def save_contigs(tupple, file_name):
    """Sort un fichier avec la liste tupple"""
    file = open(file_name, 'w+')
    for i in range(len(tupple)):
        file.write('>contig_' + str(i) + ' len=' + str(tupple[i][1]) + '\n' + str(fill(tupple[i][0])) + '\n')
    file.close()


def std(list_val):
    """Take list of values and return standard deviation"""
    return st.stdev(list_val)


def path_average_weight(graph, path):
    """Take a graph and a path and return average weigth"""
    new_G = graph.subgraph(path)
    wei = []
    for arretes in new_G.edges(data=True):
        wei.append(arretes[2]['weight'])
    mean_wei = st.mean(wei)
    return mean_wei


def remove_paths(graph, path, delete_entry_node, delete_sink_node):
    """Take graph and path and remove nodes
    delete_entry_node pour indiquer si les noeuds d’entrée seront supprimés
    et delete_sink_node pour indiquer si les noeuds de sortie seront supprimés
    Return clean graph"""
    new_G = graph
    for i in range(len(path)):
        new_G.remove_nodes_from(path[i][1:-1])
        if delete_entry_node == True:
            new_G.remove_node(path[i][0])
        if delete_sink_node == True:
            new_G.remove_node(path[i][-1])
    return new_G


def select_best_path(graph, path, path_len, path_wei,
                     delete_entry_node=False, delete_sink_node=False):
    """Take graph, path, path length, path weigth and delete entry or sink node
    Function return clean graph
    By default delete entry or sink node is set as False"""
    new_G = graph
    random.seed(9001)
    best_wei_indice = []
    for i,e in enumerate(path_wei):
        if e == max(path_wei):
            best_wei_indice.append(i)
    best_len_indice = []
    for i,e in enumerate(path_len):
        if e == max(path_len):
            best_len_indice.append(i)
    if len(best_wei_indice) == 1:
        path.remove(path[best_wei_indice[0]])
        new_G = remove_paths(new_G, path, delete_entry_node, delete_sink_node)
    elif len(best_len_indice) == 1:
        path.remove(path[best_len_indice[0]])
        new_G = remove_paths(new_G, path, delete_entry_node, delete_sink_node)
    else:
        rand = random.randint(0,len(path)-1)
        path.remove(path[rand])
        new_G = remove_paths(new_G, path, delete_entry_node, delete_sink_node)
    return new_G
    

def solve_bubble(graph, ancestor_node, descendant_node):
    """Take graph, ancestor node, descendant node and return clean graph"""
    new_G = graph
    all_path = nx.all_simple_paths(graph, ancestor_node, descendant_node)
    all_path = list(all_path)
    weigth = []
    for i in range(len(list(all_path))):
        weigth.append(path_average_weight(new_G, all_path[i]))
    length = []
    for i in range(len(list(all_path))):
        length.append(len(all_path[i]))
    new_G = select_best_path(new_G, all_path, length, weigth)
    return new_G


def simplify_bubbles(graph):
    """Take graph and return it without bubble"""
    new_G = graph
    all_path = []
    bubble = []
    for node in new_G.nodes():
        pred = list(graph.predecessors(node))
        if len(pred) > 1:
            anc = nx.lowest_common_ancestor(new_G,pred[0],pred[1])
            bubble.append([anc,node])
    for i in range(len(bubble)):
        new_G = solve_bubble(new_G, bubble[i][0], bubble[i][1])
    return new_G


def solve_entry_tips(graph, entry_node):
    """Take graph and entry nodes and return graph whithout entry indesirable path"""

'''graph_1 = nx.DiGraph()
graph_1.add_weighted_edges_from([(1, 2, 10), (3, 2, 2), (2, 4, 15), (4, 5, 15)])
graph_1 = solve_entry_tips(graph_1, [1, 3])  
assert (3, 2) not in graph_1.edges()
assert (1, 2) in graph_1.edges()
graph_2 = nx.DiGraph()
graph_2.add_weighted_edges_from([(1, 2, 2), (6, 3, 2), (3, 2, 2),
                                     (2, 4, 15), (4, 5, 15)])
graph_2 = solve_entry_tips(graph_2, [1, 6])  
assert (1, 2) not in graph_2.edges()
assert (6, 3) in graph_2.edges()
assert (3, 2) in graph_2.edges()'''


def solve_out_tips(graph, sink_node):
    """Take graph and sink nodes and return graph whithout entry indesirable path"""



'''graph_1 = nx.DiGraph()
graph_1.add_weighted_edges_from([(1, 2, 15), (2, 3, 15), (4, 5, 15), (4, 6, 2)])
graph_1 = solve_out_tips(graph_1, [5, 6])  
assert (4, 6) not in graph_1.edges()
assert (4, 5) in graph_1.edges()  
graph_2 = nx.DiGraph()
graph_2.add_weighted_edges_from([(1, 2, 15), (2, 3, 15), (4, 5, 2), (4, 6, 2) , (6, 7, 2)])
graph_2 = solve_out_tips(graph_2, [5, 7])  
assert (4, 5) not in graph_2.edges()
assert (6, 7) in graph_2.edges()'''


def main():
    """Main function"""
    parser = argparse.ArgumentParser(description='Graphe de de Bruijn')
    parser.add_argument('-i', metavar='FASTQ', type=str, help='File FASTQ', required=True)
    parser.add_argument('-k', metavar='Kmer', type=int, help='Kmer size', default='21')
    #parser.add_argument('-o', metavar='Confiq', type=str, help='Config file', required=True)
    args = parser.parse_args()
    fastq_file = args.i
    kmer_size = args.k
    #config_file = args.o
    
    graph = build_graph(build_kmer_dict(fastq_file, kmer_size))
    tupple = get_contigs(graph, get_starting_nodes(graph), get_sink_nodes(graph))
    save_contigs(tupple, 'eva71_hundred_reads.fq.out')


if __name__ == '__main__':
    main()



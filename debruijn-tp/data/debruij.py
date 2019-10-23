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

def read_fastq(fasta_q):
  with open(fasta_q, 'r') as f:
    for count, line in enumerate(f, start=1):
      if (count/2)%2 == 1: #selectionner seulement les lignes du style 2,6,10 et pas 4 et 8
        yield (line.strip("\n"))


def cut_kmer(seq,k_mer):
    for i in range(len(seq)-k_mer):
        yield seq[i:i+k_mer+1]


def build_kmer_dict(fasta_q,k_mer): # finalement il renvoie ça que pour la derniere seq ! faut t'il conctner les seq?
        dic = {}
        it = read_fastq(fasta_q)

        #seq_concat = "".join(list(it)) #une seule séquence !
        #cut_kmer(seq_concat,k_mer)
        #occu_tmp = len(list(it_tmp)) #je récuère la longueur de mon it
        #dic[k_mer] = occu_tmp #j'assigne la longueur à k_mer dans mon dic

        for sequence in it:
            it_tmp = cut_kmer(sequence,k_mer) #renvoie un itér des seq de taille k_mer
            list_tmp_kmer = list(it_tmp)
            for sous_k_mer in list_tmp_kmer:
                cnt_tmp = list_tmp_kmer.count(sous_k_mer)
                if sous_k_mer in dic:
                    dic[sous_k_mer] += cnt_tmp
                else:
                      dic[sous_k_mer] = cnt_tmp

            #occu_tmp = len(list(it_tmp)) #je récuère la longueur de mon it
            #dic[k_mer] = occu_tmp #j'assigne la longueur à k_mer dans mon dic
        return dic



def build_graph(dic_kmer) :

    Grph = nx.DiGraph()

    for i , (kmer, nb) in enumerate(dic_kmer.items()):

        node1 = kmer[:-1]

        node2 = kmer[1:]


        Grph.add_edge(node1 , node2 , weight = nb)

    return Grph

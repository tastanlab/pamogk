import numpy as np


def mark_label_on_pathways(name, pid, pw_map, gene_id_list, label=1):
    """Marks given genes to the pathways

    Parameters
    ----------
    name: str
    pid: int
        patient id
    pw_map: map of networkx graphs of pathways
        patient label mapping
    gene_id_list: list of list of string
        uniprot gene id list of genes
    label: int
        the label which will be assigned to found genes in pathways - default value is 1
    """
    label_field = 'label-{}'.format(name)
    gene_ids = [uid for a in gene_id_list for uid in a]
    for pw in pw_map.values():  # for each pathway
        for n in pw.nodes():
            nd = pw.nodes[n]
            if label_field not in nd: pw.add_node(n, **{label_field: {}})
            if np.any([g in nd['uniprot-ids'] for g in gene_ids]):
                nd[label_field][pid] = label

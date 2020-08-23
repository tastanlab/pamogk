import matplotlib.pyplot as plt
import numpy as np


def draw_hist(pat_dict):
    count_mutated_genes = []
    for patient in pat_dict.keys():
        count_mutated_genes.append(len(pat_dict[patient]))
    count_mutated_genes = sorted(count_mutated_genes)
    np.savetxt('mutated_genes.txt', np.array(count_mutated_genes))
    plt.hist(count_mutated_genes, bins=max(count_mutated_genes), edgecolor='black',
             linewidth=0.2)  # arguments are passed to np.histogram
    plt.xlabel('Mutasyona Uğramış Gen Sayısı')
    plt.ylabel('Hasta Sayısı')
    plt.savefig('a.png')
    plt.show()

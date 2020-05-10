"""
for pw in all_pathways:
    for n in pw:
        rw = random walks x RW_CNT
        vis_freq_map = freq_count_of_genes(rw)
        neigh_w[n] = calculate_weigted_neighbor_list(vis_freq_map)

for p in patients:
    for m in mutate_genes[p]:
        for n in neighbors[m]:
            w[n] += neigh_w[m][n]
"""
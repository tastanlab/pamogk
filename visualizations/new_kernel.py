from pamogk.kernels.new_kernel import *
import matplotlib.pyplot as plt

conf_def = 0.1
patient_map = read_data()

# Patient ve mutated genleri yaziyor
patients = preprocess_patient_data(patient_map)

# Pathwayler geldi graphlar ile
all_pw_map = read_pathways()

# list of neigbor of genes for all pathways
neighbor_mappings, id_mapper = get_neighbors_for_all_pathways(all_pw_map, conf_def)

#
kernels = calc_kernel_from_pathways(neighbor_mappings, patients, id_mapper)
for i in range(1):
    print(isPSD(kernels[i]))
    plt.imshow(kernels[i], cmap='hot')
    plt.savefig(f'{i}.png')

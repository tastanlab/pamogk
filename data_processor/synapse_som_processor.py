'''
    This file is to process som file which has type of .maf
    It is enough to call processOneCancerSomatic with file location as argument to function:
        import synapse_som_processor as ssp
        ssp.processOneCancerSomatic("folder/fileName.maf")
        This returns 2d numpy array with columns GeneName, EntrezId, Patient Id

    To write the report to a file:
        import synapse_som_processor as ssp
        output = ssp.processOneCancerSomatic("folder/fileName.maf")
        ssp.writeToFile(output, fileLocation)
'''

import csv
import numpy as np
import argparse
import matplotlib.pyplot as plt
import os

#CancerTypes = ["BLCA","BRCA","COAD","GBM","HNSC","KIRC","LAML","LUAD","LUSC","OV","READ","UCEC"]
CancerTypes = ["BLCA","COAD","GBM","HNSC","LAML","LUAD","LUSC","OV","READ","UCEC"]
parser = argparse.ArgumentParser(description='Run SPK algorithms on pathways')
parser.add_argument('--somatic-data', '-r', metavar='file-path', dest='somatic_data', type=str, help='Somatic Data', default='../data/kirc_data/kirc_somatic_mutation_data.csv')

args = parser.parse_args()
def processOneCancerSomatic(fileLoc):
    filename =fileLoc
    delimiter = "\t"
    startRow = 2
    dataArray = []

    with open(filename, "r",encoding="utf8", errors='ignore') as file:
        csv_reader = csv.reader(file, delimiter=delimiter)
        for i in range(startRow-1):
            next(csv_reader)
        for row in csv_reader:
            #print(row[0])
            if(len(row)>1):
                tempList = [row[0]]+[row[1]]+[row[15]]
                dataArray.append(tempList)

    output = np.array(dataArray)
    prune_ind = np.zeros(len(output))
    prune_list = []
    for idx, data in enumerate(output):
        tmp = data[2]
        splitted = tmp.split("-")
        if '01' in splitted[3]:
            prune_ind[idx] = 1
            newVal = "-"
            newVal = newVal.join(splitted[0:3])
            data[2] = newVal
        else:
            prune_list.append(idx)

    output = np.delete(output, prune_list, axis=0)

    return np.sort(output,axis=0)

#'kirc/data/synapse_kirc_som_data.csv'
def writeToFile(report,fileLoc):
    with open(fileLoc,"w", newline='') as csvfile:
        spamwriter = csv.writer(csvfile, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        spamwriter.writerow(["Gene Name","Entrez Gene ID","Patient ID"])
        for row in report:
            spamwriter.writerow(row)

def printReport(report):
    print(len(report))
    print(len(list(set(report[:,0]))))
    print(len(list(set(report[:,2]))))

def reportAllCancerTypes(folderLoc):
    for cancer in CancerTypes:
        fileLoc = folderLoc+cancer+"/som.maf"
        report = processOneCancerSomatic(fileLoc)
        print(cancer + "# of row, # of Unique Gene, #of Unique Patient")
        printReport(report)


def read_processed_data():
    ### Real Data ###
    # process RNA-seq expression data
    patients = {}
    with open(args.somatic_data) as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            pat_id = row['Patient ID']
            ent_id = row['Entrez Gene ID']
            if pat_id not in patients:
                patients[pat_id] = set([ent_id])
            else:
                patients[pat_id].add(ent_id)

    return patients

def drawHist(pat_dict):
    count_mutated_genes = []
    for patient in pat_dict.keys():
        count_mutated_genes.append(len(pat_dict[patient]))
    count_mutated_genes = sorted(count_mutated_genes)
    np.savetxt("mutated_genes.txt",np.array(count_mutated_genes))
    plt.hist(count_mutated_genes,bins=max(count_mutated_genes),edgecolor='black', linewidth=0.2)  # arguments are passed to np.histogram
    plt.xlabel("Mutasyona Uğramış Gen Sayısı")
    plt.ylabel("Hasta Sayısı")
    plt.savefig('a.png')
    plt.show()



def ensure_dir(file_path):
    directory = os.path.dirname(file_path)
    if not os.path.exists(directory):
        os.makedirs(directory)

def writeOne(cancer):
    folderLoc = "/home/yitepeli/ForExp/"
    #cancer = "BRCA"
    outFile = "../data/"+cancer.lower()+"_data/"+cancer.lower()+"_somatic_mutation_data.csv"
    ensure_dir(outFile)
    fileLoc = folderLoc + cancer + "/som.maf"
    rep = processOneCancerSomatic(fileLoc)
    writeToFile(rep, outFile)

def writeAll():
    for cancer in CancerTypes:
        writeOne(cancer)
#patientDict = read_processed_data()
#drawHist(patientDict)
#reportAllCancerTypes(folderLoc)
writeAll()
print()
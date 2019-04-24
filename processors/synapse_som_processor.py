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

CancerTypes = ["BLCA","BRCA","COAD","GBM","HNSC","KIRC","LAML","LUAD","LUSC","OV","READ","UCEC"]
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

    return output

#'kirc/data/synapse_kirc_som_data.csv'
def writeToFile(report,fileLoc):
    with open(fileLoc,"w", newline='') as csvfile:
        spamwriter = csv.writer(csvfile, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
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

folderLoc = "/home/yitepeli/ForExp/"
#reportAllCancerTypes(folderLoc)
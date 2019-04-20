import csv
import numpy as np

'''
    This file is to process clincal file
    It is enough to call processOneClinicalSomatic with file location as argument to function:
        import synapse_som_processor as ssp
        ssp.processOneCancerSomatic("folder/fileName.maf")
        This returns 2d numpy array with columns PatientId, AliveOrDead(0 or 1), days_to_last_followup or days_to_death(alive or dead), Stage

    To write the report to a file:
        import processors.synapse_clinical_processor as clinic
        output = clinical.processOneClinicalSomatic("folder/fileName.maf")
        ssp.writeToFile(output, fileLocation)



'''

CancerTypes = ["BLCA","BRCA","COAD","GBM","HNSC","KIRC","LAML","LUAD","LUSC","OV","READ","UCEC"]
def processOneClinicalSomatic(fileLoc):
    filename =fileLoc
    delimiter = "\t"
    startRow = 2
    dataArray = []

    with open(filename, "r", encoding="utf8", errors='ignore') as file:
        csv_reader = csv.reader(file, delimiter=delimiter)
        for i in range(startRow - 1):
            header = next(csv_reader)
        for row in csv_reader:
            dataArray.append(row)

    if "tumor_stage" in header:
        TSStr = "tumor_stage"
    elif "ajcc_neoplasm_disease_stage" in header:
        TSStr = "ajcc_neoplasm_disease_stage"
    elif "gynecologic_tumor_grouping_figo_stage" in header:
        TSStr = "gynecologic_tumor_grouping_figo_stage"
    else:
        TSStr = "#"
    dataIndexes = [header.index("#"), header.index("vital_status"), header.index("days_to_death"), header.index("days_to_last_followup"), header.index(TSStr)]

    '''
    for idx, col in enumerate(header):
        if col == "#":
            dataIndexes[0] = idx
        if col == "tumor_stage":
            dataIndexes[4] = idx
        if col == "vital_status":
            dataIndexes[1] = idx
        if col == "days_to_death":
            dataIndexes[2] = idx
        if col == "days_to_last_followup":
            dataIndexes[3] = idx

    '''

    dataArrayProcessed = []
    for row in dataArray:
        tmp = row[dataIndexes[0]]
        splitted = tmp.split("-")
        if (len(splitted)>3) and ('01' in splitted[3]):
            newVal = "-"
            newVal = newVal.join(splitted[0:3])
            if row[dataIndexes[1]] == "LIVING":
                temp_days = row[dataIndexes[3]]
                row3Value = 0
            else:
                temp_days = row[dataIndexes[2]]
                row3Value = 1
            stage = "NA"
            if row[dataIndexes[4]] != "NA":
                checker = row[dataIndexes[4]]
                if ("iv" in checker) or ("IV" in checker):
                    stage = 4
                elif ("iii" in checker) or ("III" in checker):
                    stage = 3
                elif ("ii" in checker) or ("II" in checker):
                    stage = 2
                elif ("i" in checker) or ("I" in checker):
                    stage = 1
                '''
                listFrom = ["Stage I", "Stage II", "Stage III", "Stage IV"]
                if row[dataIndexes[4]] in listFrom:
                    indforTo = listFrom.index(row[dataIndexes[4]])
                    listTo = [1, 2, 3, 4]
                    stage = listTo[indforTo]
                '''
            tempList = [newVal] + [row3Value] + [temp_days] + [stage]
            dataArrayProcessed.append(tempList)


    return np.array(dataArrayProcessed)

#'kirc/data/synapse_kirc_som_data.csv'
def writeToFile(report,fileLoc):
    with open(fileLoc,"w", newline='') as csvfile:
        spamwriter = csv.writer(csvfile, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        for row in report:
            spamwriter.writerow(row)

def printReport(report):
    if len(report) != 0:
        print(len(report))
        print(len(list(set(report[:,0]))))
        aliveCount = np.count_nonzero(report[:, 1] == '0')
        deadCount = np.count_nonzero(report[:, 1] == '1')
        stg1Count = np.count_nonzero(report[:, 3] == '1')
        stg2Count = np.count_nonzero(report[:, 3] == '2')
        stg3Count = np.count_nonzero(report[:, 3] == '3')
        stg4Count = np.count_nonzero(report[:, 3] == '4')
        print(aliveCount)
        print(deadCount)
        print(stg1Count)
        print(stg2Count)
        print(stg3Count)
        print(stg4Count)



def reportAllCancerTypes(folderLoc):
    for cancer in CancerTypes:
        fileLoc = folderLoc+cancer+"/clinical"
        report = processOneClinicalSomatic(fileLoc)
        print(cancer + "# of row, # of Unique Gene, #of Unique Patient")
        printReport(report)
        deneme = "ax"

folderLoc = "/home/yitepeli/ForExp/"
reportAllCancerTypes(folderLoc)

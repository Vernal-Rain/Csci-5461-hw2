import csv

rnaseq = 'HiSeqV2.txt'
microarray = 'HT_HG-U133A.txt'
clinical = 'ov_tcga_clinical_data.tsv'


def make_table(file):
    f = open(file, 'r')
    table = []
    for line in f:
        table.append(line.split())
    f.close()
    return table


def filter(data, func):
    filtered = []#sample numbers only
    for i in range(1, len(data)):
        if func(i) == True:
            filtered.append(data[i][2])
    return filtered


def make_line(data, i):
    string = ''
    for j in range(1, len(data)):
        string += (data[j][i] + ' ')
    return string[:-1] + '\n'


def data_processing(rnaseq, microarray, clinical):
    tsv_file = open(clinical)
    read_tsv = csv.reader(tsv_file, delimiter='\t')
    data = []
    for row in read_tsv:
        data.append(row)
    tsv_file.close()

    survival_time = data[0].index('Overall Survival (Months)')
    survival_status = data[0].index('Overall Survival Status')
    #Deceased && survival < 36 mo
    group1 = filter(data, lambda x: data[x][survival_status] == '1:DECEASED' and (float(data[x][survival_time]) < 36 if data[x][survival_time] != 'NA' else False))
    #survival > 36 mo
    group2 = filter(data, lambda x: float(data[x][survival_time]) > 36 if data[x][survival_time] != 'NA' else False)

    seq = make_table(rnaseq)
    f = open('SeqData.txt', 'w')
    f.write('sampleID' + ' groupID ')
    f.write(make_line(seq, 0))

    for i in range(1, len(seq[0])):
        if seq[0][i] in group1:
            f.write(seq[0][i] + ' 1 ')
            f.write(make_line(seq, i))
        if seq[0][i] in group2:
            f.write(seq[0][i] + ' 2 ')
            f.write(make_line(seq, i))

    f.close()

    array = make_table(microarray)
    f = open('ArrayData.txt', 'w')
    f.write('sampleID' + ' groupID ')
    f.write(make_line(array, 0))

    for i in range(1, len(array[0])):
        if array[0][i] in group1:
            f.write(array[0][i] + ' 1 ')
            f.write(make_line(array, i))
        if array[0][i] in group2:
            f.write(array[0][i] + ' 2 ')
            f.write(make_line(array, i))
    return

if __name__ == '__main__':
    data_processing(rnaseq, microarray, clinical)

from DataProcessing import make_table
from scipy import stats
from numpy import array


def non_zero(array):
    for thing in array:
        if thing != 0:
            return True
    return False


def my_bonferroni(file):
    data = make_table(file)
    bf_values = {}
    num_genes = len(data[0])-2
    for i in range(2, len(data[0])):
        gene_name = data[0][i]
        temp1 = []
        temp2 = []
        for j in range(1, len(data)):
            try:
                if data[j][1] == '1':
                    temp1.append(float(data[j][i]))
                elif data[j][1] == '2':
                    temp2.append(float(data[j][i]))
            except IndexError:
                pass
        if non_zero(temp1) or non_zero(temp2):
            temp1 = array(temp1)
            temp2 = array(temp2)
            t_test = stats.ttest_ind(temp1, temp2, nan_policy='omit')
            if t_test[1] < (0.05/num_genes):
                bf_values[gene_name] = t_test[1]
    bf_values = dict(sorted(bf_values.items(), key=lambda item: item[1]))
    print('Selected genes (p < 0.05): ', len(bf_values.keys()))
    for thing in bf_values:
        print(thing, bf_values[thing])
    return bf_values


if __name__ == '__main__':
    print('RNA seq')
    my_bonferroni('SeqData.txt')
    print('\nMicroarray')
    my_bonferroni('ArrayData.txt')

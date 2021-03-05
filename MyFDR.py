from DataProcessing import make_table
from scipy import stats
from numpy import array


def non_zero(array):
    for thing in array:
        if thing != 0:
            return True
    return False


def t_test(file):
    data = make_table(file)
    t_values = {}
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
            if t_test[1] < 0.05:
                t_values[gene_name] = t_test[1]
    sorted_p_values = dict(sorted(t_values.items(), key=lambda item: item[1]))
    return (sorted_p_values, num_genes)


def my_fdr(file, num_select=0):
    p_vals, n = t_test(file)
    if num_select == 0:
        for i in range(n):
            selected = list(p_vals.items())[:(i+1)]
            p_i = float(selected[-1][1])
            fdr = 1.0 * n * p_i / (i+1)
            if fdr > 0.05:
                print('Genes selected: ', i)
                print('Upper bound of fdr: ', fdr)
                return fdr
    selected = list(p_vals.items())[:num_select]
    p_i = float(selected[-1][1])
    fdr = 1.0 * n * p_i / num_select
    print('Genes selected: ', num_select)
    print('Upper bound of fdr: ', fdr)
    return fdr


if __name__ == '__main__':
    nums = [20, 50, 100, 200]
    print('RNA seq')
    for i in nums:
        my_fdr('SeqData.txt', i)
    print('\nMicroarray')
    for i in nums:
        my_fdr('ArrayData.txt', i)

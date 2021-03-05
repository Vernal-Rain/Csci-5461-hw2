from scipy import stats
from numpy import array
import matplotlib.pyplot as plt
from DataProcessing import make_table


def non_zero(array):
    for thing in array:
        if thing != 0:
            return True
    return False


def stat_tests(data):
    t_test_values = {}
    wilcoxon_values = {}
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
            wilcoxon = stats.ranksums(temp1, temp2)
            t_test_values[gene_name] = t_test[1]
            wilcoxon_values[gene_name] = wilcoxon[1]
    return (t_test_values, wilcoxon_values)


def top_10(sorted_dict):
    top_10 = list(sorted_dict.items())[:10]
    for i in top_10:
        print(i)
    return


def filter(dict):
    filtered = {}
    for i in dict.keys():
        if dict[i] < 0.05:
            filtered[i] = dict[i]
    return filtered


def sort(tuple):
    t_test, wilcoxon = tuple
    t_test = filter(t_test)
    wilcoxon = filter(wilcoxon)
    t_test = dict(sorted(t_test.items(), key=lambda item: item[1]))
    wilcoxon = dict(sorted(wilcoxon.items(), key=lambda item: item[1]))
    print('Total genes selected (t-test, p < 0.05):', len(t_test))
    print('Top 10 (t-test)')
    top_10(t_test)
    print('\nTotal genes selected (wilcoxon, p < 0.05):', len(wilcoxon))
    print('Top 10 (wilcoxon)')
    top_10(wilcoxon)
    return (t_test, wilcoxon)


def plot(dict):
    plt.hist(dict.values(), bins=50)
    plt.show()


def report_data(file):
    data = stat_tests(make_table(file))
    print(file)
    t_test, wilcoxon = data
    plot(t_test)
    plot(wilcoxon)
    data = sort(data)

def get_data(file):
    data = stat_tests(make_table(file))
    data = sort(data)
    return data


if __name__ == '__main__':
    report_data('SeqData.txt')
    print('\n')
    report_data('ArrayData.txt')

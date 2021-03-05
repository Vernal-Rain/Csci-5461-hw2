from StatTests import *


def select_n(dict, n: int):
    ls = list(dict.items())[:n]
    selected = []
    for thing in ls:
        selected.append(thing[0])
    return selected


if __name__ == '__main__':
    x = [i*10 for i in range(1, 101)]
    y = []

    t, w = stat_tests(make_table('SeqData.txt'))
    t = filter(t)
    seq = dict(sorted(t.items(), key=lambda item: item[1]))
    t, w = stat_tests(make_table('ArrayData.txt'))
    t = filter(t)
    array = dict(sorted(t.items(), key=lambda item: item[1]))

    for i in x:
        count = 0
        seq_n = select_n(seq, i)
        array_n = select_n(array, i)
        for thing in seq_n:
            if thing in array_n:
                count += 1
        y.append(count)
    plt.scatter(x, y)
    plt.show()

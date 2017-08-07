import sys

# sums all columns in biom table into one column
def collapse_biom(biom,out):
    with open(biom) as b:
        temp = b.readline()
        temp = temp.split(']}}],"columns":')[0]
        temp = temp.split(']],"rows": [{"id": "')
        rows = temp[1]
        data = temp[0].split(',"data": [[')[1]
    # Data is in the next format:
    # x1, x2, x3
    # x1, x2 are matrix coordinates
    # x3 is number of reads which is always 1
    # as in biom each column is a single read and each row is a single taxon
    data = data.split('],[')
    rows = rows.split(']}},{"id": "')
    # Go through matrix and add up all values for each taxon
    values = {}
    prev_line = 0
    last_column = 0
    for el in data:
        temp = el.split(',')
        now_line = int(temp[0])
        if now_line != prev_line:
            now_column = int(temp[1])
            # Assign as many count to a taxon as many reads (columns) were asssigned to it
            # (Every column may have no more than 1 non zero value which is 1)
            num = now_column - last_column
            values[rows[prev_line]] = num
            last_column = now_column
        prev_line = now_line
    values[rows[prev_line]] = num
    # print(values)
    lines = []
    with open(out, 'a') as t:
        for el in values:
            temp = el.split('", "metadata": {"taxonomy": [')
            lines.append(temp[0] + '\t' + str(values[el]) + '\t' + temp[1] + '\n')
        # print(lines)
        lines.sort(key = lambda x: int(x.split('\t')[1]), reverse=True)
        for line in lines:
            t.write(line)


collapse_biom(sys.argv[1],sys.argv[2])
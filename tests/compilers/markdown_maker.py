
'''
adjusts COMPILED/RAN cells in compability.csv to be prettier checkboxes
when data is ultimately displayed as a markdown table (e.g. through
https://donatstudios.com/CsvToMarkdownTable), along with other
cosmetic changes.
'''

# this script is a quick hack - please don't judge it :)

fn_in = 'compatibility.csv'
fn_out = 'compatibility_reformatted.csv'

str_0 = ''
str_1 = '<ul><li>[x] </li>'

with open(fn_in, 'r') as f:
    data = f.read() + '\n'

# adjust COMPILED/RAN columns
data = data.replace(',0,', ',%s,'  % str_0)
data = data.replace(',1,', ',%s,'  % str_1)
data = data.replace(',0\n',',%s\n' % str_0)
data = data.replace(',1\n',',%s\n' % str_1)

# adjust LANGUAGE column
data = data.replace('\nc,', '\nC,')
data = data.replace('\ncpp,', '\nC++,')

# shave gcc from compiler versions
data = data.replace(',gcc/', ',')

# shave gpu/cuda from nvcc versions
data = data.replace(',gpu/cuda/', ',')

# shave intel-compilers from icc versions
data = data.replace(',intel-compilers/', ',')

# change case of wrapper compilers
data = data.replace(',MPICC,', ',mpicc,')
data = data.replace(',NVCC,', ',nvcc,')

# shave of wrapper module's corresponding c/c++ compiler label
data = '\n'.join([
    ','.join([
        item.split('__')[0] for item in line.split(',')
    ]) for line in data.split('\n')])

with open(fn_out, 'w') as f:
    f.write(data)


'''
adjusts COMPILED/RAN cells in compability.csv to be prettier checkboxes
when data is ultimately displayed as a markdown table (e.g. through
https://donatstudios.com/CsvToMarkdownTable)
'''

fn_in = 'compatibility.csv'
fn_out = 'compatibility_md.csv'

str_0 = ''
str_1 = '<ul><li>[x] </li>'

with open(fn_in, 'r') as f:
    data = f.read() + '\n'

data = data.replace(',0,', ',%s,'  % str_0)
data = data.replace(',1,', ',%s,'  % str_1)
data = data.replace(',0\n',',%s\n' % str_0)
data = data.replace(',1\n',',%s\n' % str_1)

with open(fn_out, 'w') as f:
    f.write(data)

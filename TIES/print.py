import pprint
import ast

with open('result.dat', 'r') as file:
    data = file.read()

data_dict = ast.literal_eval(data)

pprint.pprint(data_dict)

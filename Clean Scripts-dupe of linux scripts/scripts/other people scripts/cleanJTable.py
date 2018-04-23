# convert the table results to a dataframe
# have option to save it as a csv or similar file type

import pandas as pd

jTab = plm.table()
# Note that the headers will be in random sequence due to the nature of dictionaries
headers = list(jTab[0].keys()) #you can also hardcode in the headers, or take in parameters from the user input
# ie. take values from custom=[('', []), ]
dataList = []
for row in jTab:
    newRow = []
    for key in headers:
        newRow += [row[key]]
    dataList += [newRow]

jdf =pd.DataFrame(dataList, columns=headers)
# you can remove the column that has the key '' here. Or never pass it to start (see the header assignment)

# write to a file. I think to_csv will work even if the user doesn't put in the csv extension, though you can also force the extension type
if save==True:
    # assume filename is provided. Too lazy to write this part out
    jdf.to_csv(filename, sep = ',', index = False)

## Examples of midway results:
# >>> headers
# ['latitude', '', 'filtered_height', 'Model', 'longitude']

# >>> dataList
# [[49.135, [], nan, 3800.4956, -114.113], [49.185, [], 2822.0, 2849.8267, -114.087], [49.175, [], 3023.0, 2654.6406, -114.09], [49.165, [], nan, 3610.3496, -114.093], ...

# >>> jdf
#    latitude      filtered_height        Model  longitude
#0     49.135  []              NaN  3800.495605   -114.113
#1     49.185  []             2822  2849.826660   -114.087
#2     49.175  []             3023  2654.640625   -114.090
#3     49.165  []              NaN  3610.349609   -114.093
#4     49.156  []              NaN  3581.041260   -114.096

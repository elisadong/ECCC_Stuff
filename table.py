import pandas as pd
import json

def makeJson(values, headers):
    json = []
    for row in values:
        jsonDict = {}
        for header in headers:
            jsonDict[header] = row[headers.index(header)]
        json += jsonDict
    return json


def table(self, option='filtered_height', model=True, output='json', full=False,custom=[], save='',fileName=None):
# custom should be a liust of dictionaries, where the the dict val for each key is as long as the total number of lon/lat values

# writing preliminary headers this way so it's easy to add mroe options if provided in list format
headers = ['Longitude', 'Latitude'] + [option]
# update the headers with the custom informatioln, the following should get the only key for each dictionary and add it do the headers list
# this should get skipped if custom is an empty list
for field in custom:
    headers += list[field.keys()]

# assuming only one option, otherwise rewrite to deal with the list case v
option = arg_to_list(option)

lons = self.longitude
lats = self.latitude

dataList = []
if model:
    for lon in lons:
        ind = lons.index(lon)
        lonVal = lons[ind]
        latVal = lats[ind]
        optVal = option[ind]
        rowVals = [lonVal, latVal, optVal]
        if custom != []:
            cVals = []
            for c in custom:
                cKey = c.keys()
                cVals += c[cKey]
            rowVals += cVals
        dataList += [rowVals]

df = pd.DataFrame(dataList, columns =headers)
if save ==True:
    if output == 'csv':
        if not filename:
            filename =raw_input('Enter CSV filename: ')
        df.to_csv(filename,sep=',',index=False)
    elif output == 'json':
        js = makeJson(headers, dataList)
        if not filename:
            filename = raw_input('Enter desired JSON filename: ')

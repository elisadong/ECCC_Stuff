# plot scatter plots for time stamps 0-23
# plot by every fuel type
# csvDates should exist as folders within csvFolder and are expected to be found in the csvFile path

import os
import pandas as pd
import matplotlib.pyplot as plt

# csvDate = '20170912'
csvDates = ['20170901', '20170902', '20170903', '20170904', '20170905', '20170906','20170907','20170908', '20170909', '20170910', '20170911','20170912','20170913','20170914','20170915']
csvFolder = '/fs/site2/dev/eccc/aq/r2/jac001/FireWork/CFFEPS/From_Kerry/20171114/outputs.persistence/FST/'

for date in csvDates:
    print('Now processing date: {}...'.format(date))
    csvFile = os.path.join(csvFolder, date, 'CFFEPS_{}.csv'.format(date))

    cdf = pd.read_csv(csvFile)
    cHeaders = list(cdf.columns.values)

    fuelInd = 10
    zPlumeInd = 36
    UTCInd = 17
    cdf['times'] = cdf[cHeaders[UTCInd]].map(lambda x: int(x[-5:-3]))

    for fuel, group in cdf.groupby(cdf[cHeaders[fuelInd]]):
        plt.figure()
        fig=group.plot(kind='scatter',x='times',y=cHeaders[zPlumeInd],title=fuel).get_figure()
        plt.xlabel('Time Step')
        axes=plt.gca()
        axes.set_xlim(left=0,right=23)
        saveName = '{}_{}_scatterPlot.png'.format(date,fuel.strip())
        fig.savefig(saveName)
        plt.close('all')

print('Done!')

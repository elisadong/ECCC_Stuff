import os
import rpnpy.librmn.all as rmn

## Test Constants:
defaultDir=os.getcwd()
defaultIp1=76696048
defaultSpc='TO3'
defaultFile='2010071000_001'


## getConc (level,spc,fileName) returns an array of concentration data with dimensions of the file shape attirbute
def getConc(level=defaultIp1, spc=defaultSpc, fileName=defaultFile):
    # may need full directory path to the file, thus may need to move the path where this program is taking place
    fileID = rmn.fstopenall(fileName, rmn.FST_RO) # opens file within the directory, else may need path
    dataKey = rmn.fstinf(fileID, nomvar=spc, ip1=level)['key'] # get key for the matching data
    # Note may need to turn this into a loop if it matches more than one
    dataRec = rmn.fstluk(dataKey)
    concData = dataRec['d']
    return concData, dataKey


>>> keyList = rmn.fstinl(fileOpen,nomvar=defaultSpc, ip1
... =defaultIp1)
>>> keyList
[6311937]

>>> rec = rmn.fstluk(int(keyList[0]))
Read(999) TO3  P  GM_DEV         772   642     1  363996800    76696048       1     0      300       12  f 16  Z 60193 56056     1     0


>>> O3=rmn.fstinf(fileOpen,nomvar=defaultSpc,ip1=defaultIp1)
>>> O3
{'shape': (772, 642, 1), 'key': 6311937}

>>> O3Key=rmn.fstinf(fileOpen,nomvar=defaultSpc,ip1=defaultIp1)
>>> O3Key
{'shape': (772, 642, 1), 'key': 6311937}

>>> O3Key=O3Key['key']
>>> O3Key
6311937

>>> O3Key
6311937
>>> O3Meta=rmn.fstprm(O3Key)


>>> O3Meta
{'deet': 300, 'shape': (772, 642, 1), 'typvar': 'P ', 'lng': 143566, 'ni': 772, 'nj': 642, 'nk': 1, 'swa': 110332609, 'datyp': 134, 'xtra1': 363997700, 'xtra2': 0, 'xtra3': 0, 'ip2': 1, 'ip3': 0, 'nbits': 16, 'key': 6311937, 'ubc': 0, 'npas': 12, 'ig4': 0, 'ig3': 1, 'ig2': 56056, 'ig1': 60193, 'datev': 363997700, 'nomvar': 'TO3 ', 'ip1': 76696048, 'dltf': 0, 'etiket': 'GM_DEV      ', 'grtyp': 'Z', 'dateo': 363996800}


>>> O3Rec=rmn.fstluk(O3Key)
Read(999) TO3  P  GM_DEV         772   642     1  363996800    76696048       1     0      300       12  f 16  Z 60193 56056     1     0
>>> O3Rec['d']
array([[ 27.0703125 ,  27.0703125 ,  27.0703125 , ...,  31.35546875,
         31.35546875,  31.35546875],
       [ 27.0703125 ,  27.0703125 ,  27.0703125 , ...,  31.35546875,
         31.35546875,  31.35546875],
       [ 27.0703125 ,  27.0703125 ,  27.01171875, ...,  31.36328125,
         31.42578125,  31.42578125],
       ...,
       [ 38.28515625,  38.28515625,  38.08984375, ...,  56.9296875 ,
         57.234375  ,  57.3125    ],
       [ 38.28515625,  38.28515625,  38.28515625, ...,  57.2578125 ,
         57.23046875,  57.234375  ],
       [ 38.28515625,  38.28515625,  38.28515625, ...,  57.23828125,
         57.23046875,  57.23046875]], dtype=float32)
>>> ppO3=pprint.pformat(O3Rec)
>>> pp03
Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
NameError: name 'pp03' is not defined
>>> ppO3
"{'d': array([[ 27.0703125 ,  27.0703125 ,  27.0703125 , ...,  31.35546875,\n         31.35546875,  31.35546875],\n       [ 27.0703125 ,  27.0703125 ,  27.0703125 , ...,  31.35546875,\n         31.35546875,  31.35546875],\n       [ 27.0703125 ,  27.0703125 ,  27.01171875, ...,  31.36328125,\n         31.42578125,  31.42578125],\n       ..., \n       [ 38.28515625,  38.28515625,  38.08984375, ...,  56.9296875 ,\n         57.234375  ,  57.3125    ],\n       [ 38.28515625,  38.28515625,  38.28515625, ...,  57.2578125 ,\n         57.23046875,  57.234375  ],\n       [ 38.28515625,  38.28515625,  38.28515625, ...,  57.23828125,\n         57.23046875,  57.23046875]], dtype=float32),\n 'dateo': 363996800,\n 'datev': 363997700,\n 'datyp': 134,\n 'deet': 300,\n 'dltf': 0,\n 'etiket': 'GM_DEV      ',\n 'grtyp': 'Z',\n 'ig1': 60193,\n 'ig2': 56056,\n 'ig3': 1,\n 'ig4': 0,\n 'ip1': 76696048,\n 'ip2': 1,\n 'ip3': 0,\n 'key': 6311937,\n 'lng': 143566,\n 'nbits': 16,\n 'ni': 772,\n 'nj': 642,\n 'nk': 1,\n 'nomvar': 'TO3 ',\n 'npas': 12,\n 'shape': [772, 642],\n 'swa': 110332609,\n 'typvar': 'P ',\n 'ubc': 0,\n 'xtra1': 363997700,\n 'xtra2': 0,\n 'xtra3': 0}"

O3Rec['d'][0]
array([ 27.0703125 ,  27.0703125 ,  27.0703125 ,  27.0703125 ,
        27.0703125 ,  27.0703125 ,  27.0703125 ,  27.0703125 ,
        27.0703125 ,  27.0703125 ,  27.0703125 ,  27.0703125 ,
        27.0703125 ,  27.0703125 ,  27.0703125 ,  27.0703125 ,
        27.07421875,  26.94140625,  26.8984375 ,  28.59765625,
        30.06640625,  29.9140625 ,  29.8125    ,  29.82421875,
        29.82421875,  29.82421875,  29.82421875,  29.82421875,
        29.82421875,  29.8203125 ,  29.8203125 ,  29.82421875,
        29.83203125,  29.83984375,  29.78125   ,  29.671875  ,
        29.63671875,  29.6328125 ,  29.63671875,  29.640625  ,
        29.640625  ,  29.640625  ,  29.64453125,  29.5078125 ,
        29.4296875 ,  31.28515625,  33.0234375 ,  32.875     ,
        32.74609375,  32.7578125 ,  32.7578125 ,  32.7578125 ,
        32.7578125 ,  32.7578125 ,  32.7578125 ,  32.7578125 ,
        32.7578125 ,  32.7578125 ,  32.7578125 ,  32.7578125 ,
        32.7578125 ,  32.7578125 ,  32.7578125 ,  32.7578125 ,
        32.7578125 ,  32.7578125 ,  32.7578125 ,  32.7578125 ,
        32.76171875,  32.66796875,  32.59375   ,  33.91015625,
        35.23828125,  35.1484375 ,  35.04296875,  35.05078125,
        35.046875  ,  35.109375  ,  35.05078125,  34.57421875,
        34.26171875,  34.1953125 ,  34.20703125,  34.22265625,
        34.22265625,  34.22265625,  34.22265625,  34.22265625,
        34.22265625,  34.22265625,  34.22265625,  34.22265625,
        34.22265625,  34.22265625,  34.22265625,  34.19921875,
        34.171875  ,  34.5390625 ,  34.94140625,  34.921875  ,
        34.88671875,  34.890625  ,  34.890625  ,  34.890625  ,
        34.890625  ,  34.890625  ,  34.890625  ,  34.890625  ,
        34.890625  ,  34.890625  ,  34.890625  ,  34.890625  ,
        34.890625  ,  34.890625  ,  34.890625  ,  34.890625  ,
        34.890625  ,  34.890625  ,  34.890625  ,  34.88671875,
        34.98828125,  34.8125    ,  34.046875  ,  33.66015625,
        33.7265625 ,  34.078125  ,  34.23828125,  34.203125  ,
        34.203125  ,  34.203125  ,  34.203125  ,  34.203125  ,
        34.203125  ,  34.203125  ,  34.203125  ,  34.203125  ,
        34.203125  ,  34.203125  ,  34.203125  ,  34.203125  ,
        34.203125  ,  34.203125  ,  34.203125  ,  34.203125  ,
        34.203125  ,  34.203125  ,  34.203125  ,  34.19921875,
        34.1953125 ,  34.22265625,  34.40234375,  34.8984375 ,
        35.1015625 ,  35.03125   ,  35.03125   ,  35.03515625,
        35.03515625,  35.03515625,  35.03515625,  35.03515625,
        35.03125   ,  35.08984375,  35.109375  ,  34.51953125,
        33.9609375 ,  33.82421875,  33.82421875,  33.859375  ,
        33.859375  ,  33.859375  ,  33.859375  ,  33.859375  ,
        33.859375  ,  33.85546875,  33.85546875,  33.8671875 ,
        33.8671875 ,  33.8359375 ,  33.68359375,  33.49609375,
        33.4921875 ,  33.51171875,  33.51171875,  33.51171875,
        33.51171875,  33.51171875,  33.51171875,  33.51171875,
        33.51171875,  33.51171875,  33.51171875,  33.51171875,
        33.51171875,  33.51171875,  33.51171875,  33.51171875,
        33.51171875,  33.51171875,  33.5078125 ,  33.53125   ,
        33.6328125 ,  33.12109375,  32.3359375 ,  32.3515625 ,
        32.25390625,  31.80078125,  31.01171875,  30.58203125,
        30.640625  ,  30.65234375,  30.65234375,  30.65234375,
        30.65234375,  30.65234375,  30.65234375,  30.65234375,
        30.65234375,  30.65234375,  30.65234375,  30.65234375,
        30.65234375,  30.65234375,  30.65234375,  30.65234375,
        30.65234375,  30.65234375,  30.65234375,  30.65234375,
        30.6484375 ,  30.6484375 ,  30.65234375,  30.65625   ,
        30.640625  ,  30.58203125,  30.48828125,  30.46875   ,
        30.4765625 ,  30.4765625 ,  30.56640625,  30.3671875 ,
        29.7109375 ,  29.39453125,  29.3359375 ,  29.36328125,
        29.37890625,  29.375     ,  29.375     ,  29.375     ,
        29.375     ,  29.375     ,  29.375     ,  29.375     ,
        29.375     ,  29.375     ,  29.375     ,  29.375     ,
        29.375     ,  29.375     ,  29.38671875,  29.3828125 ,
        29.3203125 ,  29.1015625 ,  28.9453125 ,  28.96484375,
        28.9765625 ,  28.9765625 ,  28.9765625 ,  28.9765625 ,
        28.9765625 ,  28.9765625 ,  28.9765625 ,  28.9765625 ,
        28.96875   ,  29.125     ,  29.00390625,  27.70703125,
        26.828125  ,  26.6484375 ,  26.671875  ,  26.72265625,
        26.71875   ,  26.71875   ,  26.71875   ,  26.71875   ,
        26.71875   ,  26.71484375,  26.71484375,  26.890625  ,
        26.3984375 ,  24.90234375,  24.5078125 ,  24.6953125 ,
        24.6953125 ,  24.69140625,  24.69140625,  24.69140625,
        24.69140625,  24.69140625,  24.69140625,  24.69140625,
        24.69140625,  24.69140625,  24.69140625,  24.69140625,
        24.69140625,  24.6875    ,  24.69140625,  24.7890625 ,
        24.51953125,  23.8125    ,  23.515625  ,  23.46875   ,
        23.5       ,  23.51171875,  23.51171875,  23.51953125,
        23.5078125 ,  23.41796875,  23.22265625,  23.15625   ,
        23.18359375,  23.18359375,  23.18359375,  23.18359375,
        23.18359375,  23.18359375,  23.18359375,  23.18359375,
        23.18359375,  23.18359375,  23.18359375,  23.18359375,
        23.18359375,  23.18359375,  23.18359375,  23.18359375,
        23.18359375,  23.18359375,  23.18359375,  23.171875  ,
        23.1875    ,  23.265625  ,  23.30859375,  23.31640625,
        23.3125    ,  23.29296875,  23.2890625 ,  23.36328125,
        23.6953125 ,  24.109375  ,  24.11328125,  24.06640625,
        24.0703125 ,  24.0703125 ,  24.0703125 ,  24.0703125 ,
        24.0703125 ,  24.0703125 ,  24.0703125 ,  24.0703125 ,
        24.0703125 ,  24.0703125 ,  24.0703125 ,  24.0703125 ,
        24.0703125 ,  24.0703125 ,  24.0703125 ,  24.0703125 ,
        24.0703125 ,  24.0703125 ,  24.078125  ,  24.05078125,
        24.0078125 ,  23.99609375,  23.98046875,  23.95703125,
        23.94140625,  23.984375  ,  24.1875    ,  24.88671875,
        25.31640625,  25.234375  ,  25.2109375 ,  25.21484375,
        25.21484375,  25.21484375,  25.21484375,  25.21484375,
        25.21484375,  25.21484375,  25.21484375,  25.21484375,
        25.21484375,  25.21484375,  25.21484375,  25.21484375,
        25.21484375,  25.21484375,  25.21875   ,  25.2734375 ,
        25.04296875,  24.60546875,  24.5703125 ,  24.61328125,
        24.61328125,  24.61328125,  24.5859375 ,  24.55078125,
        24.65625   ,  25.22265625,  26.25      ,  26.4921875 ,
        26.35546875,  26.35546875,  26.359375  ,  26.359375  ,
        26.359375  ,  26.359375  ,  26.359375  ,  26.359375  ,
        26.359375  ,  26.359375  ,  26.359375  ,  26.359375  ,
        26.359375  ,  26.2890625 ,  26.65625   ,  27.39453125,
        27.66796875,  27.71484375,  27.67578125,  27.66796875,
        27.671875  ,  27.671875  ,  27.671875  ,  27.671875  ,
        27.671875  ,  27.67578125,  27.65234375,  27.59765625,
        27.7890625 ,  28.015625  ,  27.94921875,  28.7421875 ,
        29.83984375,  29.82421875,  29.7109375 ,  29.71875   ,
        29.71875   ,  29.71875   ,  29.71875   ,  29.64453125,
        30.015625  ,  30.80859375,  31.12890625,  31.18359375,
        31.14453125,  31.1328125 ,  31.13671875,  31.13671875,
        31.13671875,  31.13671875,  31.13671875,  31.13671875,
        31.13671875,  31.13671875,  31.13671875,  31.13671875,
        31.13671875,  31.13671875,  31.13671875,  31.140625  ,
        31.10546875,  31.04296875,  30.99609375,  31.10546875,
        31.6015625 ,  33.18359375,  34.65625   ,  35.38671875,
        35.49609375,  35.40625   ,  35.40625   ,  35.41015625,
        35.41015625,  35.41015625,  35.41015625,  35.41015625,
        35.41015625,  35.41015625,  35.41015625,  35.41015625,
        35.41015625,  35.41015625,  35.41015625,  35.41015625,
        35.41015625,  35.41015625,  35.41015625,  35.41015625,
        35.296875  ,  35.6328125 ,  36.42578125,  36.7109375 ,
        36.75390625,  36.74609375,  36.7421875 ,  36.73828125,
        36.73828125,  36.7890625 ,  36.74609375,  36.49609375,
        36.5       ,  36.28125   ,  35.23828125,  34.84375   ,
        34.9921875 ,  34.9921875 ,  34.9921875 ,  34.9921875 ,
        34.9921875 ,  34.9921875 ,  34.9921875 ,  34.9140625 ,
        35.1875    ,  35.75390625,  35.82421875,  35.765625  ,
        35.765625  ,  35.76953125,  35.76953125,  35.76953125,
        35.76953125,  35.76953125,  35.76953125,  35.76953125,
        35.76953125,  35.76953125,  35.76953125,  35.76953125,
        35.76953125,  35.76953125,  35.76953125,  35.76953125,
        35.74609375,  35.75390625,  35.9453125 ,  36.109375  ,
        36.17578125,  36.2109375 ,  36.06640625,  35.7734375 ,
        35.8359375 ,  35.3515625 ,  34.11328125,  33.87109375,
        34.03125   ,  34.03125   ,  34.02734375,  34.02734375,
        34.02734375,  34.02734375,  34.02734375,  34.02734375,
        34.02734375,  34.0234375 ,  34.00390625,  34.11328125,
        34.28515625,  34.29296875,  34.27734375,  34.27734375,
        34.27734375,  34.27734375,  34.27734375,  34.27734375,
        34.27734375,  34.27734375,  34.27734375,  34.27734375,
        34.27734375,  34.27734375,  34.27734375,  34.27734375,
        34.27734375,  34.2734375 ,  34.26953125,  34.31640625,
        34.375     ,  34.375     ,  34.3671875 ,  34.375     ,
        34.4140625 ,  34.45703125,  34.4296875 ,  34.16015625,
        34.078125  ,  34.0390625 ,  32.94921875,  32.19921875,
        32.2734375 ,  32.35546875,  32.3984375 ,  32.39453125,
        32.42578125,  32.3046875 ,  31.9375    ,  31.76171875,
        31.734375  ,  31.734375  ,  31.69921875,  31.734375  ,
        31.73828125,  31.734375  ,  31.70703125,  31.73828125,
        31.73828125,  31.73828125,  31.73828125,  31.73828125,
        31.73828125,  31.73828125,  31.73828125,  31.73828125,
        31.77734375,  31.66015625,  31.37890625,  31.32421875,
        31.35546875,  31.35546875,  31.35546875,  31.35546875,
        31.35546875,  31.35546875,  31.35546875,  31.35546875,
        31.35546875,  31.35546875,  31.35546875,  31.35546875,
        31.35546875,  31.35546875], dtype=float32)
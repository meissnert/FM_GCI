# -*- coding: utf-8 -*-
"""
Tobias Meissner, 2017

querry cancergenomeinterpreter.org using foundation medicine xml files
"""

import sys
from os.path import basename
import re
import requests
from xml.dom import minidom
import tempfile

# get commandline input
xmlFile = sys.argv[1]
email = sys.argv[2]
key = sys.argv[3]
outDir = sys.argv[4]

print(sys.argv[1:])

### read fm fml file and parse variant, alteration and rearrangement info
name = re.sub('.xml', '', basename(xmlFile))

xmldoc = minidom.parse(xmlFile)

varlist = xmldoc.getElementsByTagName('short-variant')
if len(varlist) > 0:
    var = tempfile.NamedTemporaryFile(delete=False)
    var.write('protein\n')
    for s in varlist:
        var.write(s.attributes['transcript'].value + ':p.' +  s.attributes['protein-effect'].value)
        var.write('\n')
    var.close()
    varfile=True
else: varfile=False

if(varfile):
    print('variant input: \n')
    with open(var.name, 'r') as fin:
        print fin.read()
else:
    print('No variant data detexted\n')
    
    
cnalist = xmldoc.getElementsByTagName('copy-number-alteration')
if len(cnalist) > 0:
    cna = tempfile.NamedTemporaryFile(delete=False)
    cna.write('gene\tcna\n')
    for s in cnalist:
        if s.attributes['type'].value == 'amplification':
            cna.write(s.attributes['gene'].value + '\t' + 'amp')
            cna.write('\n')
        if s.attributes['type'].value == 'loss':
            cna.write(s.attributes['gene'].value + '\t' + 'del')
            cna.write('\n')
    cna.close()
    cnafile=True
else: cnafile=False

if(cnafile):    
    print('Copy number input: \n')
    with open(cna.name, 'r') as fin:
        print fin.read()
else:
    print('No copy number data detected\n')
        
    
realist = xmldoc.getElementsByTagName('rearrangement')
if len(realist) > 0:
    rea = tempfile.NamedTemporaryFile(delete=False)
    rea.write('fus\n')
    for s in realist:
        rea.write(s.attributes['targeted-gene'].value + '__' +  s.attributes['other-gene'].value)
        rea.write('\n')
    rea.close()
    reafile=True
else: reafile=False
    
if(reafile):
    print('Rearangements input: \n')
    with open(rea.name, 'r') as fin:
        print fin.read()
else:
    print('No Rearangement data detected\n')
    
### query cancergenomeinerpreter
headers = {'Authorization': email + ' ' + key}
api = 'https://www.cancergenomeinterpreter.org/api/v1/'
#disease = xmldoc.getElementsByTagName('variant-report')[0].attributes['disease'].value
disease='CANCER'

payload = {'cancer_type': disease, 'title': name}

files = {'mutations': open(var.name, 'rb') if varfile else None,
         'cnas': open(cna.name, 'rb') if cnafile else None,
         'translocations': open(rea.name, 'rb') if reafile else None
         }
files = {k:v for k,v in files.items() if v is not None}

r = requests.post(api,
                headers=headers,
                verify=False,
                files=files,
                data=payload)
job = api + r.json()

print('Job submitted \n')
print('JOB ID: ' + r.json() + '\n')

### wait for the job to complete
import time
starttime=time.time()
r = requests.get(job, headers=headers, verify=False)
while r.json()['status']=='Running':
  r = requests.get(job, headers=headers, verify=False)
  print "Waiting for job to complete, cheking again in 15s ..."
  time.sleep(15.0 - ((time.time() - starttime) % 15.0))

### download results
payload={'action':'download'}
r = requests.get(job, headers=headers, params=payload, verify=False)
with open(outDir + '/' + name + '.zip', 'wb') as fd:
    fd.write(r._content)
    
print('done')
print('Output saved to ' + outDir + '/' + name + '.zip')

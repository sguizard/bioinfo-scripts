#! /usr/bin/env python

import sys
import argparse
import re

# Set version and licence
ver='1.0.0'
gpl='gtfEnsembl2bed12.py version ' + ver
gpl=gpl + '''
Copyright (C) 2023  Sébastien Guizard <sguizard@ed.ac.uk>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.'''


# Set options
parser = argparse.ArgumentParser(
    prog='gtfEnsembl2bed12.py', 
    description='Convert Ensembl GTF file into BED12 file.', 
    usage='gtf2bed.py --gtf genes.gtf', 
    epilog=gpl, 
    formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument('-g', '--gtf', help='GTF file to be converted in BED12 format'            , type=argparse.FileType(mode='r'), required=True)
parser.add_argument('-b', '--bed', help='Name of the new BED12 file'                          , type=str                        , required=False)
parser.add_argument('-i', '--inc', help='Include genes? Default = No'                         , type=bool                       , required=False)
parser.add_argument('-r', '--rgb', help='Annotation color. Format \'[0-255],[0-255],[0-255]\'', type=str                        , required=False)

args = parser.parse_args()

# Define functions
def initiateAnnotation(id): 
    bedData[id]                = {}
    bedData[id]['seq']         = ''
    bedData[id]['start']       = int(0)
    bedData[id]['end']         = int(0)
    bedData[id]['name']        = ''
    bedData[id]['score']       = int(0)
    bedData[id]['strand']      = ''
    bedData[id]['thickStart']  = int(0)
    bedData[id]['thickEnd']    = int(0)
    bedData[id]['itemRgb']     = args.rgb
    bedData[id]['blockCount']  = int(0)
    bedData[id]['blockStarts'] = {}

def gene2bed(id): 
    bedStart      = start - 1
    bedBlockSizes = end - bedStart

    bedData[id]['seq']                   = seq
    bedData[id]['start']                 = bedStart
    bedData[id]['end']                   = end
    bedData[id]['name']                  = id
    bedData[id]['score']                 = score
    bedData[id]['strand']                = strand
    bedData[id]['thickStart']            = bedStart
    bedData[id]['thickEnd']              = end
    bedData[id]['blockCount']            = 1
    bedData[id]['blockStarts'][bedStart] = bedBlockSizes

def updateTranscript(id):
    bedStart = start - 1

    bedData[id]['seq']         = seq
    bedData[id]['start']       = bedStart
    bedData[id]['end']         = end
    bedData[id]['name']        = id
    bedData[id]['score']       = score
    bedData[id]['strand']      = strand
    bedData[id]['thickStart']  = 0
    bedData[id]['thickEnd']    = 0
    bedData[id]['blockStarts'] = {}

def addExon(id):
    bsize = end - (start - 1)
    bedStart = start - 1
    bedData[id]['blockCount'] += 1
    bedData[id]['blockStarts'][bedStart] = bsize

def addCDS(id):
    if bedData[id]['thickStart'] == 0:
        bedData[id]['thickStart'] = start -1
    elif bedData[id]['thickStart'] > (start -1):
        bedData[id]['thickStart'] = start -1

    if bedData[id]['thickEnd'] == 0:
        bedData[id]['thickEnd'] = end
    elif bedData[id]['thickEnd'] < end:
        bedData[id]['thickEnd'] = end

def getBed12Line(id):
    bedLine = bedData[key]['seq'] + "\t"
    bedLine = bedLine + str(bedData[id]['start']) + "\t"
    bedLine = bedLine + str(bedData[id]['end']) + "\t"
    bedLine = bedLine + bedData[id]['name'] + "\t"
    bedLine = bedLine + bedData[id]['score'] + "\t"
    bedLine = bedLine + bedData[id]['strand'] + "\t"
    bedLine = bedLine + str(bedData[id]['thickStart']) + "\t"
    bedLine = bedLine + str(bedData[id]['thickEnd']) + "\t"
    bedLine = bedLine + bedData[id]['itemRgb'] + "\t"
    bedLine = bedLine + str(bedData[id]['blockCount']) + "\t"

    blockSizeLine   = ''
    blockStartsLine = ''

    for bstart in sorted(bedData[id]['blockStarts'].keys()):
        blockSizeLine = blockSizeLine + str(bedData[id]['blockStarts'][bstart]) + ","
        blockStartsLine = blockStartsLine + str(bstart - bedData[id]['start']) + ","
    
    bedLine = bedLine + blockSizeLine[:-1] + "\t" + blockStartsLine[:-1] + "\n"

    return bedLine


# Check --rgb option
if args.rgb != None:
    rgbL=args.rgb.split(sep=',')
    if len(rgbL) != 3: 
        print("ERROR: Expecting 3 integers between 0 and 255. Check --rgb format please")
        sys.exit()
    
    try:
        int(rgbL[0])
    except:
        print("ERROR: First element seems to not be an integer. Please correct --rgb option.")
        sys.exit()
    
    try:
        int(rgbL[1])
    except:
        print("ERROR: Second element seems to not be an integer. Please correct --rgb option.")
        sys.exit()
    
    try:
        int(rgbL[2])
    except:
        print("ERROR: Third element seems to not be an integer. Please correct --rgb option.")
        sys.exit()

    r = range(-1,256)
    if (int(rgbL[0]) not in r):
        print("ERROR: " + rgbL[0] + " not in range 0-255. Please correct --rgb option.")
        sys.exit()

    if (int(rgbL[1]) not in r):
        print("ERROR: " + rgbL[1] + " not in range 0-255. Please correct --rgb option.")
        sys.exit()

    if (int(rgbL[2]) not in r):
        print("ERROR: " + rgbL[2] + " not in range 0-255. Please correct --rgb option.")
        sys.exit()
else: 
    args.rgb = '0,100,255'

# Set output file name
bed = ''
features2keep = ['gene', 'mRNA', 'exon', 'CDS', 'transcript']

if args.bed == None: 
    bed = re.sub(r'(GTF|gtf|GFF3|gff3|gff)$', 'bed', args.gtf.name)
else: 
    bed = args.bed.name()


print("Running gtf2bed.py version ", ver)
print("==> New bed file: ", bed)
print("==> Reading GTF file")
gtf_file = open(args.gtf.name, mode='r')
bedData = {}

for line in gtf_file:
    if re.match('^#', line): # Skip comments
        continue

    fields=line.rstrip('\r\n').split(sep='\t')

    seq       = fields[0]
    source    = fields[1]
    feature   = fields[2]
    start     = int(fields[3])
    end       = int(fields[4])
    score     = fields[5]
    strand    = fields[6]
    frame     = fields[7]
    attribute = fields[8]

    if feature not in features2keep:
        continue

    if feature == 'gene' and args.inc:
        res=re.search('gene_id "([^;]+)";', fields[8])
        id=res.group(1)

        initiateAnnotation(id)
        gene2bed(id)

        continue

    if feature == 'mRNA' or feature == 'transcript':
        res=re.search('transcript_id "([^;]+)";', fields[8])
        id=res.group(1)

        initiateAnnotation(id)
        updateTranscript(id)

        continue

    if feature == 'exon' or feature == 'CDS':
        res=re.search('transcript_id "([^;]+)";', fields[8])
        parent=res.group(1)

        if feature == 'exon':
            addExon(parent)
            continue

        if feature == 'CDS':
            addCDS(parent)
            continue

print("==> Sort annotations by start")
bedDataKeys=sorted(
    bedData, 
    key=lambda x: (
        bedData[x]['seq'], 
        bedData[x]['start']))

bedFile=open(bed, 'w')

print("==> Writing BED12 file")
for key in bedDataKeys:
    bedFile.write(getBed12Line(key))

print("==> Done")

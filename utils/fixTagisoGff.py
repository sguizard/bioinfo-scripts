#! /usr/bin/env python

import sys
import argparse
import re

# Set version and licence
ver='1.0.0'
gpl='fixTagisoGff.py version ' + ver
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
    prog='fixTagisoGff.py', 
    description='Fix Tagiso GFF file to comply with GFF standard.', 
    usage='fixTagisoGff.py --gff tagiso.gtf', 
    epilog=gpl, 
    formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument('-g', '--gff', help='Tagsio GFF file to be fixed', type=argparse.FileType(mode='r'), required=True)
parser.add_argument('-o', '--out', help='Name of Tagiso GFF', type=str, required=False)

args = parser.parse_args()

# Define functions
def initiateGene(gid):
    if gid not in geneData:
        print("initiateGene: Adding " + gid)
        geneData[gid]               = {}
        geneData[gid]['transcript'] = {}


def initiateTranscript(gid, tid):
    if tid not in geneData[gid]['transcript']:
        print("initiateTranscript: Adding " + tid)
        geneData[gid]['transcript'][tid]         = {}
        geneData[gid]['transcript'][tid]['exon'] = {}


def addGene(id, line, seq, start): 
    # print("addGene: Adding " + id)

    initiateGene(id)

    toAdd = '\tID=' + id + ';gene_id'
    line  = re.sub(r'\tgene_id', toAdd, line)

    geneData[id]['line']  = line
    geneData[id]['seq']   = seq
    geneData[id]['start'] = start


def addTranscript(id, line, start): 
    att_fields = attribute.rstrip('\r\n').split(sep=';')
    tid = re.search('transcript_id=([^=;]+)', att_fields[1]).group(1)  # transcript/exon records
    # print("addTranscript: Adding " + tid)

    initiateGene(id)
    initiateTranscript(id, tid)

    toAdd = '\tID=' + tid + ';Parent='+ id + ';gene_id'
    line  = re.sub(r'\tgene_id'     , toAdd     , line)
    line  = re.sub(r'\ttranscript\t', '\tmRNA\t', line)

    geneData[id]['transcript'][tid]['line']  = line
    geneData[id]['transcript'][tid]['start'] = start


def addExon(id, line, start, end): 
    att_fields = attribute.rstrip('\r\n').split(sep=';')
    tid = re.search('transcript_id=([^=;]+)', att_fields[1]).group(1)
    eid = tid + "_exon_" + str(start) + "_" + str(end)
    # print("addExon: Adding " + eid)

    initiateGene(id)
    initiateTranscript(id, tid)

    toAdd = '\tID=' + eid + ';Parent='+ tid + ';gene_id'
    line = re.sub(r'\tgene_id', toAdd, line)

    geneData[id]['transcript'][tid]['exon'][eid]          = {}
    geneData[id]['transcript'][tid]['exon'][eid]['line']  = line
    geneData[id]['transcript'][tid]['exon'][eid]['start'] = start


def getLines(id):
    lines = geneData[id]['line']

    transcriptKeys = sorted(
        geneData[id]['transcript'], 
        key=lambda x: geneData[id]['transcript'][x]['start'])
    
    for t in transcriptKeys:
        lines += geneData[id]['transcript'][t]['line']

        exonKeys = sorted(
            geneData[id]['transcript'][t]['exon'], 
            key=lambda x: geneData[id]['transcript'][t]['exon'][x]['start'])
        
        for e in exonKeys:
            lines += geneData[id]['transcript'][t]['exon'][e]['line']

    return lines



# Set output file name
fixed_gff = ''
features2keep = ['gene', 'transcript', 'exon']

if args.out:
    fixed_gff = args.out
else:
    fixed_gff = re.sub(r'.(GFF3|gff3|GFF|gff|gtf|GTF)$', '_fixed.gff3', args.gff.name)


print("Running fixTagisoGff.py version ", ver)
print("==> New Tagiso file: ", fixed_gff)
print("==> Reading GFF file")
gff_file = open(args.gff.name, mode='r')

geneData       = {}

for line in gff_file:
    if re.match('^#', line): # Skip comments
        continue

    line = re.sub(r'; ', ';', line)
    line = re.sub(r' ' , '=', line)
    line = re.sub(r'"' , '' , line)

    fields = line.rstrip('\r\n').split(sep='\t')

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

    res=re.search('gene_id=([^=;]+);', fields[8])
    id=res.group(1)

    if feature == 'exon':
        addExon(id, line, start, end)
        continue

    if feature == 'transcript':
        addTranscript(id, line, start)
        continue

    if feature == 'gene':
        addGene(id, line, seq, start)
        continue

print("==> Sort genes by start")
geneDataKeys = sorted(
    geneData, 
    key = lambda x: (
        geneData[x]['seq'], 
        geneData[x]['start']) )

outFile = open(fixed_gff, 'w')

print("==> Writing Fixed GFF")
for key in geneDataKeys:
    outFile.write(getLines(key))

print("==> Done")

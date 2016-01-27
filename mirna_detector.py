#!/usr/bin/python
# -*- coding: utf-8 -*-
import sys
import re
import time 
from random import shuffle

expression = '((hsa|mmu|\b)?-?(miRNA-|miR|\(miR\)|micoRNA|hsa-let|let-|microRNA-|micro ribonucleic acid)(-|\s|\d)?((-|\w|\*|\/)*\d+(-|\w|\/)*)+\*?)'
extraction = '(hsa-(mir|let)-\d+(-|\w|\/)*)+'

mirnas = []
with open('data/mirna.txt') as f:

    content = f.readlines()
    for line in content:
        values = line.split('\t')
        mirna = values[2]
        if mirna.find('hsa') != -1:
            mirnas.append(mirna)
            values = values[3].split(';')
            if values[0]!='':
                mirnas+=values


def refine(results):
    parsed = []
    for result in results:
        value = re.search(expression, result, flags=re.I)
        if value:
            parsed.append(value.group())
    
    return parsed


def filterMirnas(results):
    global mirnas

    return [result for result in results if (result in mirnas or re.search(extraction, result, flags=re.I))]


def normalize(candidate):

    #---------
    if not re.search(expression, candidate, flags=re.I):
        return (candidate, '', '')
    #---------
    

    value = re.search(r'\d+', candidate)
    if not value:
        return (candidate,'','')

    if candidate[0] == '-':
        return (candidate, '', '')

    base = value.group()
    namePrefix = 'hsa-mir-'
    if re.search('let', candidate, re.I):
        namePrefix = 'hsa-let-'

    baseIndex = candidate.find(base) + len(base)
    nameSufix = candidate[baseIndex:]

    return (namePrefix, base, nameSufix)


def splitter(results):
    results = [result for result in results.split('/') if len(result)>0]
    # extract the mirna number base aka mir-(\d+)
    namePrefix, base, nameSufix = normalize(results[0])
    splitted = []

    splitted.append(namePrefix+base+nameSufix)

    for result in results[1:]:
        
        if not result[0].isdigit() and len(result)==1:
            splitted.append(namePrefix+base+result)
        else:
            np, b, ns = normalize(result)
            if np == result:
                value = namePrefix+result
                value = value.replace('--', '-')
                splitted.append(value)
            else:
                splitted.append(np+b+ns)
    return refine(splitted)


def expand(sentence, result, limit):
    value = result[0]
    
    if sentence[result[2]] == ',' or sentence[result[2]:result[2]+5]==' and ' if len(sentence)>result[2]+5 else False:
        expanded = sentence[result[2]:limit].replace(', ','/')
        expanded = re.sub(r'\s?and ', '/', expanded)
        spaceIndex = expanded.find(' ')
        value += expanded        
    return value


def validate(sentence):
    detected = {"sentence":sentence, "detectedMirnas":[]}

    if type(sentence) == list:
        sentence = ' '.join(sentence)
    p = re.compile(expression, flags=re.I|re.M)

    results = []
    for m in p.finditer(sentence):
        results.append((m.group(), m.start(), m.end()))

    parsedResults = []
    for index, result in enumerate(results):
        lastIndex = 0
        if index+1<len(results):
            lastIndex = results[index+1][1]
        else:
            lastIndex = len(sentence)
        extractValue = expand(sentence, result, lastIndex)

        values = splitter(extractValue)
        values = filterMirnas(values)
        if len(values)>0:
            for miRNA in values:
                detected["detectedMirnas"].append({
                    "mirna": miRNA,
                    "origin": extractValue
                })
        parsedResults += values

    return detected


def simpleTest(sentence):
    print validate('miR-26a/b')
    print validate('let-7')
    print validate('hsa-mir-33')
    print validate('hsa')
    print validate('Hepatocytes contain abundant miR-122 and damage of hepatocytes caused by inflammation due to virus infection or cancer would be expected to release significant amount of this miRNA into the circulation.')
    results = validate("""
Perhaps the most interesting result was the upregulation of let-7f upon treatment with letrozole, as it inhibits the expression of the aromatase gene ().
In a study to identify estrogen-regulated miRNAs, a group of miRNAs was shown to have elevated expression upon combined treatment of tamoxifen and exemestane, including miR-21, miR-181b, miR-26a/b, miR-27b, and miR-23b.
Expression of miR-101 is thus important for estrogen-independent growth mediated by the PI3K/Akt pathway, which could be important for resistance to both tamoxifen and AIs.
Upon analysis of clinical data from The Cancer Genome Atlas, it was clear that let-7c and HER2 were inversely correlated in luminal A breast tumors and that low expression of the let-7c/miR-99a/miR-125b cluster was associated with worse overall survival compared with patients who had high expression of this cluster.
When compared with LTEDaro cells and cells treated only with testosterone, letrozole-resistant cells overexpressed miR-128a suggesting miR-128a plays a role in AI resistance.
all my microRNAs hhv6b-miR-Ro6-4-3p and miR-155
""")
    


def bigTest():
    sentences = ''
    with open('tests/sentences.txt') as thefile:
        sentences = thefile.readlines()

    shuffle(sentences)
    for sentence in sentences:
        print '\n' + '*'*50
        print sentence
        print '-    results: ',
        print validate(sentence)
        time.sleep(5)


def read_input(file):
    for sentence in file:
        yield sentence


def main():
    for sentence in read_input(sys.stdin):
        print validate(sentence),
        print ',',


if __name__ == '__main__':
    main()

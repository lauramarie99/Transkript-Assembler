#!/usr/bin/env python

import sys, ast, os
import networkx as nx
from collections import namedtuple

from sqlalchemy import false
from parse_graph_list_commented_Arbeitsdatei import nodepath_to_transcript
from PathEnumeration import activeBinPathEnumeration3

def groupPairedBins(pairedBins):
    PairedBinT = namedtuple('GroupPairedBinT', 'leftExons rightExons')
    groupedPairedBins = []
    successor = []
    rightExonList = []
    counter = 0
    for bin in pairedBins:
        if counter < 1:
            successor = bin.leftExons
            rightExonList.append(bin.rightExons)
            counter = counter + 1
            continue
        if successor == bin.leftExons:
            rightExonList.append(bin.rightExons)
        else:
            groupedPairedBins.append(PairedBinT(leftExons=successor, rightExons=rightExonList))
            successor = bin.leftExons
            rightExonList = [bin.rightExons]
    return groupedPairedBins
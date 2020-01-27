import numpy as np
import copy
import random

def printList (list1):

    for l in list1:
        print(l)


def getAllKmersFromSingleDNAString (text, k): # verified
    kmers = list()
    for i in range(len(text) - k + 1):
        kmer1 = text[i:i+k]
        if not kmer1 in kmers:
            kmers.append(kmer1)
    return kmers

def getAllKmersFromSingleString_withDupes (text, k):
    kmers = list()
    for i in range(len(text) - k + 1):
        kmer1 = text[i:i + k]
        kmers.append(kmer1)
    return kmers


def composition( text, k ):
    kmers = getAllKmersFromSingleString_withDupes (text, k)
    kmers.sort()
    return kmers

def getPrefix (pattern):
    prefix = pattern[:len(pattern)-1]
    return prefix


def getSuffix (text):
    suffix = text[1:len(text)]
    return suffix


def getLastChar (text):
    suffix = text[len(text)-1]
    return suffix


def reconstructStringFromPreOrderedKmers (reads):
    outString = ""
    for i in range (len(reads)):
        #print(i)
        if i == 0:
            outString = outString + reads[i]
        else:
            lastChar = getLastChar (reads[i])
            #print(lastChar)
            outString = outString + lastChar
    return outString

def makeGraphString(string1, string2):
    outString = string1 + " -> " + string2
    return outString


def getKMERObjectsList (kmers):
    kmersObjs = list()

    for text in kmers:
        prefix = getPrefix (text)
        suffix = getSuffix (text)
        km = Kmer(text, prefix, suffix)
        kmersObjs.append(km)
    return kmersObjs

def getKmerObjectsFromSingleString (text, k):
    kmers = getAllKmersFromSingleDNAString(text, k)
    kmersObjs = getKMERObjectsList(kmers)
    return kmersObjs

def findOverlapGraphFromStrings (texts):

    kmers = list()

    for text in texts:
        prefix = getPrefix (text)
        suffix = getSuffix (text)
        km = Kmer(text, prefix, suffix)
        kmers.append(km)

    graphStrings = list()
    for kmer1 in kmers:
        for kmer2 in kmers:
            if kmer1.seq != kmer2.seq:
                if kmer1.prefix == kmer2.suffix:
                    graphString = makeGraphString(kmer2.seq, kmer1.seq)
                    if not graphString in graphStrings:
                        graphStrings.append(graphString)
                if kmer2.prefix == kmer1.suffix:
                    graphString = makeGraphString(kmer1.seq, kmer2.seq)
                    if not graphString in graphStrings:
                        graphStrings.append(graphString)
    graphStrings.sort()
    return graphStrings




def addToDebrujinGraphDict (debrujinDict, inPrefix, outPrefix):
    if not inPrefix in debrujinDict:
        outy_List = list()
        outy_List.append(outPrefix)
        debrujinDict[inPrefix] = outy_List
    else:
        outy_List = debrujinDict[inPrefix]
        if not outPrefix in outy_List:
            outy_List.append(outPrefix)
        debrujinDict[inPrefix] = outy_List
    return debrujinDict


def addToDebrujinGraphDict_WithDuplicateOutPrefix (debrujinDict, inPrefix, outPrefix):
    if not inPrefix in debrujinDict:
        outy_List = list()
        outy_List.append(outPrefix)
        debrujinDict[inPrefix] = outy_List
    else:
        outy_List = debrujinDict[inPrefix]
        outy_List.append(outPrefix)
        debrujinDict[inPrefix] = outy_List
    return debrujinDict


def get_InKmerOutKmerObjs_From_KmerObjs (kmersObjs):
    inKmersAndoutKmers = list()
    for kmer1 in kmersObjs:  # match all graph edges
        for kmer2 in kmersObjs:
            if kmer1.seq != kmer2.seq:
                if kmer1.prefix == kmer2.suffix:
                    inny_outy = InKmerOutKmer(kmer2.seq, kmer1.seq, kmer2.prefix, kmer1.prefix)
                    inKmersAndoutKmers.append(inny_outy)
                if kmer2.prefix == kmer1.suffix:
                    inny_outy = InKmerOutKmer(kmer1.seq, kmer2.seq, kmer1.prefix, kmer2.prefix)
                    inKmersAndoutKmers.append(inny_outy)
    for kmer in kmersObjs:
        if kmer.prefix == kmer.suffix:
            inny_outy = InKmerOutKmer(kmer.seq, kmer.seq, kmer.prefix, kmer.prefix)
            inKmersAndoutKmers.append(inny_outy)
    return inKmersAndoutKmers

#
#
# def deduplicateAndSortOutPrefix (debrujinDict):
#     for key1 in debrujinDict.keys():
#         outList2 = list()
#         outy_List = debrujinDict[key1]
#         for o in outy_List:
#             if not o in outList2:
#                 outList2.append(o)
#         outList2.sort()
#         debrujinDict[key1] = outList2
#     return debrujinDict
#


def getDeBrujinGraphFromString (text, k):
    kmersObjs = getKmerObjectsFromSingleString(text, k)
    inKmersAndoutKmers = get_InKmerOutKmerObjs_From_KmerObjs(kmersObjs)
    inPrefixDict = dict()
    for edge in inKmersAndoutKmers:     # collapse duplicate in kmer prefixes
        inPrefixDict = addToDebrujinGraphDict(inPrefixDict, edge.inPrefix, edge.outPrefix)
    # woopsy. get last kmer prefix and suffix
    lastKmer = text[-k:]
    prefix = getPrefix(lastKmer)
    suffix = getSuffix(lastKmer)
    inPrefixDict = addToDebrujinGraphDict(inPrefixDict, prefix, suffix)
    return inPrefixDict

def makeStringFromList_withComma (list1):
    outString = ""
    for l in list1:
        outString = outString + l + ","
    if outString.endswith(","):
        outString = outString[:len(outString)-1]
    return outString


def makeStringsFromInOutObjDict (inPrefixDict):
    outstrings = list()
    for key1 in inPrefixDict.keys():
        outy_List = inPrefixDict[key1]
        outyString = makeStringFromList_withComma(outy_List)
        oString = makeGraphString(key1, outyString)
        outstrings.append(oString)
    return outstrings

def makeStringsFromInOutObjDict_WithSort (inPrefixDict):
    outstrings = list()
    for key1 in inPrefixDict.keys():
        outy_List = inPrefixDict[key1]
        outy_List.sort()
        outyString = makeStringFromList_withComma(outy_List)
        oString = makeGraphString(key1, outyString)
        outstrings.append(oString)
    return outstrings




def makeDeBrujinGraphFromKmers (kmers):
    deBrujinDict = dict()
    for kmer in kmers:
        prefix = getPrefix(kmer)
        suffix = getSuffix(kmer)
        deBrujinDict = addToDebrujinGraphDict_WithDuplicateOutPrefix(deBrujinDict, prefix, suffix)
    return deBrujinDict




def parseGraphLines_toDictAndList (lines):

    totalEdgeCount = 0
    graphList = list()

    for line in lines:
        l = line.split(" -> ")
        inEdge = l[0].strip()
        outEdge = l[1].strip()
        # 24 -> 22,34,820,927
        endNodes = list()
        if not "," in outEdge:
           string1 = inEdge + " -> " + outEdge
           graphList.append(string1)
           endNodes.append(outEdge)
           totalEdgeCount = totalEdgeCount + 1
        else:
            #print(line)
            oe = outEdge.split(',')
            #print(oe)
            for o in oe:
                endNodes.append(o)
                string1 = inEdge + " -> " + o
                #print(string1)
                graphList.append(string1)
                totalEdgeCount = totalEdgeCount + 1

    graph = dict()
    for g in graphList:
        startNode, endNode = getStartAndEndNodes(g)
        if startNode in graph:
            endNodes = graph[startNode]
            endNodes.append(endNode)
            graph[startNode] = endNodes
        else:
            endNodes = list()
            endNodes.append(endNode)
            graph[startNode] = endNodes
    return graph, totalEdgeCount, graphList


def eulerCyclerDumberer (graphList, graphDict):
    hasUntravelled = True

    unusedStartEdges = graphList.copy()
    unusedStartDict = copy.deepcopy(graphDict)

    #startEdge = bestFirstStart (graphDict, unusedStartEdges)
    startEdge = random.choice(unusedStartEdges)

    startNode, endNode = getStartAndEndNodes (startEdge)

    bestPath =list()
    unusedStartEdges, unusedStartDict = removeNewStartEdgeFromLists(startEdge, unusedStartEdges, unusedStartDict)

    count = 0
    path = list()
    #print("initial start edge = " + startEdge)
    while hasUntravelled:
        #print("\n\n NEW CYCLE")
        random.seed(count)
        count = count + 1
        path = list()
        untravelledEdges = copy.deepcopy(graphList)
        untravelledDict = copy.deepcopy(graphDict)
        path.append(startEdge)
        atEnd = False
        while not atEnd:
            if startEdge in untravelledEdges:
                untravelledEdges.remove(startEdge)

            startNode, endNode = getStartAndEndNodes(startEdge)
            isZero, nextEdge =  getNextEdge (untravelledDict, endNode)
            if not isZero:
                #print("currentEdge = " + startEdge + "\tnextEdge = " + nextEdge + "\t\tcurrent end node = " + endNode)

                nextStart, nextEnd = getStartAndEndNodes(nextEdge)

                untravelledDict = adjustTravelledDict (endNode, nextEnd, untravelledDict)
                if nextEdge in untravelledEdges:
                    untravelledEdges.remove(nextEdge)
                if startEdge in untravelledEdges:
                    untravelledEdges.remove(startEdge)

                path.append(nextEdge)
                startEdge = nextEdge

            else:
                atEnd = True
        if len(path) > len(graphList):
            hasUntravelled = False
        else:
            if len(path) > len(bestPath):
                bestPath = path
            #print("path = " + str(path))
            #print("len path = " + str(len(path)) + "\tlen graph = " + str(len(graphList)) + "\t len best path = " + str(len(bestPath)))


            # find new node, start again
            #startEdge = bestFirstStart(graphDict, unusedStartEdges)
            if len(unusedStartEdges) > 0:
                startEdge = random.choice(unusedStartEdges)
            else:
                unusedStartEdges = graphList.copy()
                unusedStartDict = copy.deepcopy(graphDict)
                startEdge = random.choice(unusedStartEdges)
            startNode, endNode = getStartAndEndNodes(startEdge)
            unusedStartEdges, unusedStartDict = removeNewStartEdgeFromLists(startEdge, unusedStartEdges, unusedStartDict)






    print(path)
    path = path[:len(path) - 1]
    print(path)
    pathString = pathStringMaker_forContigOverlap (path)
    return pathString




def pathStringMaker_forContigOverlap (path):

    pathString = ""

    pathString = path[0] + " -> "

    for i in range(len(path)):
        if i > 0:

            startNode, endNode = getStartAndEndNodes (path[i])
            pathString = pathString + endNode + " -> "

    if pathString.endswith(" -> "):
        pathString = pathString[:len(pathString)-4].strip()
    return pathString


def removeNewStartEdgeFromLists (startEdge, unusedList, unusedDict):

    if startEdge in unusedList:
        unusedList.remove(startEdge)


    startNode, endNode = getStartAndEndNodes(startEdge)

    endNodes = unusedDict[startNode]
    #print("startEdge = " + startEdge + "\tstart node = " + str(startNode) + "\tend node = " + str(endNode) + "\t endNodes = " + str(endNodes))
    endNodes.remove(endNode)

    unusedDict[startNode] = endNodes

    #print("new endNodes= " + str(unusedDict[startNode]))

    return unusedList, unusedDict



def adjustTravelledDict (startNode, endNode, untravelledDict):
    endNodes = untravelledDict[startNode]
    endNodes.remove(str(endNode))
    untravelledDict[startNode] = endNodes
    return untravelledDict


def getValidEdgeFromList (endNode, unusedStartEdges):
    usableStartEdges = list()
    for edge in unusedStartEdges:
        startNode = getStartAndEndNodes(edge)
        if startNode == endNode:
            usableStartEdges.append(edge)
    if len(usableStartEdges) > 0:
        startEdge = random.choice(usableStartEdges)
    else:
        startEdge = random.choice(unusedStartEdges)
    return startEdge






def getNextEdge (graphDict, endNode):
    isZero = False
    nextEdge = ""
    nextNodes = graphDict[endNode]
    if len(nextNodes) == 0:
        isZero = True
    else:
        nextEnd = random.choice(nextNodes)
        nextEdge = endNode + " -> " + nextEnd
    return isZero, nextEdge



def getStartAndEndNodes (edge):

    e = edge.split("->")

    startNode = e[0].strip()
    endNode = e[len(e)-1].strip()
    return startNode, endNode






############################# File Parsers ###############################################


def parse21_file (filename):
    k = -1
    text = ""
    file = open(filename)
    i = 0
    for line in file:
        if i == 0:
            k = int(line.strip())
        if i  == 1:
            text = line.strip()
        i = i + 1
    file.close()
    return k, text

def parse22_file (filename):
    texts = list()
    file = open(filename)
    for line in file:
        if len(line.strip()) > 0:
            texts.append(line.strip())
    return texts



########################### output ##############################################

def output22_file (filename, outString):
    file = open(filename, 'w')
    file.write(outString)
    file.close()

def output23_file (filename, outlist):
    file = open(filename, 'w')
    for outString in outlist:
        file.write(outString + "\n")
    file.close()


###################### Main things ###############################################



def doProblem_25():
    kmers = parse22_file("/home/jacklyn/Downloads/rosalind_ba3e.txt")
    debrujingraph = makeDeBrujinGraphFromKmers (kmers)
    outstrings = makeStringsFromInOutObjDict_WithSort (debrujingraph)
    output23_file("/home/jacklyn/PycharmProjects/CSE282/HW1/25.txt", outstrings)
    print("done!")

def doProblem_26():
    lines = parse22_file ("/home/jacklyn/Downloads/rosalind_ba3f.txt")
    graph, totalEdgeCount, graphList = parseGraphLines_toDictAndList (lines)

    print(graph)
    print(graphList)


    eulerPath = eulerCyclerDumberer (graphList, graph)
    print(eulerPath)
    output22_file('/home/jacklyn/PycharmProjects/CSE282/HW1/26.txt', eulerPath)
    print("done!")



doProblem_26()

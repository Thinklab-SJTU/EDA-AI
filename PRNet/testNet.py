import MST

def parsenet(filepath):
    f=open(filepath,encoding='utf-8')
    a=f.readline()[2:-2].split('], [')
    netlist=[]
    for netstring in a:
        net=[int(index) for index in netstring.split(', ')]
        netlist.append(net)
    f.close()
    return netlist

def parseresult(filepath):
    f = open(filepath, 'r')
    for line in f:
        result = eval(line)
    for item in result:
        item.append(0)
    f.close()
    return result


def getpinpairs_f(netfile,resultfile):
    nets = parsenet(netfile)
    coords = parseresult(resultfile)
    # print(nets)
    # print(coords)
    netcoords = []
    for net in nets:
        allpins = [coords[i] for i in net]
        twopins = []
        for i in range(len(allpins)):
            for j in range(i + 1, len(allpins)):
                twopins.append([tuple(allpins[i]), tuple(allpins[j])])
        mst = MST.generateMST(twopins)
        netcoords.append(mst)

    return netcoords

def getpinpairs():
    nets=getpairs()
    coords=getcoords()
    #print(nets)
    #print(coords)
    netcoords=[]
    for net in nets:
        allpins=[coords[i] for i in net]
        twopins=[]
        for i in range(len(allpins)):
            for j in range(i+1,len(allpins)):
                twopins.append([tuple(allpins[i]),tuple(allpins[j])])
        mst=MST.generateMST(twopins)
        netcoords.append(mst)

    return netcoords

def allpinpairs(pinpairs):
    allpairs=[]
    netnum=[]
    for net in pinpairs:
        netnum.append(len(net))
        for pinpair in net:
            allpairs.append(pinpair)
    return allpairs,netnum
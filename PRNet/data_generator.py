import os
import random
class Generator:
    def __init__(self):
        self.rootdir=''
        self.path=''

    def setpath(self,root='',path=''):
        self.rootdir=root
        self.path=path

    def generate(self,parameters):
        pass


def dic_generate(width,height,layer,vcap,hcap,minwidth,minspace,viaspace,twidth,theight,redcap,numnet,pinnumub):
    parameter={}
    parameter["gridSize"]=[width,height,layer]
    parameter["verticalCapacity"]=[0.0,vcap,0.0]
    parameter["horizontalCapacity"]=[hcap,0.0,hcap]
    parameter["minWidth"]=minwidth
    parameter["minSpacing"]=minspace
    parameter["viaSpacing"]=viaspace
    parameter["Origin"]=[0.0,0.0]
    parameter["tileWidth"]=twidth
    parameter["tileHeight"]=theight
    parameter["reducedCapacity"]=redcap
    parameter["numNet"]=numnet

    netinfo=[]
    for i in range(numnet):
        net={}
        net["netName"]="A%d"%(i+1)
        net["netID"]=i+1
        pinnum=random.randrange(2,pinnumub)
        net["numPins"]=pinnum
        net["minWidth"]=1.0
        for j in range(pinnum):
            pinx=random.randrange(width*twidth)
            piny=random.randrange(height*theight)
            pinz=0
            net[str(j+1)]=[pinx,piny,pinz]
        netinfo.append(net)

    parameter["netInfo"]=netinfo
    return parameter

def dic_generate_2(width,height,layer,vcap,hcap,minwidth,minspace,viaspace,twidth,theight,redcap):
    parameter={}
    parameter["gridSize"]=[width,height,layer]
    parameter["verticalCapacity"]=[vcap]
    parameter["horizontalCapacity"]=[0.0]
    parameter["minWidth"]=minwidth
    parameter["minSpacing"]=minspace
    parameter["viaSpacing"]=viaspace
    parameter["Origin"]=[0.0,0.0]
    parameter["tileWidth"]=twidth
    parameter["tileHeight"]=theight
    parameter["reducedCapacity"]=redcap

    pv=parameter["verticalCapacity"]
    for i in range(layer-1):
        if i%2==0:
            pv.append(0.0)
        else:
            pv.append(vcap)
    ph=parameter["horizontalCapacity"]
    for i in range(layer-1):
        if i%2==0:
            ph.append(hcap)
        else:
            ph.append(0.0)

    return parameter

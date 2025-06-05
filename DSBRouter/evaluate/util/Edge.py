from util.point import *
class Edge():
    def __init__(self, point1:tuple, point2:tuple,of,wl):
        self.point1 = point1
        self.point2 = point2
        self.of = of
        self.wl = wl
        if self.point1[0] == self.point2[0]:
            self.direction = 'V'
        else:
            self.direction = 'H'

    def GetEdge(self):
        pass
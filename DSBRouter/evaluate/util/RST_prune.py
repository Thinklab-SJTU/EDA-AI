import torch
import concurrent.futures
import time
from util.Edge import *
import numpy as np
from typing import Dict, Set, List
from itertools import combinations
import heapq
'''
Algorithm prototype
'''
class LType():
    def __init__(self,point1:tuple,steinerpoint:tuple,point2:tuple,tensor_:torch.Tensor,alpha,beta)->None:
        assert \
        (point1[0] == steinerpoint[0] or \
        point1[1] == steinerpoint[1] or \
        point2[0] == steinerpoint[0] or \
        point2[1] == steinerpoint[1]), "Steiner point must align with point1 or point2 either horizontally or vertically."

        C, H, W = tensor_.shape
        tensor_np = tensor_.cpu().numpy()
        channel_R = tensor_np[0]
        # horizontal capacity
        channel_G = tensor_np[1]
        # vertical capacity
        channel_B = tensor_np[2]
        if point1[0] == steinerpoint[0]:
            y1, y2 = min(point1[1], steinerpoint[1]), max(point1[1], steinerpoint[1])
            x = steinerpoint[0]
            of = channel_B[x, y1:y2].sum()-1
            wl = y2-y1
            self.edge1 = Edge(point1,steinerpoint,-of,wl)
        elif point1[1] == steinerpoint[1]:
            x1, x2 = min(point1[0], steinerpoint[0]), max(point1[0], steinerpoint[0])
            y = steinerpoint[1]
            of = channel_G[x1:x2, y].sum()-1
            wl = x2-x1
            self.edge1 = Edge(point1,steinerpoint,-of,wl)
        else:
            raise ValueError("Invalid alignment between point1 and steinerpoint.")
        if steinerpoint[0] == point2[0]:
            y1, y2 = min(steinerpoint[1], point2[1]), max(steinerpoint[1], point2[1])
            x = steinerpoint[0]
            of = channel_B[x, y1:y2].sum()-1
            wl = y2-y1
            self.edge2 = Edge(steinerpoint,point2,-of,wl)
        elif steinerpoint[1] == point2[1]:
            x1, x2 = min(steinerpoint[0], point2[0]), max(steinerpoint[0], point2[0])
            y = steinerpoint[1]
            of = channel_G[x1:x2, y].sum()-1
            wl = x2-x1
            self.edge2 = Edge(steinerpoint,point2,-of,wl)
        else:
            raise ValueError("Invalid alignment between steinerpoint and point2.")

        self.point1 = point1
        self.steinerpoint = steinerpoint
        self.point2 = point2
        self.total_of = self.edge1.of + self.edge2.of
        self.total_wl = self.edge1.wl + self.edge2.wl
        self.avg_of = self.total_of / self.total_wl
        if alpha is not None and beta is not None:
            self.alpha = alpha
            self.beta = beta
            self.score = self.alpha * self.total_of + self.beta * self.total_wl

    def __repr__(self):
        return (f"LType(edge1={self.edge1}, edge2={self.edge2}, "
                f"total_of={self.total_of}, total_wl={self.total_wl})")

    def __eq__(self, other):
        if not isinstance(other, LType):
            return False
        return (
                {self.edge1.point1, self.edge1.point2, self.edge2.point2} ==
                {other.edge1.point1, other.edge1.point2, other.edge2.point2}
        )

    def __hash__(self):
        points = sorted([self.edge1.point1, self.edge1.point2, self.edge2.point2], key=lambda p: (p[0], p[1]))
        return hash(tuple(points))

class LTypeSet():
    def __init__(self,serial_number)->None:
        self.LTypes = {}
        self.serial_number = serial_number
        self.Ltypes = []

class EdgesSet():
    def __init__(self,serial_number,alpha=None,beta=None):
        self.serial_number = serial_number
        self.edges = {}
        self.edgeset = set()
        self.oftoegdes = {}
        self.wltoedges = {}
        self.scorestoedges = {}
        self.avgoftoedges = {}
        self.alpha = None
        self.beta = None
        if alpha is not None and beta is not None:
            self.alpha = alpha
            self.beta = beta

    def initialize_edges(self,tensor_:torch.Tensor)->None:
        C,H,W = tensor_.shape
        tensor_np = tensor_.cpu().numpy()
        channel_R = tensor_np[0]
        #horizontal capacity
        channel_G = tensor_np[1]
        #vertical capacity
        channel_B = tensor_np[2]
        # points = (channel_R >= 5 / 6).nonzero(as_tuple=True)
        # print(f'points:{points}')
        for x in range(W):
            y_indices = np.where(channel_R[x, :] >=2/3)[0]
            if len(y_indices) < 2:
                continue 
            # print(f'y_indices:{y_indices},x:{x}')
            for i in range(len(y_indices) - 1):
                y1, y2 = y_indices[i], y_indices[i + 1]
                p1 = (x, y1)
                p2 = (x, y2)
                key = tuple(sorted([p1, p2]))  
                if key not in self.edgeset:
                    self.edgeset.add(key)
                    of = channel_G[x, y1:y2].sum() - 1
                    wl = y2 - y1
                    edge = Edge(p1, p2, of=-of, wl=wl)
                    self.edges[key] = edge
                    self.oftoegdes[key] = -of
                    self.wltoedges[key] = wl
                    self.avgoftoedges[key] = -of / wl
                    if self.alpha is not None and self.beta is not None:
                        self.scorestoedges[key] = self.alpha * of + self.beta * wl

        for y in range(H):
            x_indices = np.where(channel_R[:, y] >=2/3)[0]
            # print(f'x_indices:{x_indices},y:{y}')
            if len(x_indices) < 2:
                continue 
            for i in range(len(x_indices) - 1):
                x1, x2 = x_indices[i], x_indices[i + 1]
                p1 = (x1, y)
                p2 = (x2, y)
                key = tuple(sorted([p1, p2]))
                if key not in self.edgeset:
                    self.edgeset.add(key)
                    of = channel_B[x1:x2, y].sum() - 1
                    wl = x2 - x1
                    edge = Edge(p1, p2, of=-of, wl=wl)
                    self.edges[key] = edge
                    self.oftoegdes[key] = -of
                    self.wltoedges[key] = wl
                    self.avgoftoedges[key] = -of / wl
                    if self.alpha is not None and self.beta is not None:
                        self.scorestoedges[key] = self.alpha * of + self.beta * wl


    def addEdge(self,edge:Edge)->None:
        p1 = (edge.point1[0],edge.point1[1])
        p2 = (edge.point2[0],edge.point2[1])
        key = tuple(sorted([p1,p2]))
        if key not in self.edgeset:
            self.edgeset.add(key)
            self.edges[key] = edge
            self.oftoegdes[key] = edge.of
            self.wltoedges[key] = edge.wl
            self.avgoftoedges[key] = edge.of/edge.wl
            if self.alpha is not None and self.beta is not None:
                self.scorestoedges[key] = self.alpha * edge.of + self.beta * edge.wl

class union_findSet():
    def __init__(self,serial_number,parent):
        self.serial_number = serial_number
        self.pointprototype = set()
        if parent is None:
            self.parent = {}
            self.parent_copy = {}
        else:
            self.parent = parent
            self.parent_copy = parent
        self.components = {}

    def union(self, x:tuple, y:tuple)->None: 
        rootX = self.find(x)
        rootY = self.find(y)
        if rootX != rootY:
            self.parent[rootX] = rootY
    def find(self, x:tuple)->tuple:
        if self.parent[x] != x:
            self.parent[x] = self.find(self.parent[x])
        return self.parent[x]

    def union_copy(self, x:tuple, y:tuple)->None:
        rootX = self.find_copy(x)
        rootY = self.find_copy(y)
        if rootX != rootY:
            self.parent_copy[rootX] = rootY
    def find_copy(self, x:tuple)->tuple:
        if self.parent_copy[x] != x:
            self.parent_copy[x] = self.find_copy(self.parent_copy[x])
        return self.parent_copy[x]

    def initialize_union_find(self,tensor_:torch.Tensor)->None:
        tensor__ = tensor_[0]
        for i in range(tensor__.shape[0]):
            for j in range(tensor__.shape[1]):
                if tensor__[i, j] >= 2/3:
                    # point = Point(i, j)
                    point = (i,j)
                    self.parent[point] = point 
                    self.parent_copy[point] = point
                    if (i,j) not in self.pointprototype:
                        self.pointprototype.add((i,j))
        return

    def AddPoint(self,x:int,y:int)->None:
        if (x,y) in self.pointprototype:
            return
        point = (x,y)
        self.parent[point] = point
        self.parent_copy[point] = point
        self.pointprototype.add((x,y))

    def merge_edges_from_edges_set(self, edge_set: EdgesSet) -> None:
        for edge in edge_set.edges.values():
            # print(f'point1:{edge.point1}')
            # print(f'point1:{edge.point2}')
            self.union(edge.point1, edge.point2)
        components: Dict[tuple, Set[tuple]] = {}
        for point in self.parent.keys():
            root = self.find(point)
            if root not in components:
                components[root] = set()
            components[root].add(point)
        self.components = components

class Slover():
    def __init__(self, args):
        self.args = args
        self.unionfindsets = {}
        self.edgesets = {}
        self.LTypes = {}
        self.conditions = {}
        self.mode = self.args['mode']
        self.alpha = self.args['alpha'] if self.mode == 'mix' else None
        self.beta = self.args['beta'] if self.mode == 'mix' else None
        self.selected_edges = torch.full((self.args['batch_size'], 3, 64, 64), -1.0,dtype=torch.float32).to("cpu")
        for num in range(self.args['batch_size']):
            self.unionfindsets[num] = union_findSet(num,None)
            self.edgesets[num] = EdgesSet(num,self.alpha,self.beta)
            self.conditions[num] = set()
            self.LTypes[num] = LTypeSet(num)
        self.PinMaps = None
        self.conditionMaps = None
        self.SteinerMaps = None

        self.final_edges = set()

    def SetTarget(self,bg,pins,conditions,steinerpoints)->None:
        if pins!=None:
            self.PinMaps = pins
            self.PinMaps = self.PinMaps.to("cpu")
        if conditions!=None:
            self.conditionMaps = conditions
            self.conditionMaps = self.conditionMaps.to("cpu")
            # print(f'conditionMaps:{bg.shape}')
            self.selected_edges[:,1,:,:] = bg[:,1,:,:]
            self.selected_edges[:,2,:,:] = bg[:,2,:,:]
        if steinerpoints!=None:
            self.SteinerMaps = steinerpoints

    def InitialPointsInitialize(self,id:int)->str:
        time_ = time.time()
        for i in range(self.conditionMaps[id].shape[0]):
            for j in range(self.conditionMaps[id].shape[1]):
                if self.conditionMaps[id][i][j] >= 2/3:
                    # print(f'conditionMaps:{(i,j)}')
                    if (i,j) not in self.conditions[id]:
                        self.conditions[id].add((i, j))
        # print(f'conditions:{self.conditions[id]}')
        time__ = time.time()
        return f'Initial Points initialized {id} done, cost:{time__-time_}'

    def UnionSetInitialize(self,id:int)->str:
        time_ = time.time()
        self.unionfindsets[id].initialize_union_find(self.PinMaps[id])
        time__ = time.time()
        return f'union {id} done, cost:{time__-time_}'

    def EdgesInitialize(self,id:int)->str:
        time_ = time.time()
        self.edgesets[id].initialize_edges(self.PinMaps[id])
        time__ = time.time()
        return f'edges initialize {id} done, cost:{time__ - time_}'

    def PointUnion(self,id:int)->str:
        time_ = time.time()
        self.unionfindsets[id].merge_edges_from_edges_set(self.edgesets[id])
        time__ = time.time()
        return f'point union {id} done, cost:{time__ - time_}'

    # def LTypesSetInitialize(self,id:int)->str:
    #     time_ = time.time()
    #     self.LTypes[id].initialize_LTypes(self.PinMaps[id])
    #     time__ = time.time()
    #     return f'point union {id} done, cost:{time__ - time_}'

    def SteinerPointAdd(self, id: int) -> None:
        """
        生成路由连接所有组件，直到所有点在同一个连通分量中。
        """
        assert self.PinMaps[id].shape == (3, 64, 64)
        R_channel = self.PinMaps[id][0].cpu().numpy()
        points = np.argwhere(R_channel >= 2/3)
        point_list = [(int(p[0]), int(p[1])) for p in points]
        # print(f'point_list:{point_list}')
        # unique_L_Types: Set[LType] = set()
        if len(self.unionfindsets[id].components) == 1:
            return
        while len(self.unionfindsets[id].components) > 1:
            # print(f'components:{self.unionfindsets[id].components}')
            # print(f'parent:{self.unionfindsets[id].parent}')
            for p1, p2 in combinations(point_list, 2):
                if p1[0] == p2[0] or p1[1] == p2[1]:
                    continue
                steinerpoint1 = (p1[0], p2[1])
                steinerpoint2 = (p2[0], p1[1])
                if (p1[0], p2[1]) not in self.unionfindsets[id].pointprototype:
                    root_p1 = self.unionfindsets[id].find(p1)
                    root_p2 = self.unionfindsets[id].find(p2)
                    if root_p1 != root_p2:
                        try:
                            L1 = LType(point1=p1, steinerpoint=steinerpoint1, point2=p2, tensor_=self.PinMaps[id], alpha=self.alpha, beta=self.beta)
                            # unique_L_Types.add(L1)
                            self.unionfindsets[id].AddPoint(steinerpoint1[0], steinerpoint1[1])
                            self.unionfindsets[id].union(p1, steinerpoint1)
                            self.unionfindsets[id].union(steinerpoint1, p2)

                            self.edgesets[id].addEdge(L1.edge1)
                            self.edgesets[id].addEdge(L1.edge2)

                            self.unionfindsets[id].components = {}
                            for point in self.unionfindsets[id].parent.keys():
                                root = self.unionfindsets[id].find(point)
                                if root not in self.unionfindsets[id].components:
                                    self.unionfindsets[id].components[root] = set()
                                self.unionfindsets[id].components[root].add(point)
                            self.unionfindsets[id].pointprototype.add((steinerpoint1[0], steinerpoint1[1]))

                        except ValueError:
                            pass

                if (p2[0], p1[1]) not in self.unionfindsets[id].pointprototype:
                    root_p1 = self.unionfindsets[id].find(p1)
                    root_p2 = self.unionfindsets[id].find(p2)
                    if root_p1 != root_p2:
                        try:
                            L2 = LType(point1=p1, steinerpoint=steinerpoint2, point2=p2, tensor_=self.PinMaps[id], alpha=self.alpha, beta=self.beta)

                            self.unionfindsets[id].AddPoint(steinerpoint2[0], steinerpoint2[1])
                            self.unionfindsets[id].union(p1, steinerpoint2)
                            self.unionfindsets[id].union(steinerpoint2, p2)

                            self.edgesets[id].addEdge(L2.edge1)
                            self.edgesets[id].addEdge(L2.edge2)

                            self.unionfindsets[id].components = {}
                            for point in self.unionfindsets[id].parent.keys():
                                root = self.unionfindsets[id].find(point)
                                if root not in self.unionfindsets[id].components:
                                    self.unionfindsets[id].components[root] = set()
                                self.unionfindsets[id].components[root].add(point)
                            self.unionfindsets[id].pointprototype.add((steinerpoint2[0], steinerpoint2[1]))

                        except ValueError:
                            pass

    def SteinerPointAdd_mix(self, id: int) -> None:
        assert self.PinMaps[id].shape == (3, 64, 64)
        R_channel = self.PinMaps[id][0].cpu().numpy() 

        points = np.argwhere(R_channel >= 2 / 3)
        point_list = [(int(p[0]), int(p[1])) for p in points]
        # print(f'point_list:{point_list}')
        # print(f'self.unionfindsets[id].components:{self.unionfindsets[id].components}')
        if len(self.unionfindsets[id].components) <= 1:
            return

        while len(self.unionfindsets[id].components) > 1:
            candidate_LTypes = set()
            components_list = list(self.unionfindsets[id].components.values())
            for idx1 in range(len(components_list)):
                for idx2 in range(idx1 + 1, len(components_list)):
                    comp1 = components_list[idx1]
                    comp2 = components_list[idx2]
                    for p3 in comp1:
                        for p4 in comp2:
                            if p3[0] == p4[0] or p3[1] == p4[1]:
                                continue
                            steiner_candidates = [
                                (p3[0], p4[1]),
                                (p4[0], p3[1])
                            ]

                            # 构造并验证LType
                            for steiner_point in steiner_candidates:
                                if steiner_point in self.unionfindsets[id].pointprototype:
                                    continue
                                try:
                                    L = LType(point1=p3, steinerpoint=steiner_point, point2=p4,
                                              tensor_=self.PinMaps[id], alpha=self.alpha, beta=self.beta)
                                    candidate_LTypes.add(L)
                                except ValueError:
                                    continue
            # print(f'candidate_LTypes:{candidate_LTypes}')
            if not candidate_LTypes:
                break
            if self.mode == 'avg_of':
                best_LType = min(candidate_LTypes, key=lambda lt: lt.avg_of)
            elif self.mode == 'of':
                best_LType = min(candidate_LTypes, key=lambda lt: lt.total_of)
            elif self.mode == 'wl':
                best_LType = min(candidate_LTypes, key=lambda lt: lt.total_wl)
            elif self.mode == 'mix':
                best_LType = min(candidate_LTypes, key=lambda lt: lt.score)

            sp = best_LType.steinerpoint
            self.unionfindsets[id].AddPoint(sp[0], sp[1])

            self.unionfindsets[id].union(best_LType.point1, sp)
            self.unionfindsets[id].union(sp, best_LType.point2)

            self.edgesets[id].addEdge(best_LType.edge1)
            self.edgesets[id].addEdge(best_LType.edge2)
            # print(f'self.edgesets[id].edges:{self.edgesets[id].edges}')

            updated_components = {}
            for point in self.unionfindsets[id].parent.keys():
                root = self.unionfindsets[id].find(point)
                if root not in updated_components:
                    updated_components[root] = set()
                updated_components[root].add(point)

            self.unionfindsets[id].components = updated_components

            self.unionfindsets[id].pointprototype.add((sp[0], sp[1]))
            self.LTypes[id].Ltypes.append(best_LType)

    def RouteGenerate(self, id: int) -> None:
        # print(f'conditions:{self.conditions[id]}')
        point_prototype_points = list(self.conditions[id])
        if len(point_prototype_points) == 1:
            # print(f'RouteGenerate {id}: Only one point to connect.')
            # print(f'point_prototype_points:{point_prototype_points}')
            self.selected_edges[id][0][point_prototype_points[0][0]][point_prototype_points[0][1]] = 1
            return
        if len(point_prototype_points) == 0:
            # print(f'RouteGenerate {id}: no point to connect.')
            # print(f'point_prototype_points:{point_prototype_points}')
            # self.selected_edges[id][0][point_prototype_points[0][0]][point_prototype_points[0][1]] = 1
            return

        heap = []
        if self.mode == 'of':
            for key, of in self.edgesets[id].oftoegdes.items():
                heap.append((of, key))  # 使用负值以实现最大堆效果
        elif self.mode == 'wl':
            for key, wl in self.edgesets[id].wltoedges.items():
                heap.append((wl, key))
        elif self.mode == 'avg_of':
            for key, avgof in self.edgesets[id].avgoftoedges.items():
                # print(f'key:{key},avgof:{avgof}')
                heap.append((avgof, key))
        elif self.mode == 'mix':
            for key, score in self.edgesets[id].scorestoedges.items():
                # print(f'key:{key},avgof:{score}')
                heap.append((score, key))
        else:
            raise ValueError(f"Unknown mode: {self.mode}")
        heapq.heapify(heap)
        # print(f'heap:{heap}')
        selected_edges_array = -1 * np.ones_like(self.PinMaps[id].cpu().numpy(), dtype=np.int32)

        def all_connected(union_find: union_findSet, points: List[tuple]) -> bool:
            roots = set()
            for point in points:
                roots.add(union_find.find_copy(point))
                if len(roots) > 1:
                    # print(f'roots:{roots}')
                    return False
            return True
        self.final_edges = set()
        while not all_connected(self.unionfindsets[id], point_prototype_points):
            # print('RouteGenerate')
            if not heap:
                raise ValueError(f"RouteGenerate {id}: cannot connect all nodes, not enough edges to pick.")

            metric, key = heapq.heappop(heap)
            p1_coords, p2_coords = key
            root1 = self.unionfindsets[id].find_copy(p1_coords)
            root2 = self.unionfindsets[id].find_copy(p2_coords)

            if root1 != root2:
                self.unionfindsets[id].union_copy(p1_coords, p2_coords)
                self.final_edges.add(key)

                if p1_coords[0] == p2_coords[0]:
                    x = p1_coords[0]
                    y_start, y_end = sorted([p1_coords[1], p2_coords[1]])
                    selected_edges_array[0, x, y_start:y_end + 1] = 1 
                elif p1_coords[1] == p2_coords[1]:
                    y = p1_coords[1]
                    x_start, x_end = sorted([p1_coords[0], p2_coords[0]])
                    selected_edges_array[0, x_start:x_end + 1, y] = 1  
                else:
                    raise ValueError(
                        f"RouteGenerate {id}: selected edge is not horizontal or vertical: {p1_coords} - {p2_coords}")

        self.selected_edges[id][0] = torch.tensor(selected_edges_array[0], dtype=torch.float32)

    def Parallelization(self,func_name:str)->None:
        func_map = {
            'UnionSetInitialize': self.UnionSetInitialize,
            'EdgesInitialize': self.EdgesInitialize,
            'PointUnion': self.PointUnion,
            # 'LTypesSetInitialize': self.LTypesSetInitialize,
            'SteinerPointAdd': self.SteinerPointAdd,
            'RouteGenerate': self.RouteGenerate
        }
        func = func_map.get(func_name)
        if func is None:
            raise ValueError(f"Unknown function name: {func_name}")

        start_time = time.time()

        with concurrent.futures.ProcessPoolExecutor(max_workers=self.args['batch_size']) as executor:
            pool = [executor.submit(func, i) for i in range(self.args['batch_size'])]

        for p in concurrent.futures.as_completed(pool):
            print(p.result())
        end_time = time.time()
        print(f'total cost: {end_time - start_time}')

    def Parallelization_tasks(self,id)->(torch.Tensor,set):
        self.InitialPointsInitialize(id)
        self.unionfindsets[id].initialize_union_find(self.PinMaps[id])
        self.edgesets[id].initialize_edges(self.PinMaps[id])
        # print(f'{id} edges:{self.edgesets[id].edgeset}')
        # print(f'{id} before parent:{self.unionfindsets[id].parent}')
        self.unionfindsets[id].merge_edges_from_edges_set(self.edgesets[id])
        # self.SteinerPointAdd(id)
        self.SteinerPointAdd_mix(id)
        # print(f'{id} after parent:{self.unionfindsets[id].parent}')
        # if self.mode == 'of':
        #     self.SteinerPointAdd_of(id)
        # self.SteinerPointAdd(id)
        self.RouteGenerate(id)
        return self.selected_edges[id],self.final_edges




    def getResult(self)->torch.Tensor:
        return self.selected_edges


if __name__ == '__main__':
    pass
import matplotlib.pyplot as plt
import networkx as nx
import time
import matplotlib.pyplot as plt
from matplotlib import patches
import math

class TangleManager:
    """
    Manages the generation, transformation, drawing, and creation of Adjacency Graphs
    of Tangles
    """
    Adjacency_Graphs = {}
    # Coordinate constants
    north, south, east, west = (0, 1), (0, -1), (1, 0), (-1, 0)
    NE, NW, SE, SW = (1, 1), (-1, 1), (1, -1), (-1, -1)


    def __init__(self):
        pass
    @staticmethod
    def _plus(p1, p2): return (p1[0] + p2[0], p1[1] + p2[1])
    @staticmethod
    def _minus(p1, p2): return (p1[0] - p2[0], p1[1] - p2[1])
    @staticmethod
    def _rotate(p): return (-p[1], p[0])
    @staticmethod
    def _reflect(p): return (-p[0], p[1])
   
    @staticmethod
    def add_vertex(Tangle, point, nbrs):
        new_T = {pt: set(ns) for pt, ns in Tangle.items()}
        new_T[point] = set(nbrs)
        for nbr in nbrs:
            new_T[nbr].add(point)
        return new_T
    
    @staticmethod
    def remove_vertex(Tangle, point):
            new_T = {pt: set(ns) for pt, ns in Tangle.items()}
            for nbr in new_T[point]:
                new_T[nbr].remove(point)
            new_T.pop(point)
            return new_T
    @staticmethod
    def get_canonical(tangle):
        """Returns a unique hashable tuple for any equivalent Tangle."""
        def get_variants(T):
            curr = T
            for _ in range(2):
                for _ in range(4):
                    yield curr
                    curr = {TangleManager._rotate(p): {TangleManager._rotate(n) for n in nbrs} for p, nbrs in curr.items()}
                curr = {TangleManager._reflect(p): {TangleManager._reflect(n) for n in nbrs} for p, nbrs in curr.items()}

        variants = []
        for v in get_variants(tangle):
            min_x, min_y = min(p[0] for p in v), min(p[1] for p in v)
            nodes = []
            for pt in sorted(v.keys()):
                new_pt = (pt[0] - min_x, pt[1] - min_y)
                new_nbrs = tuple(sorted((n[0] - min_x, n[1] - min_y) for n in v[pt]))
                nodes.append((new_pt, new_nbrs))
            variants.append(tuple(nodes))
        return min(variants)
    
    @staticmethod   
    def canonical_to_tangle(canonical_label):
        """
        Converts a canonical label form of a tangle to the dictionary form of the tangle 
        """
        tangle_dict = {}
        # Iterate through each node entry in the canonical tuple
        for point, neighbors in canonical_label:
            # Convert neighbor tuples back into sets for O(1) lookups
            tangle_dict[point] = set(neighbors)
        return tangle_dict
    
    @staticmethod
    def generate_adj_graph(start_tangle, max_depth):
        """Builds the adjacency graph up to a specified depth if not already constructed"""
        if TangleManager.get_canonical(start_tangle) in TangleManager.Adjacency_Graphs:
            return TangleManager.Adjacency_Graphs[TangleManager.get_canonical(start_tangle)]
        start_time = time.time()
        adj_graph = nx.Graph()
        canonical_map = {}
        start_can = TangleManager.get_canonical(start_tangle)
        
        adj_graph.add_node(start_can)
        canonical_map[start_can] = start_tangle
        
        current_level = {start_can}
        visited = {start_can}

        for d in range(max_depth):
            print(f"Depth {d}: Found {len(visited)} unique tangles... ({time.time()-start_time:.2f}s)")
            next_level = set()
            for can_label in current_level:
                curr_tangle = canonical_map[can_label]
                
                moves = TangleManager._get_all_moves(curr_tangle)
                
                for move in moves:
                    move_can = TangleManager.get_canonical(move)
                    if move_can not in visited:
                        visited.add(move_can)
                        canonical_map[move_can] = move
                        next_level.add(move_can)
                    adj_graph.add_edge(can_label, move_can)
            
            if not next_level: break
            current_level = next_level
        TangleManager.Adjacency_Graphs[TangleManager.get_canonical(start_tangle)] = (adj_graph, canonical_map)
        return adj_graph, canonical_map

    @staticmethod
    def _get_all_moves(tangle):
        """Internal Function that returns all possible tangles accessible from one move on the original tangle"""

        def Reflection_Moves(Tangle):
            moves = []
            for point, neighbors in Tangle.items():
                if len(neighbors) == 1:
                    neighbor = next(iter(neighbors))
                    direction = TangleManager._minus(point, neighbor)
                    # Try switching axes (reflection-like move on a grid)
                    possible_dirs = [TangleManager.east, TangleManager.west] if direction in [TangleManager.north, TangleManager.south] else [TangleManager.north, TangleManager.south]
                    for d in possible_dirs:
                        target = TangleManager._plus(neighbor, d)
                        if target not in Tangle:
                            moves.append(TangleManager.remove_vertex(TangleManager.add_vertex(Tangle, target, [neighbor]), point))
            return moves

        def Rotation_Moves(Tangle):
            moves = []
            offsets = [(TangleManager.north, TangleManager.east, TangleManager.NE), (TangleManager.north, TangleManager.west, TangleManager.NW), (TangleManager.south, TangleManager.east, TangleManager.SE), (TangleManager.south, TangleManager.west, TangleManager.SW)]
            for point, neighbors in Tangle.items():
                for d1, d2, diag in offsets:
                    p1, p2 = TangleManager._plus(point, d1), TangleManager._plus(point, d2)
                    if p1 in neighbors and p2 in neighbors:
                        p_diag = TangleManager._plus(point, diag)
                        if p_diag not in Tangle:
                            moves.append(TangleManager.add_vertex(Tangle, p_diag, [p1, p2]))
                        elif len(Tangle[p_diag]) == 2 and p1 in Tangle[p_diag] and p2 in Tangle[p_diag]:
                            moves.append(TangleManager.remove_vertex(Tangle, p_diag))
            return moves
        return Rotation_Moves(tangle) + Reflection_Moves(tangle)
    
    @staticmethod
    def draw_dual_graph(Tangle, color="purple"):
        "Uses Matplotlib to draw skeleton graph of Tangle"
        ax = plt.gca()
        plt.axis("off")
        ax.set_aspect('equal')
        for point in Tangle:
            x_val = point[0]
            y_val = point[1]
            ax.add_patch(plt.Circle(point,.05,color=color,fill = True))
            for nbr in Tangle[point]:
                x0 = nbr[0]
                y0 = nbr[1]
                plt.plot([x_val,x0], [y_val, y0], color = color)
    @staticmethod
    def draw_tangle(tangle, color = "purple", thickness = "normalized"):
        "Uses Matplotlib to draw skeleton graph of Tangle"
        if thickness == "normalized":
            thickness = 30/(max(max([p[1] for p in tangle]), max([p[0] for p in tangle]))+1)
        epsilon = 1
        radius = 1/math.sqrt(2)
        ax = plt.gca()
        plt.axis("off")
        ax.set_aspect('equal')
        plt.xlim([min([p[0] for p in tangle]) - .9*radius,max([p[0] for p in tangle]) + .9*radius ])
        plt.ylim([min([p[1] for p in tangle]) - .9*radius,max([p[1] for p in tangle]) + .9*radius ])
        for point in tangle:
            x_val = point[0]
            y_val = point[1]
            for dir, angles in [(TangleManager.north, (45,135)), (TangleManager.south, (-135,-45)), (TangleManager.east, (-45,45)), (TangleManager.west, (135,-135))]:
                if TangleManager._plus(point, dir) not in tangle[point]:
                    ax.add_patch(patches.Arc(point, radius,radius, theta1 = angles[0]-epsilon, theta2 = angles[1]+epsilon, color = color, linewidth = thickness))
                else:
                    rot_dir = (-dir[1], dir[0])
                    up = TangleManager._plus(point, rot_dir)
                    left = TangleManager._plus(point, dir)
                    diag = TangleManager._plus(up, dir)
                    if not( up in tangle[point] and diag in tangle[up] and diag in tangle[left]):
                        out_point = (point[0] + 1/2*dir[0] + 1/2*rot_dir[0], point[1] + 1/2*dir[1] + 1/2*rot_dir[1])
                        ax.add_patch(patches.Arc(out_point, radius,radius,angle = -90, theta1 = angles[0]-epsilon, theta2 = angles[1]+epsilon, color = color, linewidth = thickness))



    @staticmethod
    def shortest_path( Tangle1, Tangle2, Adjacency_Graph):
        "Returns shortest path from Tangle1 to Tangle2 in Adjacency_graph"
        path = nx.shortest_path(Adjacency_Graph, TangleManager.get_canonical(Tangle1), TangleManager.get_canonical(Tangle2))
        return [TangleManager.canonical_to_tangle(can_tang) for can_tang in path]

    @staticmethod
    def line(n):
        "creates a tangle that is a horizontal line n long"
        def add_vertex(Tangle, point, nbrs):
            new_T = {pt: set(ns) for pt, ns in Tangle.items()}
            new_T[point] = set(nbrs)
            for nbr in nbrs:
                new_T[nbr].add(point)
            return new_T
        T = {(0, 0): set()}
        for i in range(n - 1):
            T = add_vertex(T, (0, i+1), [(0, i)])
        return T
    @staticmethod
    def L_shape(n,m):
        "creates a tangle that is an L-shape that is n wide and m tall, the class of this tangle is m + n - 1"
        def add_vertex(Tangle, point, nbrs):
            new_T = {pt: set(ns) for pt, ns in Tangle.items()}
            new_T[point] = set(nbrs)
            for nbr in nbrs:
                new_T[nbr].add(point)
            return new_T
        T = {(0,0):set()}
        for i in range(n - 1):
            print(i)
            T = add_vertex(T, (i + 1, 0), [(i, 0)])
        for j in range(m - 1):
            T = add_vertex(T, (0,j + 1), [(0,j)])
        return T
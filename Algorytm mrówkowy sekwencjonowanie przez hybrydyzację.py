import random
import numpy as np
from typing import Dict, List, Tuple
import Levenshtein
import time

seq_length = 700
my_seq = ''.join(np.random.choice(('C', 'G', 'T', 'A'), seq_length))
print(my_seq)
k = 7
spektrum_idealne = seq_length - k + 1


def split(my_seq, k):
    for j in range(0, len(my_seq), 1):
        for i in range(0, len(my_seq), k):
            yield my_seq[i + j:i + j + k]


BP = 0.03*seq_length
BN = 0.03*seq_length
wr = my_seq[:k]
lista = (list(split(my_seq, k)))
lista.sort()

b = len(lista)
for a in range(b - 1, -1, -1):
    if len(lista[a]) != k:
        lista.remove(lista[a])

powtorzenia = []
unique_list = list(set(lista))
unique_list.sort()
c = len(lista) - len(unique_list)
BN = BN - c
if BN > 0:
    for c in range(1, int(BN)):
        u = random.choice(unique_list)
        powtorzenia.append(u)
        unique_list.remove(u)

print(len(unique_list))
for i in range(0, int(BP)):
    new = ''.join(np.random.choice(('C', 'G', 'T', 'A'), k))
    licz = 0
    for a in range(0, len(unique_list)):
        if unique_list[a] == new:
            licz = licz + 1
    if licz == 0:
        unique_list.append(new)
        # print(new)

unique_list.sort()
Fragment = Tuple[str, int]
Graph = Dict[str, List[Fragment]]
fragments = unique_list
coverage = []
for x in range(1, len(fragments)):
    coverage.append(40)


def generate_graph(fragments: List[str], coverage: List[int]) -> Graph:
    graph = {}
    for fragment, cov in zip(fragments, coverage):
        neighbors = []
        for neighbor, neighbor_cov in zip(fragments, coverage):
            for i in range(1, k - 1):
                if fragment[i:] == neighbor[:-i]:
                    if neighbor not in neighbors:
                        overlap = k - i
                        neighbors.append((neighbor, overlap))

        graph[fragment] = neighbors
    return graph


fragments = unique_list
graph = generate_graph(fragments, coverage)

'''
#generator losowy

def generate_path(graph: Dict[str, List[Tuple[str, int]]], start: str, length: int) -> List[str]:
    path = [start]
    current_node = start
    for i in range(length - 1):
        edges = graph[current_node]
        next_node = max(edges, key=lambda x: x[1])[0]
        path.append(next_node)
        current_node = next_node
    return path




def generate_seq(graph: Dict[str, List[Tuple[str, int]]], start: str, length: int) -> List[str]:
    path = generate_path(graph, start, length)
    sequence = path[0]
    for i in range(1, len(path)):
        if len(sequence) < seq_length:
            sequence += path[i][-seq_length + len(sequence):]
    return sequence

'''

def levenshtein_distance(seq1: str, seq2: str) -> int:
    return Levenshtein.distance(seq1, seq2)

node_indices = {}
index = 0
for node, _ in graph.items():
    node_indices[node] = index
    index += 1


num_nodes = len(node_indices)
pheromone_matrix = np.ones((num_nodes, num_nodes))


def choose_next_node(current_node: str, graph: Dict[str, List[Tuple[str, int]]], pheromone_matrix: np.ndarray,
                     node_indices: dict, alpha: float, beta: float, end_node: str) -> str:
    current_node_index = node_indices[current_node]
    prob = []
    for neighbor, overlap in graph[current_node]:
        neighbor_index = node_indices[neighbor]
        weight = (pheromone_matrix[current_node_index][neighbor_index]) * alpha * overlap * beta * (1 + overlap) * (1 + levenshtein_distance(neighbor, end_node))
        prob.append(weight)
    if sum(prob) == 0:
        raise ValueError("No available nodes to choose from")
    next_node_index =  prob.index(random.choice(prob))
    return graph[current_node][next_node_index]


def generate_apath(graph: Dict[str, List[Tuple[str, int]]], start: str, end: str, alpha: float, beta: float,
                   pheromone_matrix: np.ndarray):
    current_node = start
    path = [start]
    while len(path) < len(graph):
        next_node = choose_next_node(current_node, graph, pheromone_matrix, node_indices, alpha, beta, end)
        path.append(next_node[0])
        current_node = next_node[0]
    return path


def solve_dna_sequencing(graph: Dict[str, List[Tuple[str, int]]], start: str, end: str, num_ants: int,
                         max_iterations: int, alpha: float, beta: float, decay_rate: float, trail_value: float, pheromone_matrix: np.ndarray) -> List[
    str]:
    best_path = None
    best_score = 0
    current_path=[]
    overlap = 0
    for i in range(max_iterations):
        current_score=0
        for j in range(num_ants):
            try:
                current_path = generate_apath(graph, start, end, alpha, beta, pheromone_matrix)
                for l in range(1, len(current_path)):
                    overlap = 0
                    for m in range(len(current_path[0])-1, 1, -1):
                        if current_path[l][:-m] != current_path[l-1][m:]:
                            overlap += len(current_path[0]) - m
                            break
                current_score += overlap
            except ValueError as e:
                if str(e) == "No available nodes to choose from":
                    max_iterations += 1
                    continue
            if current_score > best_score:
                best_path = current_path
                best_score = current_score
        update_pheromone_matrix(best_path, trail_value, decay_rate, pheromone_matrix)
    return best_path


def update_pheromone_matrix(path: List[str], trail_value: float,
                            decay_rate: float, pheromone_matrix: np.ndarray):
    for i in range(len(path) - 1):
        pheromone_matrix[i][i+1] = (1 - decay_rate) * pheromone_matrix[i][i+1] + trail_value
        pheromone_matrix[i][i+1] = pheromone_matrix[i][i+1]


def generate_aseq(graph: Dict[str, List[Tuple[str, int]]], start: str, end: str, num_ants: int, max_iterations: int,
                  alpha: float, beta: float, trail_value: float, decay_rate: float, pheromone_matrix: np.ndarray, true_sequence: str) -> str:
    best_path = solve_dna_sequencing(graph, start, end, num_ants, max_iterations, alpha, beta, decay_rate, trail_value, pheromone_matrix)
    print(best_path)
    sequence = best_path[0]

    for i in range(1, len(best_path)):
        overlap = 0
        for j in range(min(len(best_path[i]), len(sequence))):
            if best_path[i][j] == sequence[- j]:
                overlap += 1
            else:
                break
        sequence += (best_path[i][overlap:])
        if len(sequence) >= len(true_sequence):
            break
    return sequence


def test_dna_sequencing(graph: Dict[str, List[Tuple[str, int]]], start: str, end: str, true_sequence: str,
                        num_ants: int, max_iterations: int, alpha: float, beta: float, trail_value: float,
                        decay_rate: float, pheromone_matrix: np.ndarray):
    start_time = time.time()
    generated_sequence = generate_aseq(graph, start, end, num_ants, max_iterations, alpha, beta, trail_value, decay_rate,  pheromone_matrix, true_sequence)
    print(generated_sequence)
    end_time = time.time()
    print("Time taken: ", end_time - start_time)
    print("Levenshtein distance: ", levenshtein_distance(generated_sequence, true_sequence))


test_dna_sequencing(graph, my_seq[:k], fragments[-1], my_seq, 50, 10, 0.01, 0.01, 0.1, 0.1, pheromone_matrix)

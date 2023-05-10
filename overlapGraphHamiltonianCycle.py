#python3

"""
Task. Given a list of error-free reads, perform the task of Genome Assembly and return the circular genome
from which they came.

Dataset. Each of 1 618 lines of the input contains a single read, that is, a string over {A, C, G, T}. The
reads are given to you in alphabetical order because their true order is hidden from you. Each read
is 100 nucleotides long and contains no sequencing errors. Note that you are not given the 100-mer
composition of the genome, i.e., some 100-mers may be missing.

Output. Output the assembled genome on a single line.
"""

from itertools import permutations
import sys
sys.setrecursionlimit(1700)

MIN_KMER = 12
SAMPLE_READS = ["AAC", "ACG", "GAA", "GTT", "TCG"]
NUM_READS = 1618

class GenomeAssembler():

    def __init__(self, reads):
        self.reads = set(reads)
        self.overlapMap, self.overlapEdges = self.overlap_map()
        self.explored = set()


    def overlap(self, x: str, y: str, min_overlap=MIN_KMER):
        """
        Return lenght of the longest suffix of `x` that matches a prefix of `y`
        if there is a match of lenght `min_overlap`. If no such overlap exists return 0.
        """
        start = 0
        while True:
            start = x.find(y[:min_overlap], start)
            if start == -1:
                return 0
            if y.startswith(x[start:]):
                return len(x) - start    
            start += 1


    def overlap_map(self):
        """
        Return dictionary of keys `suffix read` and values list of tuples  (`prefix read`, `overlap lenght`) sorted by non increasing `overlap lenght`
        """
        overlapEdges = {}
        overlaps = {read:{} for read in self.reads}
        for x,y in permutations(self.reads, 2):
            olen = self.overlap(x,y)
            if olen:
                overlapEdges[(x,y)] = olen
                overlaps[x][y] = olen
        for read in self.reads:
            overlaps[read] = sorted(overlaps[read].items(), key=lambda d: d[1], reverse=True)

        return overlaps, overlapEdges

    
    def choose_largest_overlaping_reads(self):
        """
        Return a pair of reads `x`, `y` of max overlap `largestOverlap` from a set of reads 
        """
        x, y = None, None
        largestOverlap = MIN_KMER
        for a,b in permutations(self.reads, 2):
            olen = self.overlap(a,b, largestOverlap)
            if olen > largestOverlap:
                x, y = a, b
                largestOverlap = olen
        return x, y, largestOverlap


    def is_hamiltonian_cycle(self, path) -> bool:
        """
        Return `True` if a list of reads `path` is a valid Hamiltonian cycle, otherwise return `False`
        """
        # if path includes all vertices once
        if set(path) == self.reads and len(path) == len(self.reads):

            # if all vertices are connected to their next vertex in the path by an edge
            for i in range(len(path)-1):
                edge = (path[i],path[i+1])
                if edge not in self.overlapEdges:
                    return False
                
            # if the first and the last vertex are connected by an edge
            if (path[-1], path[0]) in self.overlapEdges:
                return True
        return False

    
    def rankedOverlaps(self, vertex):
        """
        Return neighboring vertices sorted by non increasing overlap lenght 
        """
        return [x for x,_ in self.overlapMap[vertex]]

    
    def recursivehamiltonian_path(self, path:list):
        """
        Return a Hamiltonian path
        """
    
        # return a valid Hamiltonian path when possible 
        if self.is_hamiltonian_cycle(path):
            return path

        # greedily explore neighbors of last vertex
        vertex = path[-1]
        for v in  self.rankedOverlaps(vertex):

            nextPath = path.copy()

            if v not in self.explored:

                nextPath.append(v)
                self.explored.add(v)
                result = self.recursivehamiltonian_path(nextPath)
                if result is not None:
                    return result

                # backtrack if bfs found no hamiltonian cycle down this branch
                nextPath = nextPath[:-1]
                self.explored.remove(v)
        return

    
    def overhang(self, string:str) -> int:
        """
        Return max lenght `n` of any substring of `string` that `string[:n] == string[-n:]` 
        """
        n = 0
        for i in range(1,len(string)):
            left = string[:i]
            right = string[-i:]
            if left == right:
                n = i
        return n


    def gluePath(self, path:list) -> str:
        """
        Return string circular genome from a hamiltonian cycle
        """
        a, b = path[0], path[1]
        olen = self.overlapEdges[(a,b)]
        genome = a+b[olen:]
        
        for i in range(2,len(path)):
            nextRead = path[i]
            olen = self.overlap(genome, nextRead)
            genome += nextRead[olen:]

        overhang = self.overhang(genome)
        if overhang:
            return genome[:-overhang]
        else:
            return genome


reads = [input() for _ in range(NUM_READS)]
g = GenomeAssembler(reads)
first, second,_ = g.choose_largest_overlaping_reads()
g.explored.add(first)
g.explored.add(second)
path = g.recursivehamiltonian_path([first, second])
print(g.gluePath(path))
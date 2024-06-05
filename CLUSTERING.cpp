#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <cstdlib>

// Structure to represent a point in R^n
struct Point {
    std::vector<double> coordinates;
    int cluster;
};

// Structure to represent an edge between two points with weight
struct Edge {
    int src, dest;
    double weight;
};

// Structure to represent a subset for union-find
struct Subset {
    int parent;
    int rank;
};

// Function to calculate Euclidean distance between two points
double CalculateDistance(const Point& p1, const Point& p2, int dimensionality) {
    double distance = 0;
    for (int i = 0; i < dimensionality; i++) {
        distance += std::pow(p1.coordinates[i] - p2.coordinates[i], 2);
    }
    return std::sqrt(distance);
}

// Function to find the subset of an element i
int Find(std::vector<Subset>& subsets, int i) {
    if (subsets[i].parent != i)
        subsets[i].parent = Find(subsets, subsets[i].parent);
    return subsets[i].parent;
}

// Function to perform union of two sets
void UnionSets(std::vector<Subset>& subsets, int x, int y) {
    int xroot = Find(subsets, x);
    int yroot = Find(subsets, y);

    if (subsets[xroot].rank < subsets[yroot].rank)
        subsets[xroot].parent = yroot;
    else if (subsets[xroot].rank > subsets[yroot].rank)
        subsets[yroot].parent = xroot;
    else {
        subsets[yroot].parent = xroot;
        subsets[xroot].rank++;
    }
}

// Function to swap two edges
void exchange(Edge* a, Edge* b) {
    Edge temp = *a;
    *a = *b;
    *b = temp;
}

// Partition function for randomized quicksort
int partition(std::vector<Edge>& edges, int low, int high) {
    double pivot = edges[high].weight;
    int i = low - 1;
    for (int j = low; j < high; j++) {
        if (edges[j].weight <= pivot) {
            i++;
            exchange(&edges[i], &edges[j]);
        }
    }
    exchange(&edges[i + 1], &edges[high]);
    return i + 1;
}

// Random partition function
int rand_part(std::vector<Edge>& edges, int low, int high) {
    int random = low + std::rand() % (high - low + 1);
    exchange(&edges[random], &edges[high]);
    return partition(edges, low, high);
}

// Randomized quicksort function
void rand_qsort(std::vector<Edge>& edges, int low, int high) {
    if (low < high) {
        int q = rand_part(edges, low, high);
        rand_qsort(edges, low, q - 1);
        rand_qsort(edges, q + 1, high);
    }
}

// Function to find the minimum spanning tree using Kruskal's algorithm
void KruskalMST(std::vector<Point>& points, int NumPoints, int dimensionality, int k) {
    int NumEdges = NumPoints * (NumPoints - 1) / 2;
    std::vector<Edge> edges(NumEdges);

    // Calculate distances and populate edges
    int edgeIndex = 0;
    for (int i = 0; i < NumPoints; i++) {
        for (int j = i + 1; j < NumPoints; j++) {
            edges[edgeIndex] = {i, j, CalculateDistance(points[i], points[j], dimensionality)};
            edgeIndex++;
        }
    }

    // Sort edges by weight using randomized quicksort
    rand_qsort(edges, 0, NumEdges - 1);

    // Allocate memory for subsets
    std::vector<Subset> subsets(NumPoints);

    // Initialize subsets
    for (int i = 0; i < NumPoints; i++) {
        subsets[i].parent = i;
        subsets[i].rank = 0;
    }

    int edgeCount = 0;
    int index = 0;
    while (edgeCount < NumPoints - 1 - (k - 1)) { // For the k-clustering
        Edge nextEdge = edges[index++];
        int x = Find(subsets, nextEdge.src);
        int y = Find(subsets, nextEdge.dest);

        if (x != y) {
            UnionSets(subsets, x, y);
            edgeCount++;
        }
    }

    // Assign clusters
    std::vector<int> ClusRoot(NumPoints, 0);
    int x = 0; // Cluster-count
    for (int i = 0; i < NumPoints; i++) {
        if (ClusRoot[Find(subsets, i)] == 0) {
            ClusRoot[Find(subsets, i)] = ++x;
        }
        points[i].cluster = ClusRoot[Find(subsets, i)];
    }
}

int main(int argc, char* argv[]) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <input_file> <output_file>\n";
        return 1;
    }

    std::ifstream inputFile(argv[1]);
    if (!inputFile) {
        std::cerr << "Error: Unable to open input file.\n";
        return 1;
    }

    int NumPoints, dimensionality, k;
    inputFile >> NumPoints >> dimensionality >> k;

    std::vector<Point> points(NumPoints);
    for (int i = 0; i < NumPoints; i++) {
        points[i].coordinates.resize(dimensionality);
        for (int j = 0; j < dimensionality; j++) {
            inputFile >> points[i].coordinates[j];
        }
    }
    inputFile.close();

    // Find clusters using MST
    KruskalMST(points, NumPoints, dimensionality, k);

    std::ofstream outputFile(argv[2]);
    if (!outputFile) {
        std::cerr << "Error: Unable to open output file.\n";
        return 1;
    }

    // Write points with cluster names to output file
    for (const auto& point : points) {
        for (const auto& coord : point.coordinates) {
            outputFile << coord << " ";
        }
        outputFile << point.cluster << "\n";
    }
    outputFile.close();

    return 0;
}

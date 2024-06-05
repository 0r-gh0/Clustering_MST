#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <cstdlib>
#include <ctime>

using namespace std;

// A structure to represent a weighted edge in graph
struct Edge {
    int src, dest;
    double weight;
};

// A structure to represent a subset for union-find
struct Subset {
    int parent;
    int rank;
};

// Function to find the Euclidean distance between two points
double calculateDistance(const vector<double>& point1, const vector<double>& point2) {
    double sum = 0;
    for (size_t i = 0; i < point1.size(); ++i) {
        sum += pow(point1[i] - point2[i], 2);
    }
    return sqrt(sum);
}

// A function to find set of an element i (uses path compression technique)
int find(vector<Subset>& subsets, int i) {
    if (subsets[i].parent != i) {
        subsets[i].parent = find(subsets, subsets[i].parent);
    }
    return subsets[i].parent;
}

// A function that does union of two sets of x and y (uses union by rank)
void Union(vector<Subset>& subsets, int x, int y) {
    int rootX = find(subsets, x);
    int rootY = find(subsets, y);

    if (subsets[rootX].rank < subsets[rootY].rank) {
        subsets[rootX].parent = rootY;
    } else if (subsets[rootX].rank > subsets[rootY].rank) {
        subsets[rootY].parent = rootX;
    } else {
        subsets[rootY].parent = rootX;
        subsets[rootX].rank++;
    }
}

// Function to sort edges based on their weights using randomized quicksort
int partition(vector<Edge>& edges, int low, int high) {
    double pivot = edges[high].weight;
    int i = (low - 1);

    for (int j = low; j <= high - 1; j++) {
        if (edges[j].weight <= pivot) {
            i++;
            swap(edges[i], edges[j]);
        }
    }
    swap(edges[i + 1], edges[high]);
    return (i + 1);
}

void randomizedQuicksort(vector<Edge>& edges, int low, int high) {
    if (low < high) {
        int pivot = low + rand() % (high - low + 1);
        swap(edges[pivot], edges[high]);
        int pi = partition(edges, low, high);
        randomizedQuicksort(edges, low, pi - 1);
        randomizedQuicksort(edges, pi + 1, high);
    }
}

// Function to construct MST using Kruskal's algorithm
vector<Edge> KruskalMST(vector<Edge>& edges, int V) {
    vector<Edge> result; // This will store the resultant MST
    vector<Subset> subsets(V);

    for (int v = 0; v < V; ++v) {
        subsets[v].parent = v;
        subsets[v].rank = 0;
    }

    randomizedQuicksort(edges, 0, edges.size() - 1);

    int e = 0; // Number of edges in the result
    int i = 0; // Initial index of sorted edges

    while (e < V - 1 && i < edges.size()) {
        Edge nextEdge = edges[i++];

        int x = find(subsets, nextEdge.src);
        int y = find(subsets, nextEdge.dest);

        if (x != y) {
            result.push_back(nextEdge);
            Union(subsets, x, y);
            e++;
        }
    }

    return result;
}

void assignClusters(vector<int>& clusters, vector<Subset>& subsets, int numPoints) {
    int clusterId = 0;
    for (int i = 0; i < numPoints; ++i) {
        int root = find(subsets, i);
        if (clusters[root] == -1) {
            clusters[root] = clusterId++;
        }
        clusters[i] = clusters[root];
    }
}

int main(int argc, char** argv) {
    if (argc != 3) {
        cerr << "Usage: " << argv[0] << " <input_file> <output_file>" << endl;
        return 1;
    }

    ifstream inputFile(argv[1]);
    if (!inputFile) {
        cerr << "Error opening input file." << endl;
        return 1;
    }

    int numPoints, dimensionality, k;
    inputFile >> numPoints >> dimensionality >> k;

    vector<vector<double>> points(numPoints, vector<double>(dimensionality));
    for (int i = 0; i < numPoints; ++i) {
        for (int j = 0; j < dimensionality; ++j) {
            inputFile >> points[i][j];
        }
    }
    inputFile.close();

    vector<Edge> edges;
    for (int i = 0; i < numPoints; ++i) {
        for (int j = i + 1; j < numPoints; ++j) {
            double dist = calculateDistance(points[i], points[j]);
            edges.push_back({i, j, dist});
        }
    }

    vector<Edge> mst = KruskalMST(edges, numPoints);

    // Remove the k-1 largest edges from the MST to form k clusters
    sort(mst.begin(), mst.end(), [](Edge a, Edge b) {
        return a.weight > b.weight;
    });

    for (int i = 0; i < k - 1; ++i) {
        mst.pop_back();
    }

    vector<int> clusters(numPoints, -1);
    vector<Subset> subsets(numPoints);
    for (int v = 0; v < numPoints; ++v) {
        subsets[v].parent = v;
        subsets[v].rank = 0;
    }

    for (const Edge& edge : mst) {
        Union(subsets, edge.src, edge.dest);
    }

    assignClusters(clusters, subsets, numPoints);

    ofstream outputFile(argv[2]);
    if (!outputFile) {
        cerr << "Error opening output file." << endl;
        return 1;
    }

    for (int i = 0; i < numPoints; ++i) {
        for (int j = 0; j < dimensionality; ++j) {
            outputFile << points[i][j] << " ";
        }
        outputFile << clusters[i] + 1 << endl; // Cluster labels are 1-based
    }

    outputFile.close();
    return 0;
}
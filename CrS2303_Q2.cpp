#include <iostream>      // Include the standard input-output stream library
#include <vector>        // Include the vector library for using dynamic arrays
#include <cmath>         // Include the cmath library for mathematical functions
#include <cstdlib>       // Include the cstdlib library for random number generation
#include <cstdio>        // Include the cstdio library for standard input and output
#include <fstream>       // Include the file stream library for file handling

using namespace std;     // Use the standard namespace

// Structure to represent a Tree node for Union-Find
struct Tree {
    int rank;     // Rank of the tree node for Union-Find
    int parent;   // Parent of the tree node for Union-Find
    int st;       // Cluster identifier
};

// Structure to represent an Edge
struct Edge {
    int in;      // Start node of the edge
    int out;     // End node of the edge
    float wght;  // Weight of the edge
};

// Global vectors to store tree nodes and edges
vector<Tree> t;
vector<Edge> s;

// Function to partition the edges for quicksort
int partition(int i, int j) {
    int p, q, r;          // Declare indices for partitioning
    Edge temp;            // Temporary variable for swapping
    p = (rand() % (j - i + 1)) + i;  // Randomly select a pivot
    temp = s[p];          // Swap pivot with the end element
    s[p] = s[j];
    s[j] = temp;

    q = i - 1;            // Initialize the partition index
    for (r = i; r < j; r++) {   // Partition the array
        if (s[r].wght < s[j].wght) {   // Compare edge weights
            q = q + 1;
            temp = s[q];
            s[q] = s[r];
            s[r] = temp;
        }
    }
    temp = s[q + 1];      // Swap the pivot to its correct position
    s[q + 1] = s[j];
    s[j] = temp;
    return (q + 1);       // Return the partition index
}

// Function to perform quicksort on the edges
void quicksort(int i, int j) {
    int q;                // Declare the partition index
    if (i < j) {          // Check if the segment has more than one element
        q = partition(i, j);  // Partition the segment
        quicksort(i, q - 1);  // Recursively sort the left segment
        quicksort(q + 1, j);  // Recursively sort the right segment
    }
}

// Function to find the root of a tree node with path compression
int find(int i) {
    while (i != t[i].parent)  // Traverse up to the root
        i = t[i].parent;
    return (i);              // Return the root
}

// Function to perform union of two sets
void Union(int i, int j) {
    int p, q;                // Declare roots of the sets
    p = find(i);             // Find root of the first set
    q = find(j);             // Find root of the second set
    if (t[p].rank > t[q].rank) {  // Attach the smaller tree under the larger tree
        t[q].parent = p;
    } else {
        t[p].parent = q;
        if (t[p].rank == t[q].rank)  // If ranks are equal, increment the rank
            t[q].rank = t[q].rank + 1;
    }
}

// Main function
int main() {
    int i, j, l, n, k, dim, e, T, p; // Declare variables
    vector<int> set;                 // Vector to store unique clusters
    vector<vector<float>> m, d;      // Vectors to store input data and distance matrix
    float cmpt;                      // Variable to store computed distance

    string inputFileName;            // Variable for input file name
    string outputFileName;           // Variable for output file name

    // Prompt the user to enter the input file name
    cout << "Enter the name of the Input file : ";
    getline(cin, inputFileName);

    ifstream fin(inputFileName);    // Open the input file

    fin >> n >> dim >> k;           // Read the number of points, dimensions, and clusters

    m.resize(n, vector<float>(dim));  // Resize matrix to hold the input data
    for (i = 0; i < n; i++) {          // Read the input data
        for (j = 0; j < dim; j++)
            fin >> m[i][j];
    }

    fin.close();  // Close the input file

    t.resize(n);  // Resize the tree vector
    for (i = 0; i < n; i++) {  // Initialize the tree nodes
        t[i].rank = 0;
        t[i].parent = i;
        t[i].st = 0;
    }

    e = n * (n - 1) / 2;  // Calculate the number of edges
    s.resize(e);  // Resize the edges vector

    d.resize(n, vector<float>(n));  // Resize the distance matrix
    for (i = 0; i < n; i++) {       // Initialize the distance matrix
        d[i][i] = 0;
        for (j = i + 1; j < n; j++) {
            cmpt = 0.0;
            for (l = 0; l < dim; l++)
                cmpt += pow((m[i][l] - m[j][l]), 2);
            d[i][j] = sqrt(cmpt);  // Compute the distance
            d[j][i] = d[i][j];     // Symmetric distance
        }
    }

    l = 0;  // Initialize edge index
    for (i = 0; i < n; i++) {       // Create edges with weights
        for (j = i + 1; j < n; j++) {
            s[l].in = i;
            s[l].out = j;
            s[l].wght = d[i][j];
            l++;
        }
    }

    quicksort(0, l - 1);  // Sort the edges using quicksort

    i = 0;  // Initialize edge index
    T = 0;  // Initialize the number of edges in the MST
    while (T < (n - k)) {  // Loop until we have n-k clusters
        if (find(s[i].in) != find(s[i].out)) {  // Check if the edge connects different components
            Union(s[i].in, s[i].out);  // Union the components
            T = T + 1;
        }
        i = i + 1;  // Move to the next edge
    }

    set.resize(k, -1);  // Initialize the set vector
    l = 0;  // Initialize set index
    for (i = 0; i < n; i++) {       // Create unique sets
        p = find(i);
        for (j = 0; j < k; j++) {
            if (set[j] == p)
                break;
        }
        if (j == k) {  // If not found, add to set
            set[l] = p;
            l++;
        }
    }

    for (i = 0; i < n; i++) {       // Assign cluster identifiers
        p = find(i);
        for (j = 0; j < k; j++) {
            if (set[j] == p)
                break;
        }
        t[i].st = j;  // Assign the cluster number
    }

    // Prompt the user to enter the output file name
    cout << "Enter the name of the Output file : ";
    getline(cin, outputFileName);

    ofstream fout(outputFileName);  // Open the output file
    for (i = 0; i < n; i++) {       // Write the points and their clusters to the output file
        for (j = 0; j < dim; j++)
            fout << m[i][j] << " ";
        fout << "--> " << (t[i].st + 1) << endl;
    }

    fout.close();  // Close the output file

    return 0;  // Return success code
}

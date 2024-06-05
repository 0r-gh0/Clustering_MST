#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <fstream>

using namespace std;

struct Tree {
    int rank;
    int parent;
    int st;
};

struct Edge {
    int in;
    int out;
    float wght;
};

vector<Tree> t;
vector<Edge> s;

int partition(int i, int j) {
    int p, q, r;
    Edge temp;
    p = (rand() % (j - i + 1)) + i;
    temp = s[p];
    s[p] = s[j];
    s[j] = temp;

    q = i - 1;
    for (r = i; r < j; r++) {
        if (s[r].wght < s[j].wght) {
            q = q + 1;
            temp = s[q];
            s[q] = s[r];
            s[r] = temp;
        }
    }
    temp = s[q + 1];
    s[q + 1] = s[j];
    s[j] = temp;
    return (q + 1);
}

void quicksort(int i, int j) {
    int q;
    if (i < j) {
        q = partition(i, j);
        quicksort(i, q - 1);
        quicksort(q + 1, j);
    }
}

int find(int i) {
    while (i != t[i].parent)
        i = t[i].parent;
    return (i);
}

void Union(int i, int j) {
    int p, q;
    p = find(i);
    q = find(j);
    if (t[p].rank > t[q].rank) {
        t[q].parent = p;
    } else {
        t[p].parent = q;
        if (t[p].rank == t[q].rank)
            t[q].rank = t[q].rank + 1;
    }
}

int main() {
    int i, j, l, n, k, dim, e, T, p;
    vector<int> set;
    vector<vector<float>> m, d;
    float cmpt;

    string inputFileName;
    string outputFileName;

    cout << "Enter the name of the Input file : ";
    getline(cin, inputFileName);

    ifstream fin(inputFileName);

    fin >> n >> dim >> k;

    m.resize(n, vector<float>(dim));
    for (i = 0; i < n; i++) {
        for (j = 0; j < dim; j++)
            fin >> m[i][j];
    }

    fin.close();

    t.resize(n);
    for (i = 0; i < n; i++) {
        t[i].rank = 0;
        t[i].parent = i;
        t[i].st = 0;
    }

    e = n * (n - 1) / 2;
    s.resize(e);

    d.resize(n, vector<float>(n));
    for (i = 0; i < n; i++) {
        d[i][i] = 0;
        for (j = i + 1; j < n; j++) {
            cmpt = 0.0;
            for (l = 0; l < dim; l++)
                cmpt += pow((m[i][l] - m[j][l]), 2);
            d[i][j] = sqrt(cmpt);
            d[j][i] = d[i][j];
        }
    }

    l = 0;
    for (i = 0; i < n; i++) {
        for (j = i + 1; j < n; j++) {
            s[l].in = i;
            s[l].out = j;
            s[l].wght = d[i][j];
            l++;
        }
    }

    quicksort(0, l - 1);

    i = 0;
    T = 0;
    while (T < (n - k)) {
        if (find(s[i].in) != find(s[i].out)) {
            Union(s[i].in, s[i].out);
            T = T + 1;
        }
        i = i + 1;
    }

    set.resize(k, -1);
    l = 0;
    for (i = 0; i < n; i++) {
        p = find(i);
        for (j = 0; j < k; j++) {
            if (set[j] == p)
                break;
        }
        if (j == k) {
            set[l] = p;
            l++;
        }
    }

    for (i = 0; i < n; i++) {
        p = find(i);
        for (j = 0; j < k; j++) {
            if (set[j] == p)
                break;
        }
        t[i].st = j;
    }

    cout << "Enter the name of the Output file : ";
    getline(cin, outputFileName);

    ofstream fout(outputFileName);
    for (i = 0; i < n; i++) {
        for (j = 0; j < dim; j++)
            fout << m[i][j] << " ";
        fout << "--> " << (t[i].st + 1) << endl;
    }

    fout.close();

    return 0;
}

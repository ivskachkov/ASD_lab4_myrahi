#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <limits>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <sstream>

using namespace std;

const int N = 200;           // Кількість вершин
const int M = 45;            // Кількість мурах
const int EliteAnts = 10;    // Елітні мурахи
const double Alpha = 3.0;    // Вага феромону
const double Beta = 2.0;     // Вага видимості
const double Rho = 0.7;      // Коефіцієнт випаровування феромону
const int MaxIterations = 1000;

vector<vector<int>> distances(N, vector<int>(N));
vector<vector<double>> pheromone(N, vector<double>(N, 1.0));
vector<vector<double>> visibility(N, vector<double>(N));     

void generateDistances() {
    srand(time(0));
    for (int i = 0; i < N; ++i) {
        for (int j = i + 1; j < N; ++j) {
            distances[i][j] = distances[j][i] = rand() % 40 + 1; 
            visibility[i][j] = visibility[j][i] = 1.0 / distances[i][j];
        }
    }
}

void saveDistances(const string& filename) {
    ofstream outFile(filename);
    if (!outFile) {
        cerr << "Помилка: не вдалося відкрити файл для запису." << endl;
        return;
    }
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            outFile << distances[i][j];
            if (j < N - 1) outFile << ",";
        }
        outFile << endl;
    }
    outFile.close();
}

// Функція для завантаження матриці відстаней із файлу
bool loadDistances(const string& filename) {
    ifstream inFile(filename);
    if (!inFile) {
        cerr << "Помилка: файл " << filename << " не знайдено. Буде створена нова матриця." << endl;
        return false;
    }

    string line;
    for (int i = 0; i < N && getline(inFile, line); ++i) {
        stringstream ss(line);
        string value;
        for (int j = 0; j < N && getline(ss, value, ','); ++j) {
            distances[i][j] = stoi(value);
            if (i != j) {
                visibility[i][j] = 1.0 / distances[i][j];
            }
        }
    }

    inFile.close();
    return true;
}

int greedyTSP(vector<int>& bestPath) {
    vector<bool> visited(N, false);
    bestPath.clear();
    int current = 0;
    visited[current] = true;
    bestPath.push_back(current);
    int totalDistance = 0;

    for (int step = 1; step < N; ++step) {
        int nearest = -1;
        int nearestDist = numeric_limits<int>::max();

        for (int i = 0; i < N; ++i) {
            if (!visited[i] && distances[current][i] < nearestDist) {
                nearest = i;
                nearestDist = distances[current][i];
            }
        }

        visited[nearest] = true;
        totalDistance += nearestDist;
        bestPath.push_back(nearest);
        current = nearest;
    }

    totalDistance += distances[current][bestPath[0]];
    return totalDistance;
}

void updatePheromones(const vector<vector<int>>& antPaths, const vector<int>& antDistances, int bestDistance, const vector<int>& bestPath) {
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            pheromone[i][j] *= (1 - Rho);
        }
    }

    for (int k = 0; k < M; ++k) {
        double deltaPheromone = (k < EliteAnts ? 2.0 : 1.0) / antDistances[k];
        for (int i = 0; i < N - 1; ++i) {
            int from = antPaths[k][i];
            int to = antPaths[k][i + 1];
            pheromone[from][to] += deltaPheromone;
            pheromone[to][from] += deltaPheromone;
        }
    }
}

int antTour(vector<int>& path, int start) {
    vector<bool> visited(N, false);
    path.clear();
    path.push_back(start);
    visited[start] = true;

    int totalDistance = 0;
    int current = start;

    for (int step = 1; step < N; ++step) {
        double sumProb = 0.0;
        vector<double> probabilities(N, 0.0);

        for (int i = 0; i < N; ++i) {
            if (!visited[i]) {
                probabilities[i] = pow(pheromone[current][i], Alpha) * pow(visibility[current][i], Beta);
                sumProb += probabilities[i];
            }
        }

        double randomValue = (double)rand() / RAND_MAX * sumProb;
        double cumulative = 0.0;

        int next = -1;
        for (int i = 0; i < N; ++i) {
            if (!visited[i]) {
                cumulative += probabilities[i];
                if (cumulative >= randomValue) {
                    next = i;
                    break;
                }
            }
        }

        visited[next] = true;
        path.push_back(next);
        totalDistance += distances[current][next];
        current = next;
    }

    totalDistance += distances[current][path[0]];
    return totalDistance;
}

void antColonyOptimization() {
    vector<int> bestPath;
    int bestDistance = greedyTSP(bestPath);

    cout << "Greedy solution distance: " << bestDistance << endl;

    for (int iter = 0; iter < MaxIterations; ++iter) {
        vector<vector<int>> antPaths(M);
        vector<int> antDistances(M);

        for (int k = 0; k < M; ++k) {
            int start = rand() % N;
            antDistances[k] = antTour(antPaths[k], start);

            if (antDistances[k] < bestDistance) {
                bestDistance = antDistances[k];
                bestPath = antPaths[k];
            }
        }

        updatePheromones(antPaths, antDistances, bestDistance, bestPath);
        cout << "Iteration " << iter + 1 << ": Best distance = " << bestDistance << endl;
    }

    cout << "Best path distance: " << bestDistance << endl;
}

int main() {
    srand(time(0));

    const string filename = "distances.csv";

    if (!loadDistances(filename)) {
        generateDistances();
        saveDistances(filename);
    }

    antColonyOptimization();
    return 0;
}

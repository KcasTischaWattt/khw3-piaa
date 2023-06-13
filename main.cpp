#include <iostream>
#include <vector>
#include <set>
#include <random>
#include <cstdlib>
#include <ctime>
#include <chrono>
#include <fstream>
#include <string>
#include <queue>

/// Тип графа
std::string gr_type;

/// Выбор типа графа
void setGrType(int type) {
    switch (type) {
        case 0:
            gr_type = "Full";
            break;
        case 1:
            gr_type = "Semi";
            break;
        default:
            gr_type = "Empty";
    }
}

/// Класс для генерации графов
class GraphGenerator {
public:

    // Генерация полного графа
    static std::vector<std::vector<std::vector<std::pair<int64_t, int64_t>>>>
    generateCompleteGraphs(int64_t start_vertices, int64_t end_vertices, int64_t step) {
        std::vector<std::vector<std::vector<std::pair<int64_t, int64_t>>>> graphs;
        for (int64_t n = start_vertices; n <= end_vertices; n += step) {
            std::vector<std::vector<std::pair<int64_t, int64_t>>> g(n);
            for (int64_t i = 0; i < n; ++i) {
                for (int64_t j = i + 1; j < n; ++j) {
                    int64_t weight = generateRandomWeight();
                    g[i].emplace_back(j, weight);
                    g[j].emplace_back(i, weight);
                }
            }
            graphs.push_back(g);
        }
        return graphs;
    }

    // Генерация связного графа с коэффициентом плотности density
    static std::vector<std::vector<std::vector<std::pair<int64_t, int64_t>>>>
    generateConnectedGraphs(int64_t start_vertices, int64_t end_vertices, int64_t step,
                            double density) {
        std::vector<std::vector<std::vector<std::pair<int64_t, int64_t>>>> graphs;
        for (int64_t n = start_vertices; n <= end_vertices; n += step) {
            std::vector<std::vector<std::pair<int64_t, int64_t>>> g(n);
            for (int64_t i = 0; i < n - 1; ++i) {
                for (int64_t j = i + 1; j < n; ++j) {
                    if (generateRandomProbability() < density) {
                        int64_t weight = generateRandomWeight();
                        g[i].emplace_back(j, weight);
                        g[j].emplace_back(i, weight);
                    }
                }
            }
            graphs.push_back(g);
        }
        return graphs;
    }

    // Генерация разреженного графа (дерева)
    static std::vector<std::vector<std::vector<std::pair<int64_t, int64_t>>>>
    generateSparseGraphs(int64_t start_vertices, int64_t end_vertices, int64_t step) {
        std::vector<std::vector<std::vector<std::pair<int64_t, int64_t>>>> graphs;
        for (int64_t n = start_vertices; n <= end_vertices; n += step) {
            std::vector<std::vector<std::pair<int64_t, int64_t>>> g(n);
            for (int64_t i = 0; i < n - 1; ++i) {
                int64_t parent = generateRandomVertex(0, i);
                int64_t weight = generateRandomWeight();
                g[i + 1].emplace_back(parent, weight);
                g[parent].emplace_back(i + 1, weight);
            }
            graphs.push_back(g);
        }
        return graphs;
    }


private:
    // Генерация случайного веса ребра
    static int64_t generateRandomWeight() {
        static std::random_device rd;
        static std::mt19937 gen(rd());
        static std::uniform_int_distribution<> dis(1, 10);
        return dis(gen);
    }

    // Генерация случайной вершины от low до high (включительно)
    static int64_t generateRandomVertex(int64_t low, int64_t high) {
        static std::random_device rd;
        static std::mt19937 gen(rd());
        static std::uniform_int_distribution<> dis(low, high);
        return dis(gen);
    }

    // Генерация случайной вероятности
    static double generateRandomProbability() {
        static std::random_device rd;
        static std::mt19937 gen(rd());
        static std::uniform_real_distribution<> dis(0.0, 1.0);
        return dis(gen);
    }
};


// дейкстра взят с лекции
void dijkstra(std::vector<std::vector<std::pair<int64_t, int64_t>>> graph, int64_t start_vertex) {
    size_t n = graph.size();
    std::vector<int64_t> d(n, INT_MAX);
    std::vector<int64_t> p(n);
    d[start_vertex] = 0;

    std::set<std::pair<int64_t, int64_t>> q;
    q.insert(std::make_pair(d[start_vertex], start_vertex));

    while (!q.empty()) {
        int64_t v = q.begin()->second;
        q.erase(q.begin());

        for (const auto &edge : graph[v]) {
            int64_t to = edge.first;
            int64_t len = edge.second;
            if (d[v] + len < d[to]) {
                q.erase(std::make_pair(d[to], to));
                d[to] = d[v] + len;
                p[to] = v;
                q.insert(std::make_pair(d[to], to));
            }
        }
    }
}

// Флойд-Уоршелл взят с лекции
void floydWarshall(std::vector<std::vector<std::pair<int64_t, int64_t>>> graph) {
    size_t n = graph.size();
    std::vector<std::vector<int64_t>> d(n, std::vector<int64_t>(n, INT64_MAX));

    for (int64_t i = 0; i < n; ++i) {
        d[i][i] = 0;
        for (const auto &edge : graph[i]) {
            int64_t v = edge.first;
            int64_t weight = edge.second;
            d[i][v] = weight;
        }
    }

    for (int64_t k = 0; k < n; ++k) {
        for (int64_t i = 0; i < n; ++i) {
            for (int64_t j = 0; j < n; ++j) {
                if (d[i][k] < INT64_MAX && d[k][j] < INT64_MAX) {
                    d[i][j] = std::min(d[i][j], d[i][k] + d[k][j]);
                }
            }
        }
    }

    int64_t placeholder;
    for (int64_t i = 0; i < n; ++i) {
        for (int64_t j = 0; j < n; ++j) {
            if (i != j) {
                placeholder = i + j;
                if (d[i][j] == INT64_MAX) {
                    placeholder = -1;
                } else {
                    placeholder = d[i][j];
                }
            }
        }
    }
}

// Беллман Форд
void
fordBellman(const std::vector<std::vector<std::pair<int64_t, int64_t>>> &graph, int64_t source) {
    size_t numVertices = graph.size();
    std::vector<int64_t> distance(numVertices, INT64_MAX);
    distance[source] = 0;

    for (int64_t i = 0; i < numVertices - 1; ++i) {
        for (int64_t u = 0; u < numVertices; ++u) {
            for (const auto &edge : graph[u]) {
                int64_t v = edge.first;
                int64_t weight = edge.second;
                if (distance[u] != INT64_MAX && distance[u] + weight < distance[v]) {
                    distance[v] = distance[u] + weight;
                }
            }
        }
    }

    for (int64_t u = 0; u < numVertices; ++u) {
        for (const auto &edge : graph[u]) {
            int64_t v = edge.first;
            int64_t weight = edge.second;
            if (distance[u] != INT64_MAX && distance[u] + weight < distance[v]) {
                return;
            }
        }
    }
}

// A* algorithm
void aStar(const std::vector<std::vector<std::pair<int64_t, int64_t>>> &graph, int64_t start,
           int64_t goal) {
    size_t numNodes = graph.size();
    std::vector<int64_t> distance(numNodes, INT64_MAX);
    std::vector<int64_t> heuristic(numNodes, INT64_MAX);
    std::vector<int64_t> parent(numNodes, -1);
    std::vector<bool> visited(numNodes, false);

    distance[start] = 0;
    heuristic[start] = 0;

    auto compare = [&](int64_t a, int64_t b) {
        return (distance[a] + heuristic[a]) > (distance[b] + heuristic[b]);
    };

    std::priority_queue<int64_t, std::vector<int64_t>, decltype(compare)> pq(compare);
    pq.push(start);

    while (!pq.empty()) {
        int64_t currNode = pq.top();
        pq.pop();

        if (currNode == goal) {
            std::vector<int64_t> path;
            int64_t node = goal;
            while (node != -1) {
                path.push_back(node);
                node = parent[node];
            }
            std::reverse(path.begin(), path.end());
        }

        if (visited[currNode])
            continue;

        visited[currNode] = true;

        for (const auto &edge : graph[currNode]) {
            int64_t neighbor = edge.first;
            int64_t cost = edge.second;
            int64_t newDistance = distance[currNode] + cost;

            if (newDistance < distance[neighbor]) {
                distance[neighbor] = newDistance;
                heuristic[neighbor] = std::abs(goal - neighbor);
                parent[neighbor] = currNode;
                pq.push(neighbor);
            }
        }
    }
}


bool compareGraphs(const std::vector<std::vector<std::pair<int64_t, int64_t>>> &graph1,
                   const std::vector<std::vector<std::pair<int64_t, int64_t>>> &graph2) {
    size_t numEdges1 = 0;
    for (const auto &vertices : graph1) {
        numEdges1 += vertices.size();
    }

    size_t numEdges2 = 0;
    for (const auto &vertices : graph2) {
        numEdges2 += vertices.size();
    }

    return numEdges1 < numEdges2;
}

int64_t countEdges(const std::vector<std::vector<std::pair<int64_t, int64_t>>> &graph) {
    int64_t numEdges = 0;
    for (const auto &vertices : graph) {
        numEdges += vertices.size();
    }
    return numEdges;
}

int main() {
    std::ios_base::sync_with_stdio(false);
    std::cin.tie(nullptr);

    int64_t start_vertices = 10;
    int64_t end_vertices = 1010;
    int64_t step = 50;

    // Выбор выходного файла
    std::ofstream out;
    out.open("results2.csv", std::ios_base::app);

    // В данном коде сперва генерируются все графы, заполняя массивы графов, а далее уже начинается тестирование алгоритмов
    // Это вторая версия кода. В первой не было сортировки по рёбрам, и они в ней не считались. Вместо этого считались вершины

    // Генерация полных графов
    std::vector<std::vector<std::vector<std::pair<int64_t, int64_t>>>> completeGraphs = GraphGenerator::generateCompleteGraphs(
            start_vertices,
            end_vertices, step);

    // Генерация связных графов с коэффициентом плотности приблизительно 0.4-0.5
    double density = 0.4;
    std::vector<std::vector<std::vector<std::pair<int64_t, int64_t>>>> connectedGraphs = GraphGenerator::generateConnectedGraphs(
            start_vertices,
            end_vertices,
            step, density);

    // Генерация разреженных графов (деревьев)
    std::vector<std::vector<std::vector<std::pair<int64_t, int64_t>>>> sparseGraphs = GraphGenerator::generateSparseGraphs(
            start_vertices,
            end_vertices,
            step);

    std::sort(completeGraphs.begin(), completeGraphs.end(), compareGraphs);
    std::sort(connectedGraphs.begin(), connectedGraphs.end(), compareGraphs);
    std::sort(sparseGraphs.begin(), sparseGraphs.end(), compareGraphs);
    std::vector<std::vector<std::vector<std::vector<std::pair<int64_t, int64_t>>>>> all{
            completeGraphs, connectedGraphs, sparseGraphs};
    int64_t res, numEdges;
    for (int j = 0; j < 3; ++j) {
        setGrType(j);
        // Тестирование всех алгоритмов по очереди
        for (const auto &i : all[j]) {
            res = 0;
            for (int k = 0; k < 5; ++k) {
                numEdges = countEdges(i);
                auto start = std::chrono::high_resolution_clock::now();
                dijkstra(i, 0);
                auto elapsed = std::chrono::high_resolution_clock::now() - start;
                int64_t nanoseconds = std::chrono::duration_cast<std::chrono::nanoseconds>(
                        elapsed).count();
                res += nanoseconds / 5;
            }
            out << "Dijkstra" << ";Type:" << gr_type << ";Edges:" << numEdges << ";"
                << std::round(res / 1000) << std::endl;
            std::cout << "Dijkstra" << ";Type:" << gr_type << ";Edges:" << numEdges << ";"
                      << std::round(res / 1000) << std::endl;
        }
        for (const auto &i : all[j]) {
            res = 0;
            for (int k = 0; k < 5; ++k) {
                numEdges = countEdges(i);
                auto start = std::chrono::high_resolution_clock::now();
                floydWarshall(i);
                auto elapsed = std::chrono::high_resolution_clock::now() - start;
                int64_t nanoseconds = std::chrono::duration_cast<std::chrono::nanoseconds>(
                        elapsed).count();
                res += nanoseconds / 5;
            }
            out << "Floyd-Warshall" << ";Type:" << gr_type << ";Edges:" << numEdges << ";"
                << std::round(res / 1000) << std::endl;
            std::cout << "Floyd-Warshall" << ";Type:" << gr_type << ";Edges:" << numEdges << ";"
                      << std::round(res / 1000) << std::endl;
        }
        for (const auto &i : all[j]) {
            res = 0;
            for (int k = 0; k < 5; ++k) {
                numEdges = countEdges(i);
                auto start = std::chrono::high_resolution_clock::now();
                fordBellman(i, 0);
                auto elapsed = std::chrono::high_resolution_clock::now() - start;
                int64_t nanoseconds = std::chrono::duration_cast<std::chrono::nanoseconds>(
                        elapsed).count();
                res += nanoseconds / 5;
            }
            out << "Ford-Bellman" << ";Type:" << gr_type << ";Edges:" << numEdges << ";"
                << std::round(res / 1000) << std::endl;
            std::cout << "Ford-Bellman" << ";Type:" << gr_type << ";Edges:" << numEdges << ";"
                      << std::round(res / 1000) << std::endl;
        }
        for (const auto &i : all[j]) {
            res = 0;
            for (int k = 0; k < 5; ++k) {
                numEdges = countEdges(i);
                auto start = std::chrono::high_resolution_clock::now();
                aStar(i, 0, i.size() - 1);
                auto elapsed = std::chrono::high_resolution_clock::now() - start;
                int64_t nanoseconds = std::chrono::duration_cast<std::chrono::nanoseconds>(
                        elapsed).count();
                res += nanoseconds / 5;
            }
            out << "A-star" << ";Type:" << gr_type << ";Edges:" << numEdges
                << ";" << std::round(res / 1000) << std::endl;
            std::cout << "A-star" << ";Type:" << gr_type << ";Edges:" << numEdges
                      << ";" << std::round(res / 1000) << std::endl;
        }

    }

    return 0;
}
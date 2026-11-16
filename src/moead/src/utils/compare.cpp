#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include "../../headers/global_modules/generate_initial_population/generate_rSolution.h"
#include <iomanip>

using namespace std;

// Defina a função dominates
bool dominates(const Solution& sol1, const Solution& sol2) {
    return (sol1.fitness.first < sol2.fitness.first && sol1.fitness.second > sol2.fitness.second);
}

void compareSolutions(const string& fileA, const string& fileB) {
    ifstream inFileA(fileA);
    ifstream inFileB(fileB);
    
    if (!inFileA || !inFileB) {
        cerr << "Erro ao abrir os arquivos." << endl;
        return;
    }

    vector<Solution> solutionsA;
    vector<Solution> solutionsB;

    // Ler soluções do arquivo A
    double x, y;
    while (inFileA >> x >> y) {
        Solution sol;
        sol.fitness.first = x;
        sol.fitness.second = y;
        solutionsA.push_back(sol);
        cout << "Lendo solução do arquivo A: " << x << ", " << y << endl; // Debug
    }

    // Ler soluções do arquivo B
    while (inFileB >> x >> y) {
        Solution sol;
        sol.fitness.first = x;
        sol.fitness.second = y;
        solutionsB.push_back(sol);
        cout << "Lendo solução do arquivo B: " << x << ", " << y << endl; // Debug
    }

    inFileA.close();
    inFileB.close();

    // Comparar soluções de A com soluções de B
    for (const auto& solA : solutionsA) {
        int count = 0; // Contador de soluções dominadas
        for (const auto& solB : solutionsB) {
            if (dominates(solA, solB)) {
                count++;
            }
        }
        // Exibir resultado
        cout << "Solução " << solA.fitness.first << ", " << solA.fitness.second 
             << " domina " << count << " soluções do arquivo B." << endl;
    }
}

int main() {
    // Substitua pelos nomes dos arquivos conforme necessário
    string fileA = "islame.txt";
    string fileB = "GustavoJoao.txt";

    cout << fixed << setprecision(12);

    compareSolutions(fileA, fileB);

    return 0;
}

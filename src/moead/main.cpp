#include "./headers/main.h"
#include "./headers/globals.h"

#include <iostream>
#include <string>
using namespace std;

int countRevalue = 0;

BoundedParetoSet * pareto = new BoundedParetoSet();
int stop_criteria = 1000000;
string algorithm = "moead";
string instance = "A";
string root_folder = "./";
string crossover_type = "onepoint";
string mutation_type = "random";

int main(int argc, char* argv[]){

    if(argc == 2){
        instance = argv[1];
    } else if (argc > 2){
        instance = argv[1];
        root_folder = argv[2];
    }

    if(argc > 5){
        crossover_type = argv[5];
    }

    if(argc > 6){
        mutation_type = argv[6];
    }

    if(crossover_type != "onepoint" && crossover_type != "blockswap" && crossover_type != "blockcopy"){
        cerr << "ERROR: Invalid crossover type: " << crossover_type << " (use onepoint, blockswap or blockcopy)" << endl;
        return 1;
    }

    if(mutation_type != "random" && mutation_type != "worst" && mutation_type != "movement"){
        cerr << "ERROR: Invalid mutation type: " << mutation_type << " (use random, worst or movement)" << endl;
        return 1;
    }

    string path;

    int num_neighbors = 10;

    get_instance_info(argc, argv);

    int sum = 0;
    for(auto elem : turbines_per_zone)
        sum += elem;

    cout << endl;
    cout << "Number of subproblems: " << SIZE_OF_POPULATION << endl;
    cout << "Number of neighbors: " << num_neighbors << endl;
    cout << "Number of fixed turbines: " << fixd.size() << endl;
    cout << "Number of mobile turbines: " << sum << endl;
    cout << "Wind: " << wind << endl;
    cout << "Angle: " << angle << endl;
    cout << "Crossover: " << crossover_type << endl;
    cout << "Mutation: " << mutation_type << endl << endl;

    cout << "Run time:" << endl;
    
    auto population = create_initial_population(SIZE_OF_POPULATION);
    moead(population);
}
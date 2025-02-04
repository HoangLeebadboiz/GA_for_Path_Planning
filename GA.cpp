#include <iostream>
#include <vector>
#include <random>
#include <algorithm>
#include <string>
#include <fstream>
#include "include/json.hpp"
#include <map>
#include <ctime>

using json = nlohmann::json;

class PathPlanningGA {
private:
    struct Individual {
        std::vector<int> chromosome;
        double fitness;
    };

    int populationSize;
    int generations;
    double mutationRate;
    double crossoverRate;
    std::vector<Individual> population;
    std::map<std::pair<std::string, std::string>, double> costMap;
    std::vector<std::string> locationNames;
    std::string startLocation;

    // Convert location names to indices for easier manipulation
    std::map<std::string, int> locationToIndex;
    std::map<int, std::string> indexToLocation;

public:
    PathPlanningGA(int popSize = 100, int gens = 1000, double mutRate = 0.01, double crossRate = 0.8)
        : populationSize(popSize), generations(gens), mutationRate(mutRate), crossoverRate(crossRate) {
        std::srand(std::time(nullptr));
    }

    void loadCostData(const std::string& start, const std::vector<std::string>& destinations) {
        locationNames.clear();
        costMap.clear();
        locationToIndex.clear();
        indexToLocation.clear();
        
        // Add start location and destinations to locationNames
        locationNames.push_back(start);
        locationNames.insert(locationNames.end(), destinations.begin(), destinations.end());
        
        // Create index mappings
        for (size_t i = 0; i < locationNames.size(); ++i) {
            locationToIndex[locationNames[i]] = i;
            indexToLocation[i] = locationNames[i];
        }
        
        // Load costs between all pairs of locations
        for (size_t i = 0; i < locationNames.size(); i++) {
            for (size_t j = i + 1; j < locationNames.size(); j++) {
                const std::string& from = locationNames[i];
                const std::string& to = locationNames[j];
                
                double cost = getCostFromJson(from, to);
                costMap[{from, to}] = cost;
                costMap[{to, from}] = cost; // Assuming symmetric costs
            }
        }
        
        // Print loaded costs for verification
        std::cout << "Loaded costs:\n";
        for (const auto& cost : costMap) {
            std::cout << cost.first.first << " -> " 
                     << cost.first.second << ": " 
                     << cost.second << "\n";
        }
    }

    void initializePopulation(const std::string& start, const std::vector<std::string>& destinations) {
        startLocation = start;
        population.clear();
        
        // Create a base chromosome with destination indices
        std::vector<int> baseChromosome;
        for (const auto& dest : destinations) {
            baseChromosome.push_back(locationToIndex[dest]);
        }

        // Generate initial population
        for (int i = 0; i < populationSize; ++i) {
            Individual ind;
            ind.chromosome = baseChromosome;
            std::random_shuffle(ind.chromosome.begin(), ind.chromosome.end());
            ind.fitness = calculateFitness(ind.chromosome);
            population.push_back(ind);
        }
    }

    double calculateFitness(const std::vector<int>& chromosome) {
        double totalCost = 0.0;
        std::string currentLoc = startLocation;

        // Calculate total path cost
        for (int index : chromosome) {
            std::string nextLoc = indexToLocation[index];
            auto costPair = std::make_pair(currentLoc, nextLoc);
            
            if (costMap.find(costPair) != costMap.end()) {
                totalCost += costMap[costPair];
            } else {
                totalCost += 999999; // Large penalty for invalid paths
            }
            
            currentLoc = nextLoc;
        }

        return 1.0 / (totalCost + 1.0); // Convert cost to fitness (higher is better)
    }

    std::vector<int> crossover(const std::vector<int>& parent1, const std::vector<int>& parent2) {
        if (static_cast<double>(rand()) / RAND_MAX >= crossoverRate) {
            return parent1;
        }

        // Order Crossover (OX)
        int size = parent1.size();
        std::vector<int> child(size, -1);
        
        // Select random segment
        int start = rand() % size;
        int end = rand() % size;
        if (start > end) std::swap(start, end);

        // Copy segment from parent1
        for (int i = start; i <= end; ++i) {
            child[i] = parent1[i];
        }

        // Fill remaining positions with parent2's genes
        int j = 0;
        for (int i = 0; i < size; ++i) {
            if (i >= start && i <= end) continue;
            
            while (std::find(child.begin(), child.end(), parent2[j]) != child.end()) {
                j++;
            }
            child[i] = parent2[j++];
        }

        return child;
    }

    void mutate(std::vector<int>& chromosome) {
        for (size_t i = 0; i < chromosome.size(); ++i) {
            if (static_cast<double>(rand()) / RAND_MAX < mutationRate) {
                int j = rand() % chromosome.size();
                std::swap(chromosome[i], chromosome[j]);
            }
        }
    }

    std::vector<std::string> optimize() {
        for (int gen = 0; gen < generations; ++gen) {
            std::vector<Individual> newPopulation;

            // Elitism: keep the best individual
            auto elite = std::max_element(population.begin(), population.end(),
                [](const Individual& a, const Individual& b) { return a.fitness < b.fitness; });
            newPopulation.push_back(*elite);

            // Generate new population
            while (newPopulation.size() < populationSize) {
                // Tournament selection
                auto parent1 = tournament();
                auto parent2 = tournament();

                // Crossover
                auto childChromosome = crossover(parent1.chromosome, parent2.chromosome);

                // Mutation
                mutate(childChromosome);

                // Add new individual
                Individual child;
                child.chromosome = childChromosome;
                child.fitness = calculateFitness(childChromosome);
                newPopulation.push_back(child);
            }

            population = newPopulation;
        }

        // Get best solution
        auto best = std::max_element(population.begin(), population.end(),
            [](const Individual& a, const Individual& b) { return a.fitness < b.fitness; });

        // Convert solution back to location names
        std::vector<std::string> optimalPath;
        for (int index : best->chromosome) {
            optimalPath.push_back(indexToLocation[index]);
        }

        return optimalPath;
    }

private:
    Individual tournament(int tournamentSize = 5) {
        Individual best = population[rand() % populationSize];
        
        for (int i = 1; i < tournamentSize; ++i) {
            Individual contestant = population[rand() % populationSize];
            if (contestant.fitness > best.fitness) {
                best = contestant;
            }
        }
        
        return best;
    }

    double getCostFromJson(const std::string& from, const std::string& to) {
        std::string filename = "AStarResult/" + from + "To" + to + ".json";
        std::string altFilename = "AStarResult/" + to + "To" + from + ".json";
        
        std::ifstream f(filename);
        if (f.good()) {
            json data = json::parse(f);
            return data["cost"];
        }
        
        std::ifstream altF(altFilename);
        if (altF.good()) {
            json data = json::parse(altF);
            return data["cost"];
        }
        
        return 999999; // Return large cost if path not found
    }
};

// Example usage
int main() {
    PathPlanningGA ga;
    
    // Example: Start from "Bread" and visit other locations
    std::string startLocation = "Bread";
    std::vector<std::string> destinations = {"Gas", "Fruit", "Salt", "Oil", "Sugar"};

    ga.loadCostData(startLocation, destinations);
    ga.initializePopulation(startLocation, destinations);
    std::vector<std::string> optimalPath = ga.optimize();

    std::cout << "\nOptimal path from " << startLocation << ":\n";
    for (const auto& location : optimalPath) {
        std::cout << location << " -> ";
    }
    std::cout << "End\n";

    return 0;
}

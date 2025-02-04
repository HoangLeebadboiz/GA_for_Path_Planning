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

    void loadCostData() {
        std::vector<std::string> files = {
            "AStarResult/FruitToGas.json",
            "AStarResult/BreadToSalt.json",
            "AStarResult/BreadToOil.json",
            "AStarResult/BreadToGas.json",
            "AStarResult/BreadToFruit.json"
        };

        for (const auto& file : files) {
            std::ifstream f(file);
            json data = json::parse(f);
            
            // Extract locations from filename
            std::string filename = file.substr(file.find_last_of("/") + 1);
            std::string from = filename.substr(0, filename.find("To"));
            std::string to = filename.substr(filename.find("To") + 2, filename.find(".json") - filename.find("To") - 2);
            
            // Store cost
            costMap[{from, to}] = data["cost"];
            costMap[{to, from}] = data["cost"]; // Assuming symmetric costs
            
            // Add locations to the set if not already present
            if (std::find(locationNames.begin(), locationNames.end(), from) == locationNames.end()) {
                locationNames.push_back(from);
            }
            if (std::find(locationNames.begin(), locationNames.end(), to) == locationNames.end()) {
                locationNames.push_back(to);
            }
        }

        // Create index mappings
        for (size_t i = 0; i < locationNames.size(); ++i) {
            locationToIndex[locationNames[i]] = i;
            indexToLocation[i] = locationNames[i];
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
};

// Example usage
int main() {
    PathPlanningGA ga;
    ga.loadCostData();

    // Example: Start from "Bread" and visit other locations
    std::string startLocation = "Bread";
    std::vector<std::string> destinations = {"Gas", "Fruit", "Salt", "Oil"};

    ga.initializePopulation(startLocation, destinations);
    std::vector<std::string> optimalPath = ga.optimize();

    std::cout << "Optimal path from " << startLocation << ":\n";
    for (const auto& location : optimalPath) {
        std::cout << location << " -> ";
    }
    std::cout << "End\n";

    return 0;
}

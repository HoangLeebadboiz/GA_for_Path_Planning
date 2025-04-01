#include <iostream>
#include <vector>
#include <random>
#include <algorithm>
#include <string>
#include <fstream>
#include "include/json.hpp"
#include <map>
#include <ctime>
#include <cmath>
#include <filesystem>
#include <iomanip>
#include <limits>


using json = nlohmann::json;
namespace fs = std::filesystem;

// Structure to represent a good in the supermarket
struct Good {
    double x, y;
    std::string name;
    Good(double x = 0, double y = 0, const std::string& name = "") : x(x), y(y), name(name) {}
};

// Calculate Euclidean distance between two points
double distance(const Good& p1, const Good& p2) {
    return std::sqrt(std::pow(p1.x - p2.x, 2) + std::pow(p1.y - p2.y, 2));
}

class KMeans {
private:
    int k; // Number of clusters
    std::vector<Good> goods; // Input goods
    std::vector<Good> centroids; // Cluster centroids
    std::vector<int> assignments; // Cluster assignments for each good
    double inertia; // Sum of squared distances to nearest centroid
    double maxRadius; // Maximum allowed radius for each cluster
    
public:
    KMeans(int k, double radius = std::numeric_limits<double>::max()) : k(k), inertia(0.0), maxRadius(radius) {}
    
    // Add a good to the dataset
    void addGood(const Good& g) {
        goods.push_back(g);
    }
    
    // Get centroid at specified index
    Good getCentroid(int index) const {
        if (index >= 0 && index < k) {
            return centroids[index];
        }
        return Good();
    }
    
    // K-means++ initialization
    void initializeCentroidsPlusPlus() {
        std::random_device rd;
        std::mt19937 gen(rd());

        // Choose first centroid randomly
        std::uniform_int_distribution<> dis(0, goods.size() - 1);
        centroids.clear();
        centroids.push_back(goods[dis(gen)]);

        // Choose remaining centroids
        for (int i = 1; i < k; i++) {
            std::vector<double> distances(goods.size());
            double sum = 0.0;

            // Calculate distances to nearest existing centroid
            for (size_t j = 0; j < goods.size(); j++) {
                double minDist = std::numeric_limits<double>::max();
                for (const auto& centroid : centroids) {
                    minDist = std::min(minDist, distance(goods[j], centroid));
                }
                distances[j] = minDist * minDist;
                sum += distances[j];
            }

            // Choose next centroid with probability proportional to D²
            std::uniform_real_distribution<> dis_prob(0.0, sum);
            double rand_val = dis_prob(gen);
            double cumsum = 0.0;
            size_t chosen_idx = 0;

            for (size_t j = 0; j < goods.size(); j++) {
                cumsum += distances[j];
                if (cumsum >= rand_val) {
                    chosen_idx = j;
                    break;
                }
            }

            centroids.push_back(goods[chosen_idx]);
        }
    }
    
    // Assign goods to nearest centroid
    bool assignClusters() {
        bool changed = false;
        std::vector<int> newAssignments(goods.size());
        double new_inertia = 0.0;
        
        for (size_t i = 0; i < goods.size(); ++i) {
            double minDist = std::numeric_limits<double>::max();
            int closestCentroid = 0;
            
            for (int j = 0; j < k; ++j) {
                double dist = distance(goods[i], centroids[j]);
                if (dist < minDist) {
                    minDist = dist;
                    closestCentroid = j;
                }
            }
            
            // Add penalty if distance exceeds maxRadius
            if (minDist > maxRadius) {
                new_inertia += std::numeric_limits<double>::max(); // Large penalty
            } else {
                new_inertia += minDist * minDist;
            }
            
            newAssignments[i] = closestCentroid;
            if (assignments.empty() || newAssignments[i] != assignments[i]) {
                changed = true;
            }
        }
        
        assignments = newAssignments;
        inertia = new_inertia;
        return changed;
    }
    
    // Update centroid positions
    void updateCentroids() {
        std::vector<Good> newCentroids(k, Good(0, 0));
        std::vector<int> counts(k, 0);
        
        for (size_t i = 0; i < goods.size(); ++i) {
            int cluster = assignments[i];
            newCentroids[cluster].x += goods[i].x;
            newCentroids[cluster].y += goods[i].y;
            counts[cluster]++;
        }
        
        for (int i = 0; i < k; ++i) {
            if (counts[i] > 0) {
                newCentroids[i].x /= counts[i];
                newCentroids[i].y /= counts[i];
                newCentroids[i].name = "centroid" + std::to_string(i);
            }
        }
        
        centroids = newCentroids;
    }
    
    double getInertia() const {
        return inertia;
    }
    
    // Run k-means clustering
    void run(int maxIterations = 100) {
        if (goods.empty() || k <= 0) return;
        
        initializeCentroidsPlusPlus();
        
        for (int iter = 0; iter < maxIterations; ++iter) {
            bool changed = assignClusters();
            if (!changed) break;
            updateCentroids();
        }
    }
    
    // Get clustering results as array of arrays with good names
    void printClusterResults() {
        std::vector<std::vector<std::string>> clusters(k);
        for (size_t i = 0; i < goods.size(); ++i) {
            clusters[assignments[i]].push_back(goods[i].name);
        }
        
        std::cout << "\nClustering Results:\n";
        std::cout << "==================\n";
        for (size_t i = 0; i < clusters.size(); ++i) {
            std::cout << "Cluster " << i + 1 << ": ";
            for (size_t j = 0; j < clusters[i].size(); ++j) {
                std::cout << clusters[i][j];
                if (j < clusters[i].size() - 1) std::cout << ", ";
            }
            std::cout << std::endl;
        }
        std::cout << "Inertia: " << inertia << std::endl;
        std::cout << "==================\n";
    }

    // Get all goods in each cluster
    std::vector<std::vector<std::string>> getClusterContents() const {
        std::vector<std::vector<std::string>> clusters(k);
        for (size_t i = 0; i < goods.size(); ++i) {
            clusters[assignments[i]].push_back(goods[i].name);
        }
        return clusters;
    }
};

// Function to read goods data from JSON files
std::vector<Good> readGoodsFromFiles(const std::string& folderPath) {
    std::vector<Good> goods;
    
    // Check if directory exists, if not create it
    if (!fs::exists(folderPath)) {
        std::cout << "Creating directory: " << folderPath << std::endl;
        fs::create_directory(folderPath);
        
        // Create sample data
        std::vector<std::pair<std::string, std::pair<double, double>>> samples = {
            {"Water", {1.0, 2.0}},
            {"Bread", {1.5, 2.5}},
            {"Chips", {10.0, 12.0}},
            {"Soda", {9.5, 11.5}},
            {"Candy", {5.0, 5.0}}
        };
        
        for (const auto& sample : samples) {
            json j;
            j["positon"] = {sample.second.first, sample.second.second, 0.0};
            j["orientation"] = {0.0, 0.0, 0.0, 1.0};
            
            std::string filename = folderPath + "/" + sample.first + ".json";
            std::ofstream file(filename);
            file << std::setw(4) << j << std::endl;
        }
        
        std::cout << "Created sample data files in " << folderPath << std::endl;
    }
    
    try {
        for (const auto& entry : fs::directory_iterator(folderPath)) {
            if (entry.path().extension() == ".json") {
                std::ifstream file(entry.path());
                if (file.is_open()) {
                    json j;
                    file >> j;
                    
                    if (j.contains("positon") && j["positon"].is_array() && j["positon"].size() >= 2) {
                        double x = j["positon"][0];
                        double y = j["positon"][1];
                        std::string name = entry.path().stem().string();
                        goods.emplace_back(x, y, name);
                    }
                }
            }
        }
    } catch (const fs::filesystem_error& e) {
        std::cerr << "Error accessing directory: " << e.what() << std::endl;
        return goods;
    }
    
    return goods;
}

// Function to find optimal k using elbow method
int findOptimalK(const std::vector<Good>& goods, double maxRadius = std::numeric_limits<double>::max(), int maxK = 10) {
    std::vector<double> inertias;
    
    // Try different values of k
    for (int k = 1; k <= maxK; ++k) {
        KMeans kmeans(k, maxRadius);
        for (const auto& good : goods) {
            kmeans.addGood(good);
        }
        kmeans.run();
        inertias.push_back(kmeans.getInertia());
    }
    
    // Find the elbow point using the maximum curvature
    int optimalK = 1;
    double maxCurvature = 0.0;
    
    for (int i = 1; i < maxK - 1; ++i) {
        double curvature = std::abs(inertias[i-1] - 2*inertias[i] + inertias[i+1]);
        if (curvature > maxCurvature) {
            maxCurvature = curvature;
            optimalK = i + 1;
        }
    }
    
    // Print elbow curve
    std::cout << "\nElbow Curve:\n";
    std::cout << "============\n";
    for (int i = 0; i < maxK; ++i) {
        std::cout << "K=" << i+1 << ": Inertia=" << inertias[i] << std::endl;
    }
    std::cout << "============\n";
    
    return optimalK;
}

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
    std::string endLocation;
    bool hasEndLocation = false;

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

        // Calculate path through destinations
        for (int index : chromosome) {
            std::string nextLoc = indexToLocation[index];
            auto costPair = std::make_pair(currentLoc, nextLoc);
            
            if (costMap.find(costPair) != costMap.end()) {
                totalCost += costMap[costPair];
            } else {
                totalCost += 999999;
            }
            
            currentLoc = nextLoc;
        }

        // Add cost to end location if specified
        if (hasEndLocation) {
            auto finalCostPair = std::make_pair(currentLoc, endLocation);
            if (costMap.find(finalCostPair) != costMap.end()) {
                totalCost += costMap[finalCostPair];
            } else {
                totalCost += 999999;
            }
        }

        return 1.0 / (totalCost + 1.0);
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

        // Add end location to path if specified
        if (hasEndLocation) {
            optimalPath.push_back(endLocation);
        }

        return optimalPath;
    }

    void setEndPosition(const std::string& end) {
        endLocation = end;
        hasEndLocation = true;
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

// Function to select goods by name
std::vector<Good> selectGoods(const std::vector<Good>& allGoods, const std::vector<std::string>& selectedNames) {
    std::vector<Good> selectedGoods;
    for (const auto& name : selectedNames) {
        auto it = std::find_if(allGoods.begin(), allGoods.end(),
            [&name](const Good& g) { return g.name == name; });
        if (it != allGoods.end()) {
            selectedGoods.push_back(*it);
        }
    }
    return selectedGoods;
}

// Example usage
int main() {
    // 1. Read goods from Pose folder
    std::vector<Good> allGoods = readGoodsFromFiles("ProductPoseReal");
    
    // Print available goods
    std::cout << "\nAvailable goods:\n";
    for (const auto& good : allGoods) {
        std::cout << good.name << " ";
    }
    std::cout << "\n\n";
    
    // 2. Select goods to include in clustering
    std::vector<std::string> selectedGoods = {
        "Bread", "Cooker", "Fruit", "Gas", "MicrowaveOven", "Oil", "Refrigerator", "Salt", "Snack", "SoftDrink", "Sugar", "Television", "Water"
    };
    
    // Specify start and end locations
    std::string startLocation = "Start";  // Start location
    std::string endLocation = "CheckoutCounter";      // End location
    
    // Get selected goods (excluding start and end locations)
    std::vector<Good> goods = selectGoods(allGoods, selectedGoods);
    
    if (goods.empty()) {
        std::cerr << "No valid goods selected. Exiting.\n";
        return 1;
    }
    
    std::cout << "\nSelected goods for clustering:\n";
    for (const auto& good : goods) {
        std::cout << good.name << " ";
    }
    std::cout << "\n\n";
    
    // 3. Find optimal number of clusters with radius limit
    double maxClusterRadius = 5.0; // Set your desired maximum cluster radius
    int optimalK = findOptimalK(goods, maxClusterRadius);
    // int optimalK = 3; // For testing, set a fixed number of clusters    
    std::cout << "Optimal number of clusters: " << optimalK << std::endl;
    
    // 4. Run K-means clustering without radius limit
    KMeans kmeans(optimalK);  // No radius limit for fixed k
    for (const auto& good : goods) {
        kmeans.addGood(good);
    }
    kmeans.run();
    kmeans.printClusterResults();
    
    // 5. Get cluster centers and cluster contents
    std::vector<std::string> clusterCenters;
    std::vector<std::vector<std::string>> clusterContents = kmeans.getClusterContents();
    
    for (int i = 0; i < optimalK; i++) {
        Good center = kmeans.getCentroid(i);
        // Find nearest good to this center
        double minDist = std::numeric_limits<double>::max();
        std::string nearestGood;
        
        for (const auto& good : goods) {
            double dist = distance(center, good);
            if (dist < minDist) {
                minDist = dist;
                nearestGood = good.name;
            }
        }
        clusterCenters.push_back(nearestGood);
    }
    
    // 6. Run GA with start location, cluster centers, and end location
    PathPlanningGA ga;
    
    // Use cluster centers as intermediate destinations
    std::vector<std::string> destinations = clusterCenters;
    
    // Run GA with specified start and end
    ga.loadCostData(startLocation, destinations);
    ga.initializePopulation(startLocation, destinations);
    ga.setEndPosition(endLocation);
    std::vector<std::string> optimalPath = ga.optimize();

    // Create 1D array of goods arranged by clusters
    std::vector<std::string> arrangedGoods;
    
    // Map to track which cluster each good belongs to
    std::map<std::string, int> goodToCluster;
    for (size_t i = 0; i < optimalPath.size(); i++) {
        const std::string& center = optimalPath[i];
        // Find cluster index for this center
        for (size_t j = 0; j < clusterCenters.size(); j++) {
            if (clusterCenters[j] == center) {
                // Add all goods from this cluster
                std::vector<std::string> clusterGoods = clusterContents[j];
                // Randomly shuffle goods within cluster
                std::random_shuffle(clusterGoods.begin(), clusterGoods.end());
                // Add to arranged goods
                arrangedGoods.insert(arrangedGoods.end(), clusterGoods.begin(), clusterGoods.end());
                break;
            }
        }
    }

    // Print results
    std::cout << "\nPath Planning Results:\n";
    std::cout << "=====================\n";
    std::cout << "1. Optimal Path:\n";
    std::cout << startLocation << " -> ";
    for (const auto& location : optimalPath) {
        std::cout << location << " -> ";
    }
    std::cout << endLocation << "\n\n";

    std::cout << "2. Arranged Goods (by cluster order):\n";
    std::cout << "[ ";
    for (size_t i = 0; i < arrangedGoods.size(); i++) {
        std::cout << arrangedGoods[i];
        if (i < arrangedGoods.size() - 1) {
            std::cout << ", ";
        }
    }
    std::cout << " ]\n";
    std::cout << "=====================\n";

    return 0;
}

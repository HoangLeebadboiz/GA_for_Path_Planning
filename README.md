# Path Planning Genetic Algorithm

A C++ implementation of a Genetic Algorithm (GA) for optimizing path planning with multiple destinations. The program finds the optimal order to visit multiple locations given a starting point, using pre-calculated path costs stored in JSON files.

## Features

- Genetic Algorithm implementation for path optimization
- JSON-based cost data loading
- Configurable GA parameters (population size, generations, mutation rate, etc.)
- Support for multiple destinations
- Automatic cost loading from JSON files

## Table of Contents

- [Path Planning Genetic Algorithm](#path-planning-genetic-algorithm)
  - [Features](#features)
  - [Table of Contents](#table-of-contents)
  - [About](#about)
  - [Prerequisites](#prerequisites)
  - [Project Structure](#project-structure)
  - [Installation](#installation)
  - [Usage](#usage)
  - [Algorithm Details](#algorithm-details)
    - [Genetic Algorithm Components](#genetic-algorithm-components)
    - [Example Output](#example-output)
    - [Performance](#performance)

## About

This project implements a Genetic Algorithm to solve the path planning problem, similar to the Traveling Salesman Problem (TSP). Given a starting location and multiple destinations, it finds the optimal order to visit all destinations while minimizing the total path cost. The path costs between locations are pre-calculated and stored in JSON files.

## Prerequisites

- C++ compiler with C++11 support (g++ recommended)
- JSON library (nlohmann/json)
- Pre-calculated path costs in JSON format
- Git (optional, for cloning the repository)

## Project Structure

```
project_folder/
├── GA.cpp                  # Main GA implementation
├── include/
│   └── json.hpp           # JSON library header
├── AStarResult/           # Path cost JSON files
│   ├── FruitToGas.json
│   ├── BreadToSalt.json
│   ├── BreadToOil.json
│   ├── BreadToGas.json
│   └── BreadToFruit.json
└── README.md
```

## Installation

1. Clone or download the repository:
```bash
git clone <repository-url>
cd <project-folder>
```

2. Download the JSON library:
```bash
mkdir -p include
curl -o include/json.hpp https://raw.githubusercontent.com/nlohmann/json/develop/single_include/nlohmann/json.hpp
```

3. Compile the program:
```bash
g++ GA.cpp -o GA -std=c++11
```

## Usage

1. Prepare your JSON files in the `AStarResult` folder. Each file should follow the format:
```json
{
    "cost": 123.45
}
```
File naming convention: `{FromLocation}To{ToLocation}.json`

2. Run the program:
```bash
./GA
```

3. Modify the main function in GA.cpp to change start location and destinations:
```cpp
std::string startLocation = "Bread";
std::vector<std::string> destinations = {"Gas", "Fruit", "Salt", "Oil"};
```

4. Configure GA parameters (optional):
```cpp
PathPlanningGA ga(
    100,    // population size
    1000,   // generations
    0.01,   // mutation rate
    0.8     // crossover rate
);
```

## Algorithm Details

### Genetic Algorithm Components

1. **Chromosome Representation**
   - Each chromosome is a permutation of destination indices
   - Represents visit order of destinations

2. **Fitness Function**
   - Based on total path cost
   - Fitness = 1 / (totalCost + 1)
   - Higher fitness means better solution

3. **Genetic Operators**
   - Selection: Tournament selection
   - Crossover: Order Crossover (OX)
   - Mutation: Swap mutation
   - Elitism: Preserves best solution

### Example Output

```
Loaded costs:
Bread -> Gas: 100
Bread -> Fruit: 150
...

Optimal path from Bread:
Gas -> Fruit -> Salt -> Oil -> End
```

The output shows:
1. All loaded path costs between locations
2. The optimized path order starting from the initial location

### Performance

- Default parameters are optimized for small to medium-sized problems (5-15 destinations)
- Larger problems may require parameter tuning
- Typical runtime: few seconds for default parameters

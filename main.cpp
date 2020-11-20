#pragma clang diagnostic push
#pragma ide diagnostic ignored "cert-err58-cpp"

#include <iostream>
#include <fstream>
#include <cfloat>
#include <vector>
#include <random>
#include <chrono>

using namespace std;

mt19937 g_RNG(chrono::steady_clock::now().time_since_epoch().count());

// To handle all the data and generate chromosomes and values
class dataHandler {
private:
    vector<double> lowerBounds, upperBounds, channelROIs;
    vector<string> channelNames;
    double marketingBudget;
    int nChannels;
public:
    dataHandler(const vector<double> &lowerBounds, const vector<double> &upperBounds, const vector<double> &channelRoIs,
                const vector<string> &channelNames, double marketingBudget, int nChannels) :
            lowerBounds(lowerBounds), upperBounds(upperBounds), channelROIs(channelRoIs), channelNames(channelNames),
            marketingBudget(marketingBudget), nChannels(nChannels) {}

    int getNChannels() const { return nChannels; }

    double getChannelROI(int idx) const { return channelROIs[idx]; }

    double getMarketingBudget() const { return marketingBudget; };

    string getChannelName(int idx) const { return channelNames[idx]; }

    /// generates a random investment for a specific channel (0 - base) based on the channel lower and upper bounds
    double generateInvestment(int channelIdx) {
        // nextafter -> gets the value after upperBounds[channelIdx] to make the UB inclusive
        uniform_real_distribution<double> dist(lowerBounds[channelIdx], nextafter(upperBounds[channelIdx], DBL_MAX));
        return dist(g_RNG);
    }


};

class Chromosome {
private:
    dataHandler *algorithmsData = nullptr;
    vector<double> chromosomeData;
    double fitness = DBL_MIN;
    bool isFitnessUptoDate = false;

    void uniformMutate() {
        // TODO: implement uniform mutation
        // dont update fitness
    }

    void nonuniformMutate() {
        // TODO: implement non-uniform mutation
        // dont update fitness here
    }

    // fitness is the total profit
    void setFitness() {
        // TODO: calculate and set fitness here
        isFitnessUptoDate = true;
    }

public:
    Chromosome() = default;

    explicit Chromosome(dataHandler *algorithmsData) {
        this->algorithmsData = algorithmsData;

        // initialize the chromosome randomly
        chromosomeData.resize(algorithmsData->getNChannels());
        for (int i = 0; i < chromosomeData.size(); ++i) {
            chromosomeData[i] = algorithmsData->generateInvestment(i);
        }
    }

    double getData(int idx) { return chromosomeData[idx]; }

    void mutate() {
        // randomly selects between the two types
        std::uniform_real_distribution<> dist(0, 1);
        int selection = (int) round(dist(g_RNG));

        if (selection) uniformMutate();
        else nonuniformMutate();

        isFitnessUptoDate = false;
    }

    pair<Chromosome, Chromosome> combineWith(const Chromosome &other) const {
        if (this->chromosomeData == other.chromosomeData) return {*this, other};

        Chromosome ret1 = *this, ret2 = other;
        // TODO: implement 2 points cross-over and fill ret1 and ret2

        return {ret1, ret2};
    }

    double getFitness() {
        if (!isFitnessUptoDate)
            setFitness();
        return fitness;
    }

    void printData(ostream &out) {
        for (int i = 0; i < chromosomeData.size(); ++i)
            out << algorithmsData->getChannelName(i) << " -> " << this->getData(i) << "K\n";
        out << "\nThe total profit is " << this->getFitness() << "K\n";
    }

};

class GeneticAlgorithm {
private:
    dataHandler *data;
    Chromosome bestChromosome;
    vector<Chromosome> population;

    /// gets the population fitness
    vector<double> getChromosomesFitness() {
        vector<double> fitness(population.size());
        for (int i = 0; i < fitness.size(); ++i)
            fitness[i] = population[i].getFitness();

        return fitness;
    }

    void replaceOffspring(int oldIdx, Chromosome &newData) {
        population[oldIdx] = newData;

        if (newData.getFitness() > bestChromosome.getFitness())
            bestChromosome = newData;
    }

    void nextGeneration() {
        vector<double> chromosomesFitness = getChromosomesFitness();
        int firstParentIdx = 0, secondParentIdx = data->getNChannels() - 1;
        /* TODO: implement tournament selection and fill (firstParentIdx and secondParentIdx),
         * After implementing remove firstParentIdx and secondParentIdx initialization
         */


        Chromosome firstParent = population[firstParentIdx], secondParent = population[secondParentIdx];
        pair<Chromosome, Chromosome> offSprings = firstParent.combineWith(secondParent);

        offSprings.first.mutate();
        offSprings.second.mutate();

        replaceOffspring(firstParentIdx, offSprings.first);
        replaceOffspring(secondParentIdx, offSprings.second);

    }

    void logState(int stateID, ostream &loggerStream) {
        /* To become more fancy add this line of code
         * if(stateID == 1) loggerStream << "Initial State\n";
         * else loggerStream << stateID << "\n";
         */
        // TODO: log the current state to the stream
    }

public:
    GeneticAlgorithm(dataHandler *data, int populationSize) {
        this->data = data;
        population.resize(populationSize);

        population[0] = Chromosome(data);
        bestChromosome = population[0];

        for (int i = 1; i < populationSize; ++i) {
            population[i] = Chromosome(data);

            if (population[i].getFitness() > bestChromosome.getFitness())
                bestChromosome = population[i];
        }
    }

    Chromosome execute(int epochs, ostream *loggerStream = nullptr) {
        for (int currState = 1; currState <= epochs; ++currState) {
            if (loggerStream)
                logState(currState, *loggerStream);
            nextGeneration();
        }
        (*loggerStream) << "\nEND OF THE EXECUTION\n------------------------------------------\n";
        return bestChromosome;
    }

};

int main() {
    double marketingBudget, channelROI, channelLB, channelUB;
    int nChannels;
    string channelName, tempLB, tempUB;


    string outputFilename = "logger.out";
    auto *outputFile = new ofstream(); // to redirect it to cout -> `auto* outputFile = &cout`
    outputFile->open(outputFilename, ios::app); // If you want to clear the output every run make it ios::out

    cout << "Enter the marketing budget (in thousands \"K\"):\n";
    cin >> marketingBudget;

    cout << "Enter the number of marketing channels:\n";
    cin >> nChannels;

    vector<double> lowerBounds(nChannels), upperBounds(nChannels), channelROIs(nChannels);
    vector<string> channelNames(nChannels);

    cout << "Enter the name and ROI (in %) of each channel separated by space:\n";
    for (int i = 0; i < nChannels; ++i) {
        cin >> channelName >> channelROI;

        channelNames[i] = channelName;
        channelROIs[i] = channelROI;
    }

    cout << "Enter the lower (k) and upper bounds (%) of investment in each channel: (enter x if there is no bound)\n";
    for (int i = 0; i < nChannels; ++i) {
        cin >> tempLB >> tempUB;

        try {
            channelLB = stod(tempLB);
        } catch (...) {
            channelLB = 0;
        }
        try {
            channelUB = stod(tempUB);
        } catch (...) {
            channelUB = 100;
        }


        lowerBounds[i] = channelLB;
        upperBounds[i] = channelUB;
    }

    cout << "Please wait while running the GA…\n\n";

    // pass this by reference to other objects
    auto *data = new dataHandler(lowerBounds, upperBounds, channelROIs, channelNames, marketingBudget, nChannels);
    int populationSize = 20, epochs = 100; // change these parameters for different outputs

    GeneticAlgorithm algorithmModel(data, populationSize);
    Chromosome result = algorithmModel.execute(epochs, outputFile);

    cout << "The final marketing budget allocation is:\n";
    result.printData(cout);

    if(outputFile->is_open()) outputFile->close();
    return 0;

}
/* Input:
100
4
8
12
7
11
2.7 58
20.5 x
x 18
10 x
 */


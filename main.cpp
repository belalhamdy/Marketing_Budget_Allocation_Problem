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
    bool isInBounds(int idx, double val) const {
        return (lowerBounds[idx] <= val && upperBounds[idx] >= val);
    }
    //Forces val to be in bounds of the idx-th constraint.
    double clipToBounds(int idx, double val) const {
        if (val < lowerBounds[idx]) val = lowerBounds[idx];
        if (val > upperBounds[idx]) val = upperBounds[idx];
        return val;
    }

};

/*
 * Chromosome carries the investment percentage for each channel
 * When calculating fitness and representing data we multiply the investment by the marketing budget
 */
class Chromosome {
private:
    dataHandler *algorithmsData = nullptr;
    vector<double> chromosomeData;
    double fitness = DBL_MIN;
    bool isFitnessUptoDate = false;

    double getChannelReturn(int idx) {
        double investment = this->getData(idx) / 100 * algorithmsData->getMarketingBudget();
        return (investment * algorithmsData->getChannelROI(idx));
    }
    // fitness is the total profit
    void setFitness() {
        double totalFitness = 0;
        for (int i = 0; i < chromosomeData.size(); ++i)
            totalFitness += getChannelReturn(i);

        fitness = totalFitness;
        isFitnessUptoDate = true;
    }

    bool isFeasible() {
        double totalPercentages = 0;
        bool isFeasible = true;
        for (int i = 0; i < chromosomeData.size(); ++i) {
            totalPercentages += chromosomeData[i];
            isFeasible &= algorithmsData->isInBounds(i, chromosomeData[i]);
        }
        isFeasible &= (totalPercentages <= 100);
        return isFeasible;
    }

    void makeChromosomeFeasible() {
        if (!isFeasible()) {
            double totalPercentages = 0;
            for (int i = 0; i < chromosomeData.size(); ++i) {
                chromosomeData[i] = algorithmsData->clipToBounds(i, chromosomeData[i]);
                totalPercentages += chromosomeData[i];
            }
            if (totalPercentages > 100) {
                for (int i = 0; i < chromosomeData.size(); ++i) {
                    double newVal = algorithmsData->clipToBounds(i, 0.0);
                    totalPercentages -= newVal - chromosomeData[i];
                    chromosomeData[i] = newVal;
                    if (totalPercentages <= 100)
                        break;
                }
            }
            isFitnessUptoDate = false;
        }
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

        makeChromosomeFeasible();
    }

    double getData(int idx) { return chromosomeData[idx]; }

    void mutate(double probability) {
        std::uniform_real_distribution<double> doMutation(0, 1);
        if (doMutation(g_RNG) < probability) return;

        std::uniform_int_distribution<int> mutationIdx(0, chromosomeData.size());
        int i = mutationIdx(g_RNG);

        chromosomeData[i] = algorithmsData->generateInvestment(i);

        makeChromosomeFeasible();
    }

    pair<Chromosome, Chromosome> combineWith(const Chromosome &other) const {
        if (this->chromosomeData == other.chromosomeData) return {*this, other};

        Chromosome ret1 = *this, ret2 = other;

        int chromosomeSz = (int) chromosomeData.size();

        uniform_int_distribution<int> uid(0, chromosomeSz - 1);

        //With 2 point cross-over at L and R we get three parts: [0, L-1][L, R][R+1, N]
        //In the off-chance that L and R are the same, the algorithm changes to 1-point cross-over.

        int L = uid(g_RNG), R = uid(g_RNG);
        if (R == L) {
            for (int i = R + 1; i < chromosomeSz; ++i)
                std::swap(ret1.chromosomeData[i], ret2.chromosomeData[i]);
        } else {
            if (R < L) std::swap(L, R);
            for (int i = L; i <= R; ++i)
                std::swap(ret1.chromosomeData[i], ret2.chromosomeData[i]);
        }

        return {ret1, ret2};
    }

    double getFitness() {
        if (!isFitnessUptoDate)
            setFitness();
        return fitness;
    }


    void printData(ostream &out) {
        for (int i = 0; i < chromosomeData.size(); ++i) {
            out << algorithmsData->getChannelName(i) << " -> " << getChannelReturn(i) << "K\n";
        }
        out << "\nThe total profit is " << this->getFitness() << "K\n";
    }

};

class GeneticAlgorithm {
private:
    dataHandler *data;
    Chromosome bestChromosome;
    vector<Chromosome> population;
    double mutationProbability;
    bool uniformMutation;

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

    void nextGeneration(int generationIdx) {
        vector<double> chromosomesFitness = getChromosomesFitness();

        std::uniform_int_distribution<int> uid(0, data->getNChannels() - 1);
        auto Tournament = [this](int option1, int option2) {
            return population[option1].getFitness() > population[option2].getFitness() ? option1 : option2;
        };
        int firstParentIdx = Tournament(uid(g_RNG), uid(g_RNG));
        int secondParentIdx = Tournament(uid(g_RNG), uid(g_RNG));

        Chromosome firstParent = population[firstParentIdx], secondParent = population[secondParentIdx];
        pair<Chromosome, Chromosome> offSprings = firstParent.combineWith(secondParent);

        double probability;
        if (uniformMutation)
            probability = mutationProbability;
        else
            probability = mutationProbability / generationIdx;

        offSprings.first.mutate(probability);
        offSprings.second.mutate(probability);

        replaceOffspring(firstParentIdx, offSprings.first);
        replaceOffspring(secondParentIdx, offSprings.second);

    }

    void logState(int stateID, ostream &loggerStream) {
        if (stateID == 1) loggerStream << "Initial State(" << stateID << ")\n";
        else loggerStream << stateID << "\n";

        loggerStream << "Population:\n";
        for (auto &curr : population) {
            curr.printData(loggerStream);
        }
        loggerStream << "\nBest Chromosome:\n";
        bestChromosome.printData(loggerStream);
        loggerStream << "\n";
    }

public:
    GeneticAlgorithm(dataHandler *data, int populationSize, double mutationProbability = 0.1,
                     bool uniformMutation = true) {
        this->data = data;
        population.resize(populationSize);

        population[0] = Chromosome(data);
        bestChromosome = population[0];

        for (int i = 1; i < populationSize; ++i) {
            population[i] = Chromosome(data);

            if (population[i].getFitness() > bestChromosome.getFitness())
                bestChromosome = population[i];
        }
        this->mutationProbability = mutationProbability;
        this->uniformMutation = uniformMutation;
    }

    Chromosome execute(int epochs, ostream *loggerStream = nullptr) {
        if (loggerStream)
            (*loggerStream) << "\nSTART OF THE EXECUTION\n------------------------------------------\n"
                            << "Number of Epochs = " << epochs << "\n"
                            << "Population Size = " << population.size();

        for (int currState = 1; currState <= epochs; ++currState) {
            if (loggerStream)
                logState(currState, *loggerStream);

            nextGeneration(currState);
        }

        if (loggerStream)
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

    cout << "Enter the name (without any spaces) and ROI (in %) of each channel separated by space:\n";
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

    if (outputFile->is_open())
        outputFile->close();

    return 0;

}
/* Input:
100
4
TV_Advertisement 8
Google 12
Twitter 7
Facebook 11
2.7 58
20.5 x
x 18
10 x

*/


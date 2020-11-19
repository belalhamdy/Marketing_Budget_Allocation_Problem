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

int main() {
    double marketingBudget, channelROI, channelLB, channelUB;
    int nChannels;
    string channelName, tempLB, tempUB;

    string outputFilename = "logger.out";
    ofstream outputFile;
    outputFile.open(outputFilename, ios::out);
    // TODO: send the ofstream to the function that runs the algorithm to log the outputs to this file and don't forget to close the file

    cout << "Enter the marketing budget (in thousands \"K\"):\n";
    cin >> marketingBudget;

    cout << "Enter the number of marketing channels:\n";
    cin >> nChannels;

    cout << "Enter the name and ROI (in %) of each channel separated by space:\n";
    for (int i = 0; i < nChannels; ++i) {
        cin >> channelName >> channelROI;
        // TODO: Create a data structure to store these values (channelName,channelROI)
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
            channelUB = DBL_MAX;
        }

        // TODO: Create a data structure to store these values (channelLB,channelUB)
    }

    cout << "Please wait while running the GA…\n\n";
    // TODO: Do processing here
    cout << "The final marketing budget allocation is:\n";
    // TODO: Output the solution


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


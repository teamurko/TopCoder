#include <iostream>
#include <vector>
#include <string>

#include "detector.h"

using namespace std;

int main(int argc, char *argv[])
{
    int N;
    cin >> N;
    vector<int> imageData(N);
    for (int i=0;i<N;i++)
        cin >> imageData[i];
    int R, C;
    cin >> R >> C;

    OmegaDetector theDetector;
    vector<int> result = theDetector.ROI_Finder(imageData, R, C);
    cout << result.size() << endl;
    for (int i=0; i < result.size(); i++)
        cout << result[i] << endl;
    cout.flush();

    return 0;
}

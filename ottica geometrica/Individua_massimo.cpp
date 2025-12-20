#include <iostream>
#include <fstream>
#include <vector>
#include <limits>
using namespace std;

int main() {
    string filename;
    cout << "Inserisci nome file dati: ";
    cin >> filename;

    ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Errore: impossibile aprire il file\n";
        return 1;
    }

    vector<int> pixel;
    vector<double> value;

    int p;
    double v;
    while (file >> p >> v) {
        pixel.push_back(p);
        value.push_back(v);
    }
    file.close();

    int p_min, p_max;
    cout << "Inserisci pixel minimo: ";
    cin >> p_min;
    cout << "Inserisci pixel massimo: ";
    cin >> p_max;

    double max_value = -numeric_limits<double>::infinity();
    int max_pixel = -1;

    for (size_t i = 0; i < pixel.size(); ++i) {
        if (pixel[i] >= p_min && pixel[i] <= p_max) {
            if (value[i] > max_value) {
                max_value = value[i];
                max_pixel = pixel[i];
            }
        }
    }

    if (max_pixel == -1) {
        cout << "Nessun dato trovato nell'intervallo selezionato\n";
        return 1;
    }

    cout << "Massimo trovato:\n";
    cout << "Pixel = " << max_pixel << "\n";
    cout << "Valore = " << max_value << "\n";

    return 0;
}

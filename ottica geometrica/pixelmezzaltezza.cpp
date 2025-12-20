#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

double trovaMassimo(const std::vector<double>& dati) {
    if (dati.empty()) {
        throw std::invalid_argument("Il vettore è vuoto");
    }

    double maxVal = dati[0];
    for (double val : dati) {
        if (val > maxVal) {
            maxVal = val;
        }
    }
    return maxVal;
}



int main() {
    string filename = "dati2.txt";   // file esportato da ImageJ
    double frazione = 0.5;           // 50%

    vector<double> x, I;

    // ===== LETTURA FILE =====
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Errore: impossibile aprire il file\n";
        return 1;
    }

    double xp, Ip;
    while (file >> xp >> Ip) {
        x.push_back(xp);
        I.push_back(Ip);
    }
    file.close();

    if (x.size() < 3) {
        cerr << "Errore: dati insufficienti\n";
        return 1;
    }

    // ===== INTENSITÀ MASSIMA =====
    double Imax = trovaMassimo(I);
    double Ith = frazione * Imax;

    // ===== TROVA INDICI SOPRA SOGLIA =====
    int iL = -1, iR = -1;
    for (int i = 1; i < I.size(); i++) {
        if (I[i-1] < Ith && I[i] >= Ith && iL == -1) {
            iL = i;
        }
        if (I[i-1] >= Ith && I[i] < Ith) {
            iR = i - 1;
        }
    }

    if (iL == -1 || iR == -1) {
        cerr << "Errore: il profilo non interseca il 50%\n";
        return 1;
    }

    // ===== INTERPOLAZIONE LINEARE =====
    auto interp = [&](int i1, int i2) {
        if (I[i2] - I[i1] == 0) return x[i1]; // evita divisione per zero
        return x[i1] + (Ith - I[i1]) * (x[i2] - x[i1]) / (I[i2] - I[i1]);
    };

    double xL = interp(iL-1, iL);
    double xR = interp(iR, iR+1);

    double width = xR - xL;

    // ===== OUTPUT =====
    cout << "Larghezza al 50% (FWHM) = " << width << " pixel" << endl;

    return 0;
}


#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cmath>

using namespace std;

struct ParabolaFit {
    double a, b, c;  // y = ax² + bx + c
};

// Risoluzione sistema 3x3 con sostituzione diretta (metodo di Cramer)
ParabolaFit fitParabola(const vector<double>& x, const vector<double>& y) {
    if (x.size() != y.size() || x.size() < 3) {
        cerr << "Errore: Dati insufficienti per il fit" << endl;
        return {0, 0, 0};
    }

    double Sx=0, Sx2=0, Sx3=0, Sx4=0, Sy=0, Sxy=0, Sx2y=0;
    int N = x.size();

    for (int i=0; i<N; i++) {
        double xi = x[i];
        double xi2 = xi*xi;
        Sx += xi;
        Sx2 += xi2;
        Sx3 += xi2*xi;
        Sx4 += xi2*xi2;
        Sy += y[i];
        Sxy += xi*y[i];
        Sx2y += xi2*y[i];
    }

    // Costruisco matrice dei coefficienti
    double M[3][4] = {
        { Sx4, Sx3, Sx2, Sx2y },
        { Sx3, Sx2, Sx , Sxy  },
        { Sx2, Sx , N  , Sy   }
    };

    // Eliminazione gaussiana
    for (int i=0; i<3; i++) {
        // Pivoting
        double pivot = M[i][i];
        for (int j=i; j<4; j++) M[i][j] /= pivot;
        for (int k=0; k<3; k++) {
            if (k==i) continue;
            double factor = M[k][i];
            for (int j=i; j<4; j++)
                M[k][j] -= factor*M[i][j];
        }
    }

    // Soluzione
    double a = M[0][3];
    double b = M[1][3];
    double c = M[2][3];

    return {a, b, c};
}

int main() {
    string filename;
    cout << "Inserisci il nome del file con i dati: ";
    cin >> filename;

    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Errore: impossibile aprire il file" << endl;
        return 1;
    }

    vector<double> omega_f, theta;
    double w, t;
    while (file >> w >> t) {
        omega_f.push_back(w);
        theta.push_back(t);
    }

    if (omega_f.size() < 3) {
        cerr << "Errore: Troppo pochi punti nel file per il fit" << endl;
        return 1;
    }

    // Trova indice del massimo di theta
    auto it_max = max_element(theta.begin(), theta.end());
    size_t idx_max = distance(theta.begin(), it_max);
    double omega_f_max = omega_f[idx_max];

    // Determina intervallo di punti intorno al massimo
    int start = static_cast<int>(idx_max) - 6;
    int end = static_cast<int>(idx_max) + 6;

    if (start < 0) start = 0;
    if (end >= static_cast<int>(theta.size())) end = theta.size() - 1;
    
    cout << "\nIntervallo di fit:" << endl;
cout << "Primo punto: omega_f = " << omega_f[start] << ", theta = " << theta[start] << endl;
cout << "Ultimo punto: omega_f = " << omega_f[end] << ", theta = " << theta[end] << endl;


    vector<double> omega_f_fit, theta_fit;
    for (int i = start; i <= end; i++) {
        omega_f_fit.push_back(omega_f[i] - omega_f_max); // centratura
        theta_fit.push_back(theta[i]);
    }

    if (omega_f_fit.size() < 3) {
        cerr << "Errore: Troppo pochi punti per il fit" << endl;
        return 1;
    }

    cout << "\nDati usati per il fit (omega centrata):" << endl;
    for (size_t i = 0; i < omega_f_fit.size(); i++) {
        cout << omega_f_fit[i] << "\t" << theta_fit[i] << endl;
    }

    ParabolaFit fit = fitParabola(omega_f_fit, theta_fit); 
    // Scrivo i valori della parabola su file fit.dat
ofstream outfile("fit.dat");
if (!outfile.is_open()) {
    cerr << "Errore: impossibile creare il file fit.dat" << endl;
    return 1;
}

// Intervallo di omega_f per il grafico
double delta_omega_min = omega_f[start] - omega_f[idx_max];
double delta_omega_max = omega_f[end] - omega_f[idx_max];
double passo = 0.001;

// Scrivo i punti del fit riportati in ?f assoluto
for (double delta_omega = delta_omega_min; delta_omega <= delta_omega_max; delta_omega += passo) {
    double theta_fit_val = fit.a * delta_omega * delta_omega + fit.b * delta_omega + fit.c;
    double omega = delta_omega + omega_f[idx_max];  // riconversione a wf
    outfile << omega << " " << theta_fit_val << endl;
}
outfile.close();
cout << "\nFile fit.dat generato con il fit parabolico." << endl;


    cout << "\nRisultato del fit parabolico:" << endl;
    cout << "theta = " << fit.a << " w^2 + " << fit.b << " w + " << fit.c << endl;

    if (fit.a != 0) {
        double delta_omega_R = -fit.b / (2 * fit.a);
        double omega_R = delta_omega_R + omega_f_max;
        double erroreomega_R = (fabs(omega_f[end]- omega_f[start]))/(sqrt(12)); 
        cout << "\nPulsazione di risonanza stimata: w_R = " << omega_R << endl;
        cout << "L'errore di w_R e' " << erroreomega_R << endl;
        if (fit.a < 0)
            cout << "(Massimo di risonanza verificato)" << endl;
        else
            cout << "Attenzione: il fit non ha trovato un massimo!" << endl;
    } else {
        cout << "Attenzione: coefficiente a nullo, fit non valido!" << endl;
    }

    return 0;
}


#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <string>
#include <sstream>

using namespace std;

struct ParabolaFit {
    double a, b, c;      // θ = a·ω'² + b·ω' + c
    double sigma_a;      // errore su a
    double sigma_b;      // errore su b
    double sigma_c;      // errore su c
};

ParabolaFit fitParabola(const vector<double>& x, const vector<double>& y) {
    ParabolaFit res{0,0,0, 0,0,0};
    int N = x.size();
    if (N != (int)y.size() || N < 3) {
        cerr << "Errore: Dati insufficienti per il fit" << endl;
        return res;
    }

    // Calcolo delle somme
    double S0 = N;
    double Sx = 0,   Sx2 = 0,  Sx3 = 0,  Sx4 = 0;
    double Sy = 0,  Sxy = 0, Sx2y = 0;
    for (int i = 0; i < N; i++) {
        double xi = x[i];
        double xi2 = xi*xi;
        Sx   += xi;
        Sx2  += xi2;
        Sx3  += xi2*xi;
        Sx4  += xi2*xi2;
        Sy   += y[i];
        Sxy  += xi * y[i];
        Sx2y += xi2 * y[i];
    }

    // Matrice aumentata 3x4 per Gauss
    double M[3][4] = {
        {Sx4, Sx3, Sx2, Sx2y},
        {Sx3, Sx2, Sx,  Sxy },
        {Sx2, Sx,  S0,  Sy  }
    };

    // Eliminazione di Gauss per risolvere M * [a b c]^T = rhs
    for (int i = 0; i < 3; i++) {
        // pivot
        double piv = M[i][i];
        if (fabs(piv) < 1e-12) {
            cerr << "Errore: Matrice singolare" << endl;
            return res;
        }
        for (int j = i; j < 4; j++)
            M[i][j] /= piv;
        for (int k = 0; k < 3; k++) {
            if (k == i) continue;
            double fac = M[k][i];
            for (int j = i; j < 4; j++)
                M[k][j] -= fac * M[i][j];
        }
    }
    // Estrazione dei parametri
    res.a = M[0][3];
    res.b = M[1][3];
    res.c = M[2][3];

    // Calcolo chi² sui punti di fit
    double chi2 = 0;
    for (int i = 0; i < N; i++) {
        double yi_fit = res.a*x[i]*x[i] + res.b*x[i] + res.c;
        double r = y[i] - yi_fit;
        chi2 += r*r;
    }
    double var_resid = chi2 / (N - 3);  // stima della varianza dei residui

    // Costruisco la matrice normale H (3x3)
    double H00 = Sx4, H01 = Sx3, H02 = Sx2;
    double H11 = Sx2, H12 = Sx;
    double H22 = S0;

    // Calcolo determinante di H
    double detH = 
        H00*(H11*H22 - H12*H12)
      - H01*(H01*H22 - H12*H02)
      + H02*(H01*H12 - H11*H02);

    // Inversa diagonale di H (covarianze senza cov. incrociate)
    double invH00 =  (H11*H22 - H12*H12) / detH;
    double invH11 =  (H00*H22 - H02*H02) / detH;
    double invH22 =  (H00*H11 - H01*H01) / detH;

    // Deviazioni standard parametri
    res.sigma_a = sqrt(invH22 * var_resid);
    res.sigma_b = sqrt(invH11 * var_resid);
    res.sigma_c = sqrt(invH00 * var_resid);

    return res;
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
    string line;
    // Salta intestazione
    getline(file, line);

    // Leggo w, devw, t, devt, nomefile
    while (true) {
        double w, devw, t, devt;
        string fname;
        if (!(file >> w >> devw >> t >> devt >> fname)) break;
        omega_f.push_back(w);
        theta   .push_back(t);
    }
    file.close();

    int M = omega_f.size();
    if (M < 3) {
        cerr << "Errore: troppi pochi punti (" << M << ")" << endl;
        return 1;
    }

    // Cerco il massimo di theta
    auto it_max = max_element(theta.begin(), theta.end());
    int idx_max = distance(theta.begin(), it_max);
    double w_max = omega_f[idx_max];
    double t_max = *it_max;

    cout << "Massimo trovato: θ=" << t_max << " a ω=" << w_max << endl;

    // Seleziono 13 punti intorno al massimo (6 prima, 6 dopo)
    int start = max(0, idx_max - 6);
    int end   = min(M-1, idx_max + 6);

    vector<double> x_fit, y_fit;
    for (int i = start; i <= end; i++) {
        x_fit.push_back(omega_f[i] - w_max);  // centratura su w_max
        y_fit.push_back(theta[i]);
    }

    cout << "Uso " << x_fit.size() << " punti per il fit parabolico\n";

    // Eseguo il fit
    ParabolaFit fit = fitParabola(x_fit, y_fit);

    // Calcolo ω_R e la sua incertezza tramite propagazione (no covarianza)
    double a = fit.a, b = fit.b;
    double da = fit.sigma_a, db = fit.sigma_b;
    double omega_R = -b / (2*a) + w_max;
    // ∂ω_R/∂a =  b/(2 a²),  ∂ω_R/∂b = -1/(2 a)
    double dR_da =  b / (2 * a*a);
    double dR_db = -1.0 / (2 * a);
    double sigma_omega_R = sqrt( (dR_da*da)*(dR_da*da)
                              + (dR_db*db)*(dR_db*db) );

    // Stampi il risultato
    cout << "\n--- Risultati fit parabolico ---\n";
    cout << "a = " << fit.a << "con incertezza" << fit.sigma_a << "\n";
    cout << "b = " << fit.b << " con incertezza" << fit.sigma_b << "\n";
    cout << "c = " << fit.c << " con incertezza " << fit.sigma_c << "\n";

    cout << "\nPulsazione di risonanza:\n";
    cout << "omega risonanza = " << omega_R
         << " con incertezza" << sigma_omega_R << "\n";

    return 0;
}

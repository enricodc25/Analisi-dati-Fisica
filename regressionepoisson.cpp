#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>
// #include <boost/math/distributions/chi_squared.hpp>
using namespace std;

//prototipi
void interpolazionePoisson(vector<double> x,vector<double> y, int cov);
double DeviazioneStandardPosteriori(double a, double b, vector<double> x,vector<double> y);
double CoefficienteCorrelazione(vector<double> x, vector<double> y);
double chiquadroLineare(vector<double> x, vector <double> y, vector<double> sigmay, double a, double b);
// double PValueChiquadro(double chiQuadrato, int gradiLiberta); //ha bisogno dell'iinstallazione della libreria Boost nel terminale in alternativa si usa l'approssimazione funzione gamma
double p_value_chi_quadrato(double x, int k);
double gamma_incompleta(double s, double z);

int main() {
	string filename;
    vector<double> x, y;
	double xval,yval,sigma;
    double n1,n2,n3,n4;
    int cov=0;
	
	cout<<"___________REGRESSIONE LINEARE PESATA POISSONIANA RAGGI GAMMA______ "<<endl;
    //operazioni dal file e presa dati
	cout << "Inserire il nome del file di input: ";
    cin >> filename;
    ifstream fin(filename);

    if(!fin){
    	cout<<"Il file non e' stato caricato correttamente"<<endl;
    	return -1;
	}
    cout<<"Inserisci i due intervalli di valori delle X dove vuoi fare l'interpolazione:  "<<endl;
    cout<<" Da ";
    cin>>n1;
    cout<<" a "; 
    cin>>n2;
    cout<<"mentre il secondo intervallo: ";
    cout<<" Da ";
    cin>>n3;
    cout<<" a ";
    cin>>n4;

	while(fin>>xval>>yval){
		if(xval>=n1 && xval<=n2 || xval>=n3 && xval<=n4){
            x.push_back(xval);    		
            y.push_back(yval);
        }
    }
    cout<<"I due parametri a e b sono covarianti? ( inserisci 1 se si , mentre 0 se non lo sono ): "<<endl;
    cin>>cov;
    interpolazionePoisson(x,y,cov);

    return 0;
}


double DeviazioneStandardPosteriori(double a, double b, vector<double> x,vector<double> y) { 
    double N = x.size();
    double somma_quadrati_residui = 0.0;
    // Calcola la somma dei quadrati dei residui
    for (size_t i = 0; i < N; ++i) {
        double residuo = y[i] - a - b * x[i];
        somma_quadrati_residui += residuo * residuo;
    }
    // Calcola sigma posteriori
    return sqrt(somma_quadrati_residui / (N - 2));
}

double CoefficienteCorrelazione(vector<double> x, vector<double> y) {    
    const  double N = x.size();
    double somma_x = 0.0, somma_y = 0.0;
    double somma_x2 = 0.0, somma_y2 = 0.0;
    double somma_prodotti = 0.0;
    
    // Calcola le medie
    for (int i = 0; i < N; i++) {
        somma_x += x[i];
        somma_y += y[i];
    }
    double media_x = somma_x / N;
    double media_y = somma_y / N;
    
    // Calcola le somme necessarie per la formula
    for (size_t i = 0; i < N; i++) {
        double diff_x = x[i] - media_x;
        double diff_y = y[i] - media_y;
        
        somma_prodotti += diff_x * diff_y;
        somma_x2 += diff_x * diff_x;
        somma_y2 += diff_y * diff_y;
    }
    
    // Calcola il coefficiente di correlazione
    double denominatore = sqrt(somma_x2 * somma_y2);

    if (denominatore == 0.0) {
        return 0.0;
    }
    return somma_prodotti / denominatore;
}

double chiquadroLineare(vector<double> x, vector<double> y, vector<double> sigmay, double a, double b){
    double chitot=0.0;
    //I vettori x e y hanno la stessa dimensione
    for(int i=0;i<y.size();i++){
        chitot+=pow(y.at(i)-a-b*x.at(i),2)/pow(sigmay.at(i),2);
    }

    return chitot;
}

// double PValueChiquadro(double chiQuadrato, int gradiLiberta) {
//     //Bisogna aver installato nel terminale la libreria boost
//     // Crea una distribuzione chi-quadrato con il numero specificato di gradi di libertà
//     boost::math::chi_squared distribuzione(gradiLiberta);

//     // Calcola il p-value come la probabilità cumulativa complementare
//     return 1 - boost::math::cdf(distribuzione, chiQuadrato);
// }


//------------------------APPROSSIMAZIONE p-value  PER N° gradi di libertà elevati------------------------------------//
// Funzione Gamma incompleta inferiore
double gamma_incompleta(double s, double z) {
    const int max_iter = 100;
    const double epsilon = 1e-6;

    double sum = 1.0 / s;
    double term = sum;
    for (int n = 1; n < max_iter; ++n) {
        term *= z / (s + n);
        sum += term;
        if (std::fabs(term) < epsilon * sum) break;
    }
    return sum * std::exp(-z + s * std::log(z) - std::lgamma(s));
}

// Calcolo del p-value CON L'APPOSSIMAZIONE DELLA FUNZIONE GAMMA
double p_value_chi_quadrato(double x, int k) {
    double s = k / 2.0;
    double z = x / 2.0;
    double gamma_incomp = gamma_incompleta(s, z);
    double gamma_completa = std::tgamma(s);
    return 1.0 - gamma_incomp / gamma_completa;
}

//____________________________________________________________________________________________________________________________________//





void interpolazionePoisson(vector<double> x,vector<double> y, int cov){
    double spesi = 0.0;   
    double sx = 0.0;   
    double sy = 0.0;  
    double sxx = 0.0;  
    double sxy = 0.0;  
	vector<double> sigmay;

    for(int i=0;i<y.size();i++){
        sigmay.push_back(sqrt(y.at(i)));
    }

    //attenzione che nel calcolo del peso potrebbero esserci alcuni problemi se il dato della y è zero
	//calcolo somme intermedie
    for (int i = 0; i < x.size(); i++) {
        double w = 1.0 / (sigmay.at(i) * sigmay.at(i)); // peso = 1/s_i^2
        spesi  += w; //sommatoria pesi 
        sx  += x.at(i) * w; 
        sy  += y.at(i)* w;
        sxx += x.at(i) * x.at(i) * w;
        sxy += x.at(i) * y.at(i) * w;
    }

    // Calcolo del DELTA
    double delta = (spesi * sxx) - (sx * sx);

    // Stima dei parametri a e b
    double a = ( (sxx * sy) - (sx * sxy) ) / delta;
    double b = ( (spesi  * sxy) - (sx * sy ) ) / delta;

    //Calcolo delle incertezze relative ai parametri   
    double sigma_a = sqrt(sxx/ delta);
    double sigma_b = sqrt(spesi/ delta);
    double covarianza=0;
    if(cov==1){
        covarianza=-(sx)/delta;
    }


    // Stampa dei risultati
    cout << "=== RISULTATI DEL FIT LINEARE PESATO DI TIPO POISSONIANO per una retta y=a+bx ===" << endl;
    cout << "a = " << a << endl;
    cout << "b = " << b << endl;
    cout << "sigma_a = " << sigma_a <<endl;
    cout << "sigma_b = " << sigma_b << endl;
    cout << "sigma posteriori = "<<DeviazioneStandardPosteriori(a,b,x,y)<<endl;
    if(cov==1){
        cout<<"I due parametri a e b sono covarianti con una covarianza pari a: "<<covarianza<<endl;
    }
    cout<<"-------------------------TEST STATISTICI----------------------"<<endl;
    cout <<"coefficiente di correlazione r = "<<CoefficienteCorrelazione(x,y)<<endl;
    cout<<"Il X^2 con "<<(x.size()-2)<<" gradi di liberta' X^2= "<<chiquadroLineare(x,y,sigmay,a,b)<<endl;
    cout<<"Con un p-value= "<<p_value_chi_quadrato(chiquadroLineare(x,y,sigmay,a,b), (x.size()-2))<<endl<<endl;
    
}

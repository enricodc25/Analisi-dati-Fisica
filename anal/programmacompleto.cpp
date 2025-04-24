#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>
using namespace std;

//prototipi
void interpolazionePoisson(vector<double> x,vector<double> y, int cov,double& a, double& b, double& sigma_a,
    double&sigma_b, double& covarianza);
double DeviazioneStandardPosteriori(double& a, double& b, vector<double> x,vector<double> y);
double CoefficienteCorrelazione(vector<double> x, vector<double> y);
double chiquadroLineare(vector<double> x, vector <double> y, vector<double> sigmay, double& a, double& b);
// double PValueChiquadro(double chiQuadrato, int gradiLiberta); //ha bisogno dell'iinstallazione della libreria Boost nel terminale in alternativa si usa l'approssimazione funzione gamma
// double p_value_chi_quadrato(double x, int k);
// double gamma_incompleta(double s, double z);

int main() {
	string filename;
    vector<double> x, y;
	double xval,yval,sigma;
    double n1,n2,n3,n4;
    int cov=0;
    double a,b;
    double sigma_a, sigma_b, covarianza;
    vector<double> x2, y2; 
   
	
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
		if(xval>=n1 && xval<=n2 || xval>=n3 && xval<=n4 ){
			if(yval!=0){
				x.push_back(xval);    		
          			y.push_back(yval);
			}
        }
    }
    
    cout<<"I due parametri a e b sono covarianti? ( inserisci 1 se si , mentre 0 se non lo sono ): "<<endl;
    cin>>cov;
    interpolazionePoisson(x,y,cov,a,b,sigma_a, sigma_b, covarianza);

    
    
    double m1,m2;
    
    cout<<"Inserisci l'intervallo di valori delle X di cui vuoi calcolare l'integrale':  "<<endl;
    cout<<" Da ";
    cin>>m1;
    cout<<" a "; 
    cin>>m2;
    ifstream fin2(filename);
	while(fin2>>xval>>yval){
		if(xval>=m1 && xval<=m2){
				if(yval!=0){
				x2.push_back(xval);    		
          		y2.push_back(yval);
			}
        }
    }
    
    double Y = 0.0; //segnale complessivo//
    for (int i = 0; i < y2.size(); i++) {
        Y += y2[i];
    }
    cout << "L'area di Y vale:" << Y <<endl;
    
  double F= 0.0; //segnale di fondo//
    F= (m2-m1)*(a+b*((m2+m1)/2));
    cout<< "Il contributo di fondo F vale " << F <<endl; 
    
	double S=0.0; //segnale del picco//
	S= Y-F; 
	// la varianza di Y è uguale a Y stessa//
	cout << "Il segnale di picco vale " << S <<endl;
	
	double sigma2_F= 0.0; //incertezza al quadrato del segnale di fondo//
	sigma2_F= pow((m2-m1), 2)* pow(sigma_a,2)+pow((m2-m1), 2)* pow(((m2+m1)/2),2)* pow(sigma_b,2)+ 2* pow((m2-m1),2)* ((m2+m1)/2)*covarianza;
	cout << "La sigma al quadrato del segnale di fondo vale " << sigma2_F << endl;
	
	double sigma_S=0.0; 
	sigma_S= sqrt(Y+sigma2_F);
	
	double lambda=0.0;
	lambda= (S/sigma_S);
	cout << "La compatibilita' vale " << lambda <<endl;
	
	if(lambda>3){
		cout<< "La compatibilita' non e' accettabile,quindi il picco e' significativo";
	}
	else if (lambda>=2&&lambda <=3){
		cout<<"La compatibilita' e' sufficiente, quindi il picco protrebbe non essere significativo";
	}
	else if(lambda>= 1&&lambda<2){
		cout<< "La compatibilita' e' buona ";
	}
	else if(lambda<1){
		cout<< "La compatibilita' e' ottima, il picco non e' significativo ";
	}

double alpha = 1- erf (lambda / sqrt (2));
	cout << "alpha = " << alpha << endl;

	
return 0;    
}


double DeviazioneStandardPosteriori(double& a, double& b, vector<double> x,vector<double> y) { 
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
double chiquadroLineare(vector<double> x, vector<double> y, vector<double> sigmay, double& a, double& b){
    double chitot=0.0;
    //I vettori x e y hanno la stessa dimensione
    for(int i=0;i<y.size();i++){
        chitot+=pow(y.at(i)-a-b*x.at(i),2)/pow(sigmay.at(i),2);
    }

    return chitot;
}

void interpolazionePoisson(vector<double> x,vector<double> y, int cov, double& a, double& b, double& sigma_a,
    double&sigma_b, double& covarianza){
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
     a = ( (sxx * sy) - (sx * sxy) ) / delta;
     b = ( (spesi  * sxy) - (sx * sy ) ) / delta;

    //Calcolo delle incertezze relative ai parametri   
    sigma_a = sqrt(sxx/ delta);
    sigma_b = sqrt(spesi/ delta);
    covarianza=0;
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
    cout <<"coefficiente di correlazione r = "<<CoefficienteCorrelazione(x,y)<<endl;
    if(cov==1){
        cout<<"I due parametri a e b sono covarianti con una covarianza pari a: "<<covarianza<<endl;
    }
cout<<"-------------------------TEST STATISTICI----------------------"<<endl;
    cout <<"coefficiente di correlazione r = "<<CoefficienteCorrelazione(x,y)<<endl;
    cout<<"Il X^2 con "<<(x.size()-2)<<" gradi di liberta' X^2= "<<chiquadroLineare(x,y,sigmay,a,b)<<endl;
}

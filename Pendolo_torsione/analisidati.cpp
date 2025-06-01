#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>        
#include <cmath>
using namespace std;

struct ParabolaFitResult {
    double a, b, c;
};

struct Punto {
    double x;
    double y;
};

// dichiaro qui, definisco in fondo
// ParabolaFitResult fitParabola(const vector<double>& x, const vector<double>& y);
Punto trovaVertice(const ParabolaFitResult& parabola);
int estraiPositivi(const vector<double>& p,const vector<double>& t,int startIndex,vector<double>& out_p,vector<double>& out_t);
double media(vector<double> v);

int main() {
    //  vector<string> nomifile={"f0.9.txt","f0.905.txt","f0.910.txt","f0.915.txt","f0.925.txt", "f0.930.txt", "f0.935.txt", "f0.940.txt", "f0.945.txt", "f0.964.txt","f0.965.txt","f0.970.txt"
    //  ,"f0.975.txt", "f0.980.txt", "f0.06470.txt", "f0.9200.txt", "f0.9625.txt", "f0.9900.txt","f0.96450.txt", "f0.96460.txt"};
         vector<double> t, f, p, a, fa;
         string nomefile;
         cout<<"Inserisci il nome del file: ";
         cin>>nomefile;
         double M_PI= 3.14159226535897932384; 
    

    ifstream fin(nomefile);
    if (!fin) {
        cerr << "Errore: impossibile aprire il file\n";
        return 1;
    }
    double tv, fv, pv, av, fav;
    while (fin >> tv >> fv >> pv >> av >> fav) {
        t.push_back(tv);
        f.push_back(fv);
        p.push_back(pv);
        a.push_back(av);
        fa.push_back(fav);
    }
    fin.close();

    // ofstream pic("picchicandidati.txt");
    // if (!pic) {
    //     cerr << "Errore: impossibile aprire picchicandidati.txt\n";
    //     return 1;
    // }

    vector<double> bloccoP, bloccoT;   
    double idxraccolta = 0.0;
   	vector<double> xsFit, ysFit; 
    vector<double> stimax1,stimamax;
    double stimax,stimam;    
	vector <double> maxtotali;  
	while (idxraccolta<p.size()){
		 int len = estraiPositivi(p, t, int(idxraccolta), bloccoP, bloccoT); //mettere tipo int(...) significa fare un casting, ovvero "forzare" quel dato di essere di quel tipo primitivo
                                                                   			 //essendo un indice va messo int, nonostante noi lo ricaviamo come double
        if (len == 0) {
            idxraccolta += 1.0;
            continue;   
        }
		double maxVal = *max_element(bloccoP.begin(), bloccoP.end());
		maxtotali.push_back(maxVal);
		idxraccolta +=double(len);
	}
	double mediaMaxtotali = media(maxtotali);
	double sogliadeimax = mediaMaxtotali / 3.0;  // Soglia della media dei massimi
	
	
	double idx=0.0;
    while (idx < p.size()) {
        //questa funzione mi permette di prendere tutti i dati positivi fino a quando non incontra un dato negativo ( significa fine della mia curva) 
        int len = estraiPositivi(p, t, int(idx), bloccoP, bloccoT); //mettere tipo int(...) significa fare un casting, ovvero "forzare" quel dato di essere di quel tipo primitivo
                                                                    //essendo un indice va messo int, nonostante noi lo ricaviamo come double
        if (len == 0) {
            idx += 1.0;
            continue;
        }
        //l'idea del programma è andare a selezionare ogni curva positiva e predere come soglia il 70% rispetto al valore massimo globale in quell'intervallo
        double maxVal = *max_element(bloccoP.begin(), bloccoP.end()); //questa è una funzione della libreria algorithm, usando i metodi .begin() e .end() io mi assicuro di prendere dall'inizio alla fine, è una scrittura compatta
       	        
       	if (maxVal < sogliadeimax) {  // Salta picchi troppo piccoli
        idx += double(len);
        continue;
    	}
         double soglia = 0.70 * maxVal;
        cout << "Bloc #" << idx << ": max=" << maxVal //occhio, questa istruzione continua sotto, per questo nonostante sembra senza ; funziona il programma
             << " -> soglia=" << soglia << "\n"; //il barra \n è un modo derivante dal c per andare a capo, come con endl
        xsFit.clear();
        ysFit.clear(); //qua devo ripulire i vector

        for (int i = 0; i < bloccoP.size(); ++i) {
            if (bloccoP.at(i) >= soglia) {
                ysFit.push_back(bloccoP.at(i));
                xsFit.push_back(bloccoT.at(i));
            }
        } 
        
     
        stimax=media(xsFit);
        stimax1.push_back(stimax);
        stimam=media(ysFit);
        stimamax.push_back(stimam);
        //------------------------------------QUA PROVAVO A FARE IL FIT PARABOLICO---------------------------//
        // pic<<stimax<<" "<<stimam<<endl;
        // cout<<"Il picco vale: x-> "<<stimax<<" mentre la y-> "<<stimam<<endl;

        // if (xsFit.size() >= 3) {
        //     ParabolaFitResult fit = fitParabola(xsFit, ysFit);
        //     Punto vert = trovaVertice(fit);
        //     // salva vertice su file e stampalo a schermo
        //     pic << vert.x << " " << vert.y << "\n";
        //     cout << "  parabola vertice: t=" << vert.x
        //          << ", p=" << vert.y << "\n";
        // }
        idx += double(len);
    }
    // pic.close();

    //Una volta che abbiamo un vector di massimi ora prendiamo no covarianti
    vector<double> tnocorr;
    double tmedio;
    // for(auto c:stimax1){
    //     cout<<c<<" "<<endl;
    // 

    for(int k=0; k<stimax1.size()-1;k++){
        tnocorr.push_back((stimax1.at(k+1)-stimax1.at(k)));
        k=k+1;
    }

    tmedio=media(tnocorr);
    cout<< "Numero periodi: " << tnocorr.size() << endl;
    cout<<"Il periodo medio stimato complessivo (considerando periodi separati): "<<tmedio <<" del file "<<nomefile<<endl;
    cout<<"Omega per il file "<<nomefile<<" vale: "<<(2*M_PI/tmedio)<<" [unita di misura] "<<endl;
    cout<<"---------------------------------------------------------------------------------------------------------"<<endl;
    

   

    return 0;
}

//questa funzione va a caricare i vettori dei valori che stanno sopra lo zero, 
int estraiPositivi(const vector<double>& p,const vector<double>& t,int startIndex,vector<double>& out_p,vector<double>& out_t){
    out_p.clear();
    out_t.clear();
    int i = startIndex;
    // Se il primo non è positivo, niente da fare
    if (i >= p.size() || p.at(i) <= 0.0) return 0;
    int count = 0;
    while (i < p.size() && p.at(i) > 0.0) {
        out_p.push_back(p.at(i));
        out_t.push_back(t.at(i));
        ++i;  //mettere prima o dopo il simbolo ++ qua non cambia nulla, ma a livello di ragionamento prima incrementa e poi usa la variabile 
        ++count;
    }
    return count;
}

double media(vector<double> v){
    double somma=0.0;
    for(int i=0;i<v.size();i++){
        somma+=v.at(i);
    }

    return somma/v.size();
}

//-----------------------fit parabolico
// ParabolaFitResult fitParabola(const vector<double>& x, const vector<double>& y) {
//     size_t n = x.size();
//     double Sx=0, Sx2=0, Sx3=0, Sx4=0;
//     double Sy=0, Sxy=0, Sx2y=0;
//     for (size_t i = 0; i < n; ++i) {
//         double xi = x[i], yi = y[i], xi2 = xi*xi;
//         Sx   += xi;
//         Sx2  += xi2;
//         Sx3  += xi2*xi;
//         Sx4  += xi2*xi2;
//         Sy   += yi;
//         Sxy  += xi*yi;
//         Sx2y += xi2*yi;
//     }
//     double D  = Sx4*(Sx2*n - Sx*Sx)
//               - Sx3*(Sx3*n - Sx*Sx2)
//               + Sx2*(Sx3*Sx - Sx2*Sx2);
//     double Da = Sx2y*(Sx2*n - Sx*Sx)
//               - Sxy*(Sx3*n - Sx*Sx2)
//               + Sy*(Sx3*Sx - Sx2*Sx2);
//     double Db = Sx4*(Sxy*n - Sx*Sy)
//               - Sx3*(Sx2y*n - Sx*Sy)
//               + Sx2*(Sx2y*Sx - Sxy*Sx2);
//     double Dc = Sx4*(Sx2*Sy - Sx*Sx2y)
//               - Sx3*(Sx3*Sy - Sx*Sxy)
//               + Sx2*(Sx3*Sxy - Sx2*Sx2y);

//     return { Da/D, Db/D, Dc/D };
// }

// Punto trovaVertice(const ParabolaFitResult& parabola) {
//     double xv = -parabola.b / (2.0 * parabola.a);
//     double yv = parabola.a*xv*xv + parabola.b*xv + parabola.c;
//     return { xv, yv };
// }

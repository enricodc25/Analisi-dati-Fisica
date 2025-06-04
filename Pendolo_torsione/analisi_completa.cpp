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

//Prototipi
Punto trovaVertice(const ParabolaFitResult& parabola);
int estraiPositivi(const vector<double>& p,const vector<double>& t,int startIndex,vector<double>& out_p,vector<double>& out_t);
int estraiNegativi(const vector<double>& p,const vector<double>& t,int startIndex,vector<double>& out_p,vector<double>& out_t);
double media(vector<double> v);
double deviazioneStandardCampionaria(vector<double>& dati);
void interpolazionesemplice(vector<double> x, vector<double> y, double& a);

int main() {
    //  vector<string> nomifile={"f0.9.txt","f0.905.txt","f0.910.txt","f0.915.txt","f0.925.txt", "f0.930.txt", "f0.935.txt", "f0.940.txt", "f0.945.txt", "f0.964.txt","f0.965.txt","f0.970.txt"
    //  ,"f0.975.txt", "f0.980.txt", "f0.06470.txt", "f0.9200.txt", "f0.9625.txt", "f0.9900.txt","f0.96450.txt", "f0.96460.txt"};
         vector<double> t, f, p, a, fa;
         int scelta;
         string nomefile;
         cout<<"Inserisci il nome del file: ";
         cin>>nomefile;
         cout<<"Il file è di smorzamento? ( inserire 0, altrimenti un altro numero) ";
         cin>>scelta;
        //  double M_PI= 3.14159226535897932384; 
    
//CARICAMENTO DATI DEL FILE
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
        idx += double(len);    
    }
    cout<<"_____________________INIZIANO I MINIMI____________________"<<endl<<endl;

    //||||||||||________________________PARTE DI CALCOLO DEI MINIMI_________________________|||||||||//
    vector<double> xsMinFit, ysMinFit;
    vector<double> stimin1, stiminval;  // stimin1 = ascisse dei minimi, stiminval = valori dei minimi
    idx = 0.0;
    idxraccolta = 0.0;
    double minval;
    bloccoP.clear();
    bloccoT.clear();
    //controllo della soglia sensata per i minimi
    vector<double> mintotali;
    while (idxraccolta < p.size()) {
        int len = estraiNegativi(p, t, int(idxraccolta), bloccoP, bloccoT);
        if (len == 0) {
            idxraccolta += 1.0;
            continue;
        }
        double minVal = *min_element(bloccoP.begin(), bloccoP.end());
        mintotali.push_back(minVal);
        idxraccolta += double(len);
    }
 
    double mediaMintotali = media(mintotali);
    double sogliadeimin = mediaMintotali / 3.0;  // Soglia basata sulle medie dei minimi

    while (idx < p.size()) {
        int lenNeg = estraiNegativi(p, t, int(idx), bloccoP, bloccoT);
        if (lenNeg == 0) {
            idx += 1.0;
            continue;
        }
        double minVal = *min_element(bloccoP.begin(), bloccoP.end());
    //qua nel caso di smorzamento non deve esserci il controllo della soglia?
        
            if (minVal > sogliadeimin) {  
            idx += double(lenNeg);
            continue;
            }
        
        
    
        
        double sogliaMin = 0.70 * minVal;
        cout << "Bloc Neg #" << idx << ": min=" << minVal
             << " -> sogliaMin=" << sogliaMin << "\n";

        xsMinFit.clear();
        ysMinFit.clear();
        for (int i = 0; i < bloccoP.size(); ++i) {
            if (bloccoP.at(i) <= sogliaMin) {
                ysMinFit.push_back(bloccoP.at(i));
                xsMinFit.push_back(bloccoT.at(i));
            }
        }
        double stimin = media(xsMinFit);
        double stimin_y = media(ysMinFit);
        stimin1.push_back(stimin);
        stiminval.push_back(stimin_y);
        idx += double(lenNeg);
    }


    //Una volta che abbiamo un vector di massimi ora prendiamo no covarianti
    vector<double> tnocorr;
    double tmedio;
    double maxMedio=media(stimamax);
    for(int k=0; k<stimax1.size()-1;k++){
        tnocorr.push_back((stimax1.at(k+1)-stimax1.at(k)));
        k=k+1;
    }
    tmedio=media(tnocorr);

    //Ora ci occupiamo dei minimi
    double minMedio=media(stiminval);
    double dev_std = deviazioneStandardCampionaria(tnocorr);
    int n = tnocorr.size();
    double erroreomegaf  = ((2*M_PI)/ (tmedio*tmedio))*(dev_std/sqrt(n));

    //stima theta particolare: per la sua stima non posso usare la media brutale di tutti i massimi e di tutti i minimi, ma per evitare la correlazione
    //devo fare in modo di "separare la loro stima"
    double thetaMedio;
    vector<double> thetaparticolari;
    double theta;
    for(int h=0; h<stiminval.size()&&h<stimamax.size(); h++){
        theta=(stimamax.at(h)-stiminval.at(h))/2; //metto il meno perchè so che i dati contenuti qua sono negativi
        thetaparticolari.push_back(theta);
        h++;
    }
    thetaMedio=media(thetaparticolari);
    double devtheta=(deviazioneStandardCampionaria(thetaparticolari))/(sqrt(thetaparticolari.size()));

    cout<< "Numero periodi: " << tnocorr.size() << endl;
    cout<<"Il picco massimo medio vale: "<<maxMedio<<endl;
    cout<<"Il picco minimo medio vale: "<<minMedio<<endl;
    cout<<"Il periodo medio stimato complessivo(considerando periodi separati): "<<tmedio <<" del file "<<nomefile<<endl;
    cout<<"Omega per il file "<<nomefile<<" vale: "<<(2*M_PI/tmedio)<<" [unita di misura] "<<endl;
    cout << "Deviazione standard campionaria del periodo: " << dev_std << endl;
    cout << "Errore sulla media del periodo: " << dev_std/sqrt(n) << endl;
    cout << "Errore omega f è " << erroreomegaf <<endl;
    cout<<"La theta particolare (semiampiezza) con dati scorrelati vale: "<<thetaMedio*2*M_PI<<endl;
    cout<<"La deviazione standard del theta particolare medio (semiampiezza): "<<devtheta*2*M_PI<<endl;
    cout<<"---------------------------------------------------------------------------------------------------------"<<endl;


    //Analisi smorzamento
    

    if(scelta==0){
        ofstream stampamax("maxgamma.txt");
        vector<double> appog;
        vector<double> appog2;
        //nel vector thetaparticolari ho i tetha da convertire
        cout<<"-----------Vari theta particolari con rispettivo periodo: "<<endl;
        cout<<"____PER I MASSIMI: ___"<<endl;
        for(int i=0;i<stimax1.size();i++){
            cout<<"-> tempo: "<<stimax1.at(i)<<" theta: "<<stimamax.at(i)<<endl;
            //aggiorno i dati in logaritmi naturali
            appog.push_back(log(stimamax.at(i)*2*M_PI));
        }       
        // for(int i=0;i<appog.size();i++){
        //     stampamax<<stimax1.at(i)<<" "<<appog.at(i)<<endl;
        // }
        cout<<"_____PER I MINIMI: ____"<<endl;
        for(int i=0;i<stimin1.size();i++){
            cout<<"-> tempo: "<<stimin1.at(i)<<" theta: "<<stiminval.at(i)<<endl;
             appog2.push_back(log(fabs(stiminval.at(i) * 2*M_PI)));
        }
    double gammamax,gammamin;
    interpolazionesemplice(stimax1,appog,gammamax);
    interpolazionesemplice(stimin1,appog2,gammamin);

    cout<<"Il valore di gamma per i massimi: "<<gammamax<<endl;
    cout<<"Il valore di gamma per i minimi: "<<gammamin<<endl;
    //calcolo l'omega di risonanza atteso a livello teorico
    //bisogna prima esaminare se è legittimo fare la media delle gamma, soprattuto per quanto riguarda l'analisi poi degli errori
    double mediagamma=(gammamax+gammamin)/2;
    double omegarteorico=sqrt((2*M_PI/tmedio)*(2*M_PI/tmedio)-mediagamma*mediagamma);
    //bisogna confrontare questo valore con l'oomega di risonanza ricavato attraveros il fit parabolico con l'altro programma ( test gaussiano) 
    cout<<"La omega TEORICA di risonanza è pari a: "<<omegarteorico<<endl;

    }

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
    return count; //in questo modo mi tengo in memoria quanti dati ho passato in rassegna
}

int estraiNegativi(const vector<double>& p,const vector<double>& t,int startIndex,vector<double>& out_p,vector<double>& out_t){
    out_p.clear();
    out_t.clear();
    int i = startIndex;
    if (i >= p.size() || p.at(i) >= 0.0) return 0;
    int count = 0;
    while (i < p.size() && p.at(i) < 0.0) {
        out_p.push_back(p.at(i));
        out_t.push_back(t.at(i));
        ++i;  
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
//calcolo l'errore sulla media del periodo medio
// Deviazione standard campionaria
double deviazioneStandardCampionaria(vector<double>& dati) {
    int n = dati.size();
    double somma = 0.0;
    for (double valore : dati) {
        somma += valore;
    }
    double media = somma / n;
    double sommaQuadrati = 0.0;

    for (double valore : dati) {
        sommaQuadrati += (valore - (somma/n)) * (valore - (somma/n));
    }
    return sqrt(sommaQuadrati / (n - 1));
}
//INTERPOLAZIONE SEMPLICE CON INCERTEZZE SULLE Y TUTTE UGUALI
void interpolazionesemplice(vector<double> x, vector<double> y, double& b){  
    double sx = 0.0;   
    double sy = 0.0;  
    double sxx = 0.0;  
    double sxy = 0.0; 
    double n=x.size(); 
    double sigmay=deviazioneStandardCampionaria(y);

	//calcolo somme intermedie
    for (int i = 0; i < n; i++) {
        sx  += x.at(i); 
        sy  += y.at(i);
        sxx += x.at(i) * x.at(i);
        sxy += x.at(i) * y.at(i);
    }

    // Calcolo del DELTA
    double delta = n*(sxx) - (sx * sx);

    // Stima dei parametri a e b
    double a = ( (sxx * sy) - (sx * sxy) ) / delta;
    b = ( (n* sxy) - (sx * sy ) ) / delta;

    //Calcolo delle incertezze relative ai parametri   
    double sigma_a = sigmay*sqrt(sxx/ delta);
    double sigma_b = sigmay*sqrt(n/ delta);

    cout << "=== RISULTATI DEL FIT LINEARE SEMPLICE per una retta y=a+bx ===" << endl;
    cout << "a = " << a << endl;
    cout << "b = " << b << endl;
}

//ora inserisco l'interpolazione lineare che però verrà modificata ( perchè devo fare il logaritmo delle theta)
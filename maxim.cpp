#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>
#include <limits>
using namespace std;

// Functia pentru citirea problemei de programare liniara
void citire_problema(
    vector<vector<double>> & A,
    vector<double> & b,
    vector<double> & c,
    vector<int> & tipuri_restrictii
) {
    int n, m;

    cout << "Introduceti numarul de restrictii: ";
    cin >> n;

    cout << "Introduceti numarul de variabile: ";
    cin >> m;

    // Redimensionarea matricei A si vectorului b
    A.resize(n, vector<double>(m, 0));
    b.resize(n, 0);
    tipuri_restrictii.resize(n, 0);

    // Citirea coeficientilor functiei obiectiv
    cout << "Introduceti coeficientii functiei obiectiv (max z = c1*x1 + c2*x2 + ... + cm*xm):\n";
    c.resize(m, 0);
    for (int j = 0; j < m; j++) {
        cout << "c" << j + 1 << " = ";
        cin >> c[j];
    }

    // Citirea restrictiilor
    cout << "Introduceti coeficientii restrictiilor:\n";
    for (int i = 0; i < n; i++) {
        cout << "Restrictia " << i + 1 << ":\n";
        for (int j = 0; j < m; j++) {
            cout << "a" << i + 1 << j + 1 << " = ";
            cin >> A[i][j];
        }

        cout << "Tipul restrictiei (1: <=, 2: =, 3: >=): ";
        cin >> tipuri_restrictii[i];

        cout << "Termenul liber b" << i + 1 << " = ";
        cin >> b[i];
    }
}

// Functie pentru afisarea problemei in forma matriceala
void afisare_forma_matriceala(const vector<vector<double>>& A, const vector<double>& b, const vector<double>& c) {
    int n = A.size(); // Numarul de restrictii
    int m = A[0].size(); // Numarul de variabile

    cout << "\nProblema in forma matriceala:\n";

    // Afisam functia obiectiv
    cout << "Maximizeaza z = ";
    for (int j = 0; j < m; j++) {
        cout << c[j] << " * x" << j + 1;
        if (j < m - 1) cout << " + ";
    }
    cout << "\n\n";

    // Afisam matricea coeficientilor A si vectorul termenilor liberi b
    cout << "Matricea coeficientilor A si vectorul termenilor liberi b:\n";
    for (int i = 0; i < n; i++) {
        cout << "| ";
        for (int j = 0; j < m; j++) {
            cout << setw(8) << fixed << setprecision(4) << A[i][j] << " ";
        }
        cout << " | " << setw(8) << fixed << setprecision(4) << b[i] << " |\n";
    }
    cout << "\n";
}

// Functia pentru aducerea problemei la forma standard
void aducere_la_forma_standard(
    vector<vector<double>> & A,
    vector<double> & b,
    vector<double> & c,
    vector<int> & tipuri_restrictii,
    vector<vector<double>> & A_standard,
    vector<double> & b_standard,
    vector<double> & c_standard,
    vector<int> & baza_initiala
) {
    int n = A.size(); // numarul de restrictii
    int m = A[0].size(); // numarul de variabile initiale

    // Numarul de variabile de compensare (slack) si artificiale
    int nr_slack = 0;
    int nr_artificiale = 0;

    for (int i = 0; i < n; i++) {
        if (tipuri_restrictii[i] == 1) { // <=
            nr_slack++;
        } else if (tipuri_restrictii[i] == 3) { // >=
            nr_slack++;
            nr_artificiale++;
        }
    }

    // Dimensiunea noii matrice A_standard
    int m_nou = m + nr_slack + nr_artificiale;
    A_standard.resize(n, vector<double>(m_nou, 0));
    c_standard.resize(m_nou, 0);
    b_standard = b;
    baza_initiala.resize(n, -1);

    // Copiem coeficientii variabilelor initiale
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            A_standard[i][j] = A[i][j];
        }
    }

    // Copiem coeficientii functiei obiectiv pentru variabilele initiale
    for (int j = 0; j < m; j++) {
        c_standard[j] = c[j];
    }

    // Adaugam variabilele de compensare (slack) si artificiale
    int index_slack = m;
    int index_artificial = m + nr_slack;

    for (int i = 0; i < n; i++) {
        if (tipuri_restrictii[i] == 1) { // <=
            // Adaugam variabila de compensare (slack)
            A_standard[i][index_slack] = 1.0;
            baza_initiala[i] = index_slack;
            index_slack++;
        } else if (tipuri_restrictii[i] == 3) { // >=
            // Adaugam variabila de compensare (slack) cu semn negativ
            A_standard[i][index_slack] = -1.0;
            index_slack++;

            // Adaugam variabila artificiala (penalizare)
            A_standard[i][index_artificial] = 1.0;
            c_standard[index_artificial] = -1000.0; // Penalizare mare pentru variabilele artificiale
            baza_initiala[i] = index_artificial;
            index_artificial++;
        }
    }
}

// Functia pentru crearea tabelului simplex initial
vector<vector<double>> creare_tabel_simplex(
    const vector<vector<double>> & A_standard,
    const vector<double> & b_standard,
    const vector<double> & c_standard,
    const vector<int> & baza
) {
    int n = A_standard.size(); // numarul de restrictii
    int m = A_standard[0].size(); // numarul total de variabile

    // Tabelul simplex are dimensiunea (n+1) x (m+1)
    vector<vector<double>> tabel(n + 1, vector<double>(m + 1, 0));

    // Completam tabelul cu coeficientii restrictiilor si termenii liberi
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            tabel[i][j] = A_standard[i][j];
        }
        tabel[i][m] = b_standard[i]; // termenii liberi pe ultima coloana
    }

    // Calculam valorile Z_j si c_j - Z_j pentru fiecare coloana
    for (int j = 0; j < m; j++) {
        double z_j = 0;
        for (int i = 0; i < n; i++) {
            int var_baza = baza[i];
            z_j += c_standard[var_baza] * tabel[i][j];
        }
        tabel[n][j] = c_standard[j] - z_j; // c_j - z_j
    }

    // Calculam valoarea functiei obiectiv Z
    double Z = 0;
    for (int i = 0; i < n; i++) {
        int var_baza = baza[i];
        Z += c_standard[var_baza] * tabel[i][m];
    }
    tabel[n][m] = Z; // valoarea functiei obiectiv

    return tabel;
}

void afisare_tabel_simplex(
    const vector<vector<double>> & tabel,
    const vector<int> & baza,
    int iteratie,
    const vector<double> & c,
    pair<int, int> pivot = {-1, -1}
) {
    int n = tabel.size() - 1; // numarul de restrictii
    int m = tabel[0].size() - 1; // numarul de variabile + 1

    cout << "\n\nTabelul simplex - Iteratia " << iteratie << ":\n";

    // Afisam antetul tabelului
    cout << setw(8) << "C_B" << " | " << setw(8) << "Baza" << " | " << setw(8) << "Solutia" << " | ";
    for (int j = 0; j < m; j++) {
        cout << setw(8) << "a" << j + 1 << " | ";
    }
    cout << endl;

    cout << string(10 + (m + 3) * 12, '-') << endl;

    // Afisam liniile restrictiilor
    for (int i = 0; i < n; i++) {
        cout << setw(8) << c[baza[i]] << " | " << setw(8) << "a" << baza[i] + 1 << " |  ";
        cout << setw(8) << fixed << setprecision(2) << tabel[i][m] << " | ";
        for (int j = 0; j < m; j++) {
            cout << setw(8) << fixed << setprecision(2) << tabel[i][j] << " |  ";
        }
        cout << endl;
    }

    cout << string(10 + (m + 3) * 12, '-') << endl;

    // Afisam linia functiei obiectiv (c_j - z_j)
    cout << setw(8) << " " << " | " << setw(8) << " " << " | " << setw(8) << " " << " | ";
    for (int j = 0; j < m; j++) {
        cout << setw(8) << fixed << setprecision(2) << tabel[n][j] << "  |  ";
    }

    // Afisam valoarea functiei obiectiv Z
    double Z = tabel[n][m];
    cout <<"z_"<<iteratie<<"="<<setw(8) << fixed << setprecision(2) << -Z << endl;

    // Afisam pivotul
    if (pivot.first != -1 && pivot.second != -1) {
        cout << "\nPivot: Linia " << pivot.first + 1 << ", Coloana " << pivot.second + 1 << endl;
    }
}

// Functia pentru o iteratie a algoritmului simplex
// Functia pentru o iteratie a algoritmului simplex
pair<bool, pair<int, int>> iteratie_simplex(
    vector<vector<double>> & tabel,
    vector<int> & baza
) {
    int n = tabel.size() - 1; // numarul de restrictii
    int m = tabel[0].size() - 1; // numarul de variabile + 1

    // Determinam variabila care intra in baza (coloana pivot)
    int coloana_pivot = -1;
    double valoare_maxima = 0;

    for (int j = 0; j < m; j++) {
        if (tabel[n][j] > valoare_maxima) {
            coloana_pivot = j;
            valoare_maxima = tabel[n][j];
        }
    }

    // Daca nu exista variabila care sa intre in baza, solutia este optima
    if (valoare_maxima <= 0) {
        return {false, {-1, -1}}; // am terminat, solutia este optima
    }

    // Determinam variabila care iese din baza (linia pivot)
    int linia_pivot = -1;
    double raport_minim = numeric_limits<double>::max();

    for (int i = 0; i < n; i++) {
        if (tabel[i][coloana_pivot] > 0) {
            double raport = tabel[i][m] / tabel[i][coloana_pivot];
            if (raport < raport_minim) {
                raport_minim = raport;
                linia_pivot = i;
            }
        }
    }

    // Daca nu exista variabila care sa iasa din baza, problema este nemarginita
    if (linia_pivot == -1) {
        cout << "Problema este nemarginita!" << endl;
        return {false, {-1, -1}};
    }

    // Actualizam baza
    baza[linia_pivot] = coloana_pivot;

    // Calculam valoarea pivotului
    double pivot_valoare = tabel[linia_pivot][coloana_pivot];

    // Impartim linia pivot prin valoarea pivotului pentru a face pivotul egal cu 1
    for (int j = 0; j <= m; j++) {
        tabel[linia_pivot][j] /= pivot_valoare;
    }

    // Actualizam restul tabelului
    for (int i = 0; i <= n; i++) {
        if (i != linia_pivot) {
            double factor = tabel[i][coloana_pivot];
            for (int j = 0; j <= m; j++) {
                tabel[i][j] -= factor * tabel[linia_pivot][j];
            }
        }
    }

    return {true, {linia_pivot, coloana_pivot}}; // Continuam iteratia
}

















// Functia pentru rezolvarea problemei folosind metoda simplex


// Functia pentru rezolvarea problemei folosind metoda simplex
void rezolvare_simplex(
    const vector<vector<double>> & A_standard,
    const vector<double> & b_standard,
    const vector<double> & c_standard,
    vector<int> & baza
) {
    // Cream tabelul simplex initial
    vector<vector<double>> tabel = creare_tabel_simplex(A_standard, b_standard, c_standard, baza);

    // Afisam tabelul initial
    afisare_tabel_simplex(tabel, baza, 0, c_standard);

    int iteratie = 1;
    while (true) {
        auto rezultat = iteratie_simplex(tabel, baza);
        bool continua = rezultat.first;
        auto pivot = rezultat.second;

        if (!continua) break;

        // Afisam tabelul dupa fiecare iteratie
        afisare_tabel_simplex(tabel, baza, iteratie, c_standard, pivot);
        iteratie++;
    }

    // Verificam daca variabilele artificiale sunt in baza cu valori nenule
    int n = tabel.size() - 1; // numarul de restrictii
    int m = A_standard[0].size(); // numarul de variabile (fara variabilele artificiale)
    bool solutie_fezabila = true;

    for (int i = 0; i < n; i++) {
        if (baza[i] >= m && abs(tabel[i][tabel[0].size() - 1]) > 1e-6) {
            solutie_fezabila = false;
            break;
        }
    }

    if (!solutie_fezabila) {
        cout << "\nProblema nu are solutie fezabila!\n";
        return;
    }

    // Afisam solutia finala
    cout << "\nSolutia optima a fost gasita!\n";

    // Extragem si afisam valoarea optima a functiei obiectiv (fara influenta variabilelor artificiale)
    double Z = 0;
    for (int i = 0; i < n; i++) {
        int var_baza = baza[i];
        if (var_baza < m) { // ignoram variabilele artificiale
            Z += c_standard[var_baza] * tabel[i][tabel[0].size() - 1];
        }
    }

    cout << "Valoarea optima a functiei obiectiv: Z = " << Z << "\n";

    // Extragem si afisam valorile variabilelor din solutia optima
    cout << "Valorile variabilelor din solutia optima:\n";
    for (int j = 0; j < m; j++) {
        bool este_in_baza = false;
        int index_in_baza = -1;

        // Verificam daca variabila este in baza
        for (int i = 0; i < n; i++) {
            if (baza[i] == j) {
                este_in_baza = true;
                index_in_baza = i;
                break;
            }
        }

        // Daca variabila este in baza, afisam valoarea ei
        if (este_in_baza) {
            cout << "x" << j + 1 << " = " << tabel[index_in_baza][tabel[0].size() - 1] << "\n";
        } else {
            // Daca variabila nu este in baza, valoarea ei este 0
            cout << "x" << j + 1 << " = 0\n";
        }
    }

    // Afisam valorile variabilelor de compensare (slack) si artificiale
    for (int j = m; j < tabel[0].size() - 1; j++) {
        bool este_in_baza = false;
        int index_in_baza = -1;

        // Verificam daca variabila este in baza
        for (int i = 0; i < n; i++) {
            if (baza[i] == j) {
                este_in_baza = true;
                index_in_baza = i;
                break;
            }
        }

        // Daca variabila este in baza, afisam valoarea ei
        if (este_in_baza) {
            cout << "x" << j + 1 << " = " << tabel[index_in_baza][tabel[0].size() - 1] << "\n";
        } else {
            // Daca variabila nu este in baza, valoarea ei este 0
            cout << "x" << j + 1 << " = 0\n";
        }
    }
}        














// Functia principala
int main() {
    vector<vector<double>> A; // matricea coeficientilor restrictiilor
    vector<double> b, c; // vectorii de termeni liberi si coeficienti
    vector<int> tipuri_restrictii;

    // Citirea datelor problemei
    citire_problema(A, b, c, tipuri_restrictii);

    // Afisarea problemei in forma matriceala
    afisare_forma_matriceala(A, b, c);

    // Aducem problema la forma standard
    vector<vector<double>> A_standard;
    vector<double> b_standard, c_standard;
    vector<int> baza_initiala;

    aducere_la_forma_standard(A, b, c, tipuri_restrictii, A_standard, b_standard, c_standard, baza_initiala);

    // Rezolvam problema folosind metoda simplex
    rezolvare_simplex(A_standard, b_standard, c_standard, baza_initiala);

    return 0;
}
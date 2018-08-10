#ifndef __PROGTEST__
#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <cctype>
#include <climits>
#include <cmath>
#include <cfloat>
#include <cassert>
#include <unistd.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#if defined ( __cplusplus ) && __cplusplus > 199711L	/* C++ 11 */
#include <memory>
#endif /* C++ 11 */
using namespace std;
#endif /* __PROGTEST__ */

class CPolynomial
{
public:
  // default constructor
    CPolynomial();
  // operator +
  // operator -
  // operator * (polynomial, double)
    CPolynomial operator+(const CPolynomial& other) const;
    CPolynomial operator-(const CPolynomial& other) const;
    CPolynomial operator*(const CPolynomial& other) const;
    CPolynomial operator*(double a) const; //const double ne, protoze se predava hodnotou ne referenci (ale slo by to, jen se to u prim. typu nevyplati)
  // operator ==
  // operator !=
    bool operator==(const CPolynomial& other) const;
    bool operator!=(const CPolynomial& other) const;
  // operator []
    double& operator[] (size_t deg); //pro nastaveni a[1] = neco (pokud neni, vratit nulu)
    double operator[] (size_t deg) const; //nemeni polynom, pro neco = a[1]
  // operator ()
    double operator()(double x) const;
  // Degree (), returns unsigned value
    unsigned int Degree()const;
    friend ostream& operator<<(ostream& os, const CPolynomial&); //fce, ne metoda, proto prijima i ostreampr
private:
CPolynomial add(const CPolynomial& bigger, const CPolynomial& smaller) const;
  vector<double> koef; 
};

CPolynomial::CPolynomial(){
    koef.resize(1);
}

double& CPolynomial::operator[] (size_t idx){
if ( idx < 0)
    throw "index out of range";
//#dfzoloyoloswagkoef.insert(koef.begin() + idx, 0);
if (idx >= koef.size()){
    //cout << "index " << idx << " je vetsi nez aktualni velikost" << endl;
    koef.resize(idx+1);
}

  return koef[idx];
}

double CPolynomial::operator[] (size_t idx)const{
if ( idx < 0)
    throw "index out of range";
if (idx >= koef.size())
    return 0;


  return koef[idx];
}

unsigned int CPolynomial::Degree()const{
    unsigned int deg = koef.size() - 1;
    if(deg > 0){
        while(koef[deg] == 0 && deg > 0)
           deg--;
    }
    return deg;
}

double CPolynomial::operator()(double x) const{
    double sum = koef[0];
    for(size_t i = 1; i < koef.size(); i ++){
        double exp = x;
        for(size_t j = 1; j < i; j++){
            exp*=x;
            //cout << "x na " << i << " = " << exp << endl;
        }
        sum+= exp * koef[i];
       // cout << "mezisoucet: " << sum << endl;
    }
    return sum;
}

ostream& operator<<(ostream& os, const CPolynomial& poly){
    if(poly.koef[poly.Degree()] < 0)
        os << "- ";
    
    
    for(int i = poly.Degree(); i > 0 ; i--){
        if(poly.koef[i] != 0){
        
            if(abs(poly.koef[i]) != 1)
                os << abs(poly.koef[i]) <<"*";
            os << "x^" << i;
        }
        
        if(poly.koef[i-1] > 0)
            os << " + ";
        else if(poly.koef[i-1] < 0)
            os << " - ";
    
    
    }
    //absolutni clen, nulu vypiseme pouze, pokud se jedna o prazdny polynom
    if(poly.koef[0] != 0 || poly.Degree() == 0)
        os << abs(poly.koef[0]);
    
    return os;
}

CPolynomial CPolynomial::operator*(double a) const{
    CPolynomial res;
    if( a == 0 )
        return res;
    for(size_t i = 0; i < koef.size(); i++){
        res[i] = koef[i] * a;
    }
    return res;
}

CPolynomial CPolynomial::operator+(const CPolynomial& other) const{
    CPolynomial res;
    
    if(this->koef.size() > other.koef.size())
        res = add(*this, other);
    else
        res = add(other, *this);
    
    
    return res;
}

CPolynomial CPolynomial::operator-(const CPolynomial& other) const{
    CPolynomial res;
    
    if(this->koef.size() > other.koef.size()){
        for(size_t i = 0; i < this->koef.size(); i++ ){
            if(i < other.koef.size())
                res[i] = this->koef[i] - other.koef[i];
            else
                res[i] = this->koef[i];
        }
    }else{
        for(size_t i = 0; i < other.koef.size(); i++ ){
            if(i < this->koef.size())
                res[i] = this->koef[i] - other.koef[i];
            else
                res[i] = -(other.koef[i]);
        }
    }
    
    
    return res;
}

CPolynomial CPolynomial::operator*(const CPolynomial& other) const{
    CPolynomial res;
    res.koef.resize(this->Degree() + other.Degree() + 1);
    
    for(size_t i = 0; i < other.koef.size(); i++ )
        for(size_t j = 0; j < this->koef.size(); j++ )
            res.koef[i+j] += this->koef[j] * other.koef[i];
    
    return res;
}

CPolynomial CPolynomial::add(const CPolynomial& bigger, const CPolynomial& smaller) const{
    CPolynomial res;
    for (size_t i = 0; i < bigger.koef.size(); i++)
        if(i < smaller.koef.size())
            res[i] = smaller.koef[i] + bigger.koef[i];
        else
            res[i] = bigger.koef[i];
    return res;
    
}

bool CPolynomial::operator!=(const CPolynomial& other) const{
    if( *this == other )
        return false;    
    //cout << "nerovnaji se" << endl;
    return true;
}

bool CPolynomial::operator==(const CPolynomial& other) const{
    if (this->Degree() != other.Degree()){
        //cout << "Polynomy nestejneho stupne" << endl;
        return false;
    }
    
    for(unsigned int i = 0; i <= Degree(); i++){
        if( koef[i] != other.koef[i] )
                return false;
    }
    return true;
}

#ifndef __PROGTEST__
bool smallDiff (double a, double b){
    return abs(a -b) < 0.000001; 
}

bool dumpMatch (const CPolynomial & x, const vector <double>&ref){
    //hranaty zavorky + degree
    for(unsigned int i = 0; i <= x.Degree(); i++){
        if( x[i] != ref[i] )
                return false;
    }
    return true;
}

int
main (void){
    
  CPolynomial a, b, c;
  
  ostringstream out;

  a[0] = -10;
  a[1] = 3.5;
  //double m = a[1]; //jine uziti operatoru!
  a[3] = 1;
  //a[4] = 0;
  
  /*cout << m << endl;
  
  for (int i = 0; i < 5; i++)
          cout << a[i] << endl;
  
  cout << a.Degree() << endl;*/
  
 // cout << "Hodnota v bode 2: " << a(2) << endl;
 
  assert (smallDiff (a (2), 5));
  out.str ("");
    
  out << a;
   //cout << out.str() << endl;
  assert (out.str () == "x^3 + 3.5*x^1 - 10");
  a = a * -2;  
  //out1 << a;
   //cout << out1.str() << endl;
  
  assert (a.Degree () == 3 && dumpMatch (a, vector < double >
					 {
					 20.0, -7.0, -0.0, -2.0}));

  out.str ("");
  out << a;
  assert (out.str () == "- 2*x^3 - 7*x^1 + 20");
  out.str ("");
  out << b;
  
  assert (out.str () == "0");
  b[5] = -1;
  out.str ("");
  out << b;
   //cout << "os: " << out.str() << endl;
  assert (out.str () == "- x^5");
  c = a + b;
  assert (c.Degree () == 5 && dumpMatch (c, vector < double >
					 {
					 20.0, -7.0, 0.0, -2.0, 0.0, -1.0}));

  out.str ("");
  out << c;
  assert (out.str () == "- x^5 - 2*x^3 - 7*x^1 + 20");
  c = a - b;
  assert (c.Degree () == 5 && dumpMatch (c, vector < double >
					 {
					 20.0, -7.0, -0.0, -2.0, -0.0, 1.0}));

  out.str ("");
  out << c;
  assert (out.str () == "x^5 - 2*x^3 - 7*x^1 + 20");
  c = a * b;
  assert (c.Degree () == 8 && dumpMatch (c, vector < double >
					 {
					 -0.0, -0.0, 0.0, -0.0, 0.0, -20.0,
					 7.0, -0.0, 2.0}));

  out.str ("");
  out << c;
  assert (out.str () == "2*x^8 + 7*x^6 - 20*x^5");
  assert (a != b);
  b[5] = 0;
  assert (!(a == b));
  a = a * 0;
  assert (a.Degree () == 0 && dumpMatch (a, vector < double >
					 {
					 0.0}));

  assert (a == b);

  return 0;
}
#endif /* __PROGTEST__ */
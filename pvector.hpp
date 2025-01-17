#ifndef _PVECTOR_
#define _PVECTOR_
// EXERCISE create an efficient vector class in C++ with methods and overloaded operators which you can find in the TODO LIST.
// TODO LIST:
// * = higher priority, X = done
// [X] constructor (with overloading->initializer_list)
// [X] destructor (empty)
// [X] show()
// [X] overloading of =  <----
// [X] overloading of + operator for addition of vectors
// [X] overloading of - operator for substraction of vectors
// [X] overloading of += and -= operators (both with vectors and scalars)
// [X] sum()
// [X] get( ) (get i-th element)
// [X] set( , ) (set i-th element)
//
// [X] overloading of [...] operator (to use vector as C arrays, e.g. v[0]=1
// [X] comma initialization: overloading of << and , operators
// [X] "scalar time vector" through friend function
// [X] overloading of * operator for "vector times scalar" product
// [X] overloading of operator == (check whether two vectors are equal)
// [X] overloading of *= and /= (with scalars)
// [X] overloading of * operator for scalar product
// [X] overloading of ^ operator for cross product 
// [X] norm(), calculate the modulus of a vector
//
// [X] rint 
// [X] overloading of << for output (friend function)
// [X] mulcw, multiplication element wise
// [X] divcw, division element wise
// [ ] random(L), random vector inside a box of length L
// [ ] random_orient() random unit on unit sphere (marsaglia) 
#include<initializer_list> // inizializzazione vettori
#include<iostream> // input/output
#include<cmath>
#include<string>
//
#include <fstream>
#include <sstream>
#include <vector>
#include <numeric> // libreria standard che fornisce una serie di funzioni per eseguire operazioni numeriche su sequenze di dati
#include <iomanip>
#include <string>
// STRATEGIE PER RENDERE GENERICA LA CLASSE
//#define NT 3
//constexpr int NT=3;

//typedef float ntype;
//using ntype=double;

// int main(void)
// {
//    pvector vec({1,2,3});
//
// }
template <typename ntype, int NT=3>
class pvector
{
  ntype v[NT]; // dato privato
public:
  pvector() // constructor
    {
      int i;
      for(i=0; i < NT; i++)
        {
          v[i]=0;
        }
    }
  ~pvector() // destructor
    {
    }

  ntype norm() const
    {
      return sqrt((*this)*(*this));
    }

  // pvector A, B; 
  // A.mulcw(B);
  pvector mulcw(const pvector& vec)
    {
      pvector vt;
      for (int i=0; i < NT; i++)
        {
          vt[i] = v[i]*vec[i];
        }
      return vt;
    }
 
  pvector divcw(const pvector& vec)
    {
      pvector vt;
      for (int i=0; i < NT; i++)
        {
          vt[i] = v[i]/vec[i];
        }
      return vt;
    }
  
  pvector operator^(const pvector& vec) const
    {
      if (NT==3)
        {
          pvector vt;
          vt[0] = v[1]*vec[2]-v[2]*vec[1];
          vt[1] = v[2]*vec[0]-v[0]*vec[2];
          vt[2] = v[0]*vec[1]-v[1]*vec[0];
          return vt;
        }
      else
        {
          std::cout << "Cross product not defined\n";
          exit(1);
        }
    }
  // (std::cout << v) << "\n";
  friend std::ostream& operator<<(std::ostream& os, const pvector& vec)
    {
      os << "(";
      for (int i=0; i < NT; i++)
        {
          os << vec[i];
          if (i < NT-1)
           os << ","; 
        }
      os << ")";
      return os;
    }
  // multiply by scalar and assign result to vector
  pvector& operator *=(ntype s) 
    {
      for (int i=0; i < NT; i++)
        v[i] *= s;
      return (*this);
    }
 
  // divide by scalar and assign result to vector
  pvector& operator /=(ntype s)
    {
      for (int i=0; i < NT; i++)
        v[i] /= s;
      return (*this);
    }

  // prodotto scalare
  ntype operator*(const pvector& vec) const
    {
      ntype sp=0;
      for (int i=0; i < NT; i++)
        sp += v[i]*vec[i];
      return sp;
    }
 
  // prodotto vettore per scalare
  pvector operator*(ntype s) const
    {
      pvector vt;
      for (int i=0; i < NT; i++)
        vt[i] = v[i]*s;
      return vt;
    }

  // comma initialization
  // A << 1.0,2.0,3.0;
  int counter;
  pvector& operator<<(ntype val)
    {
      counter=0;
      v[0] = val;
      return (*this);
    } 

  pvector& operator,(ntype val)
    {
      counter++;
      if( counter >= NT)
        {
          std::cout << "Troppi elementi forniti\n";
          exit(1);
        }

      v[counter] = val;

      return (*this);
    }
  // =================================================

  // overloading of () operator to get/set i-th element
  // v(i) = 2.0;
  // cout << v(i);
  ntype &operator()(int i)
    {
      return v[i];
    }

  ntype operator()(int i) const
    {
      return v[i];
    }

  // ===================================================
  // overloading dell'operatore A[B] A chiama l'operatore parentesi quadre con argomento B,
  // dove B è un intero
  // A[2] = 1.0;
  ntype& operator[](int i)
    {
      return v[i];
    }

  ntype operator[](int i) const
    {
      return v[i];
    }

  // prodotto scalare vettore 
  // s*vec 
  friend pvector operator*(ntype s, const pvector& vec)
    {
      pvector vt;
      for (int i=0; i < NT; i++)
        {
          vt.v[i] = s*vec.v[i];
        }
      return vt;
    }

  pvector(std::initializer_list<ntype> list)  // overloading constructor
    {
      int c=0;
      for (auto el: list) // itero su tutti gli elementi della lista "list" e
                          // quindi "el" sarà
                          // ognuno di questi elementi 
        {
          if (c < NT)
            {
              v[c] = el;
            }
          c++;
        }
      for (;c < NT; c++) // gli elementi sono in numero di NT allora inizializzo 0
        {
          v[c]=0.0;
        }
    }

  // comparison between two vectors
  bool operator==(const pvector& V2) const
    {
      for (auto i=0; i < NT; i++)
        {
          if (V2.v[i]!=v[i])
            {
              return 0;
            }
        }
      return 1;
    }

  // (A=B)=C;
  pvector& operator=(const pvector& v2)
    {
      for (int i=0; i < NT; i++)
        {
          (*this).v[i] = v2.v[i];
        }
      return(*this);
    }

  // sum 
  // A.sum(B) questo è equivalente a scrivere A+B
  pvector sum(const pvector& v2) const
    {
      return (*this)+v2;
    }

  // get
  ntype get(int i) const
    {
      return (*this).v[i];
    }

  // set
  ntype set(int i, ntype val) 
    {
      return (*this).v[i]=val;
    }

  // operator+: (*this)+v2
  pvector operator+(const pvector& v2) const
    {
      pvector vs;
      for (int i=0; i < NT; i++)
        {
          vs.v[i] = (*this).v[i] + v2.v[i];
        }
      return vs;
    } 

  // operator-
  pvector operator-(const pvector& v2) const
    {
      pvector vs;
      for (int i=0; i < NT; i++)
        {
          vs.v[i] = (*this).v[i] - v2.v[i];
        }
      return vs;
    } 

  pvector& operator+=(const pvector& v2)
    {
      for (int i=0; i < NT; i++)
        {
          v[i] += v2.v[i];
        }
      return (*this);
    } 

  pvector& operator-=(const pvector& v2)
    {
      for (int i=0; i < NT; i++)
        {
          v[i] -= v2.v[i];
        }
      return (*this);
    } 

  void show(std::string s="") const // argomento di default 
    {
      std::cout << s << "(";
      for (int i=0; i < NT; i++)
        {
          std::cout << v[i];
          if (i < NT-1) 
            std::cout << ",";
        }
      std::cout << ")";
    }
};

template<typename ntype, int NT>
pvector<ntype,NT> rint(const pvector<ntype,NT>& vec)
{
  pvector<ntype,NT> vt;

  for (int i=0; i < NT; i++)
    {
      vt[i] = rint(vec[i]);
    }
  return vt;
}
template <typename ntype>
using pvec3=pvector<ntype,3>;
#endif

// array.hh: generic, templated array class
// see array.cc for a simple example of its use
// resizes itself dynamically

// Paul Heckbert	10 Feb 1999

#ifndef ARRAY_H
#define ARRAY_H

#include <iostream>
#include <assert.h>

template <class Type> class Array;

template <class Type>
class Array {				// Array<Type>: dynamic array
    Type *data;				// an array of length cap

  public:
    int cap;				// capacity of array
    int n;				// number currently in array
					// must satisfy 0<=n<=cap

    Array<Type>();			// create an empty array
    
    int num() const			// number of elements in array
	{return n;}
    void resize(int num);		// grow or shrink array to length num

    Type &operator[](int i);		// array indexing, for assignment

    Type &operator[](int i) const;	// array indexing, for const objects
};

template <class Type>
Array<Type>::Array() {
    cap = 8;	// default array capacity
    data = new Type[cap];
    assert(data);
    n = 0;
};

template <class Type>
void Array<Type>::resize(int num) {
    // cout << "....resize from " << cap << " to " << num << endl;
    Type *ndata = new Type[num];
    assert(ndata);
    int i;
    // copy all n, if n<num, else only copy the first num elements
    if (num<n) n = num;
    for (i=0; i<n; i++)
	ndata[i] = data[i];
    delete[] data;
    data = ndata;
    cap = num;
}

template <class Type>
Type &Array<Type>::operator[](int i) {
    // this version is used for assignment, e.g. "Array<int> a; a[i] = 7;"
    if (i<0) {
	cerr << "Array<Type>::operator[](" << i << ")" << endl;
	exit(1);
    }
    if (i>=n) {
	if (i>=cap) {
	    int newcap = cap*2;
	    if ((i+1)*3/2>newcap) newcap = (i+1)*3/2;
	    resize(newcap);
	}
	n = i+1;
    }
    return data[i];
}

template <class Type>
Type &Array<Type>::operator[](int i) const {
    // this version needed for reading elements of const arrays,
    // e.g. "printfirst(const Array<int> &a) {cout << a[0];}"
    if (i<0 || i>=n) {
	cerr << "Array<Type>::operator[](" << i << ") but n=" << n << endl;
	exit(1);
    }
    return data[i];
}

// to support the following, the caller must supply
// ostream &operator<<(ostream &s, const Type &x) {	// print item

template <class Type>
ostream &operator<<(ostream &s, const Array<Type> &a) {	// print array
    s << "array n=" << a.num() << " cap=" << a.cap << ": ";
    int i;
    for (i=0; i<a.num(); i++)
	s << "[" << i << "]=" << a[i] << ", ";
    s << "end" << endl;
    return s;
}

#endif

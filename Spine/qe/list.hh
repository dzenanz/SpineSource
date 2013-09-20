// list.hh: generic, templated singly-linked list and stack class
// see list.cc for a simple example of its use

// Paul Heckbert	28 Sep 1996 - version 1: List_item derived from Type
//					this turned out to be awkward
//			9 Feb 1999 - version 2: List_item now has "Type *obj"

#ifndef LIST_H
#define LIST_H

#include <iostream>
#include <assert.h>

template <class Type> class List;

//////////////////////////////////////////////////////////////////
// A List_item<Type> is an item in a list.
//
// Internally, it consists of a next pointer and a "Type *" pointer.

template <class Type>
class List_item {			// List_item<Type>: item in linked list
    friend class List<Type>;
    List_item<Type> *nextp;		// next in list

  public:
    Type *obj;				// object in list
    List_item<Type>(Type *x);		// construct a list item

    const List_item<Type> *next() const // next in list (NULL if end)
	{return nextp;}
    List_item<Type> *&next()            // next in list (NULL if end)
	{return nextp;}
    // if you understand this, you get a PhD
};

template <class Type>
List_item<Type>::List_item(Type *x) {	// constructor
    // make an item structure that points to x, but no more
    obj = x;
    nextp = 0;				// null for now
}

//////////////////////////////////////////////////////////////////
// A List<Type> is internally a header structure pointing to a linked list
// of List_item's of type Type.

template <class Type>
class List {			   // List<Type>: linked list or stack of items
  public:
  //?? these two should be protected, but compiler barfs
    List_item<Type> *firstp;		// first item in list
    List_item<Type> *lastp;		// last item in list

    List()                              // create an empty list
	{init();}
    void init()                         // create an empty list
	{firstp = lastp = 0;}

    // When a list is copied with operator=, we want header to be copied
    // but items to be shared (not copied).
    // To get that behavior, we do not define "List &List(List &)" or
    // "operator=(List &)" but let them default to memberwise-initialization.
    // If you want a completely new copy of a list, use copy() below.

    ~List()				// destructor
        {free_items();}
        // {cout << "~List()" << endl; free_items();}

    void free_items();			// something like a destructor
    // empties out list by deleting the List_items only

    void free_all();			// something like a destructor
    // empties out list, deleting List_items and the objects themselves

    void copy(const List<Type> &);	// copies list & items (not objects)
    void prepend(Type *);		// prepend item at beginning of list
    void append(Type *);		// append item to end of list
    void concat(List<Type> *);		// concatenate list on end of this one

    List_item<Type> *first() const	// return first in list
	{return firstp;}

    List_item<Type> *last() const	// return last in list
	{return lastp;}

    Type *remove(List_item<Type> *);	// remove one item from list, return it

    // stack stuff
    void push(Type *x)			// push item on stack
	{prepend(x);}
    void dup();				// duplicate top item
    Type *pop();			// pop item off stack
    int length(void) const;		// return length of list
};

// to support the following, the caller must supply
// ostream &operator<<(ostream &s, const Type &x) {	// print item

template <class Type>
ostream &operator<<(ostream &s, const List<Type> &l) {	// print list
    s << "list (length " << l.length() << "): ";
    List_item<Type> *p;
    for (p=l.first(); p; p=p->next())
	s << *p->obj << ", ";
    s << "end" << endl;
    return s;
}

// the routine above must be in the header file or g++ will barf when
// somebody tries to write operator<< for an instantiated list type,
// e.g. "operator<<(ostream &, const Prim &)" trying to print a List<Light>

//////////////////////////////////////////////////////////////////
// the following must be in the header file for g++.
// need not be in the header file for SGI CC (it's a smarter compiler)

template <class Type>
Type *List<Type>::pop() {			// pop item off stack
    List_item<Type> *p = firstp;
    if (!p) return 0;				// stack empty, can't pop
    firstp = p->next();
    if (!firstp) lastp = 0;
    Type *x = p->obj;
    delete p;
    return x;
}

template <class Type>
int List<Type>::length() const {                // return length of list
    List_item<Type> *p;
    int count;
    for (count=0, p=first(); p; p=p->next(), count++);
    return count;
}

template <class Type>
void List<Type>::free_items() {
    return;  // seg fault from this  -andrewb ???
    Type *x;
    do
	x = pop();
    while (x);
}

//////////////////////////////////////////////////////////////////
// Note: be careful calling free_all if you've done any
// list copying or object sharing because there will be multiple pointers
// to the same List_item or object (aliasing).

template <class Type>
void List<Type>::free_all() {
    Type *x;
    do {
	x = pop();
	if (x) delete x;
    } while (x);
}

template <class Type>
void List<Type>::copy(const List<Type> &l) {
    // copy list l into *this, copying each list item
    // first we empty out the existing list
    free_items();
    List_item<Type> *p;
    for (p=l.firstp; p; p=p->nextp)
	this->append(p->obj);
}

template <class Type>
void List<Type>::prepend(Type *x) {		// prepend item at beg. of list
    List_item<Type> *p = new List_item<Type>(x);
    assert(p);
    p->nextp = firstp;
    firstp = p;
    if (!lastp) lastp = p;
}

template <class Type>
void List<Type>::append(Type *x) {		// append item to end of list
    List_item<Type> *p = new List_item<Type>(x);
    assert(p);
    if (firstp)
	lastp = lastp->nextp = p;
    else
	firstp = lastp = p;
    p->nextp = 0;
}

template <class Type>
void List<Type>::concat(List<Type> *a) {	// concat list a on end of this
    assert(a);
    if (firstp)
	lastp->nextp = a->firstp;
    else
	firstp = a->firstp;
    lastp = a->lastp;
    a->firstp = 0;	// so delete doesn't wipe out the items
    delete a;		// delete the list a itself
}

template <class Type>
Type *List<Type>::remove(List_item<Type> *q) {	// remove item q from list
    assert(q);
    List_item<Type> *p;
    Type *t = q->obj;
    if (first()==q) {
	firstp = q->nextp;
	if (!firstp) lastp = 0;	// list is now empty
	delete q;
	return t;
    }
    for (p=first(); p; p=p->next())
	if (p->next()==q) {
	    // we found the predecessor to q, now splice the list
	    p->nextp = q->nextp;
	    if (!p->nextp) {assert(lastp==q); lastp = p;} //??nuke assert
	    delete q;
	    return t;
	}
    cerr << "can't remove nonexistent item" << endl;
    return t;
}

template <class Type>
void List<Type>::dup() {			// duplicate top item
    if (!firstp) {
	cerr << "ERROR: attempt to dup() first item in empty list" << endl;
	exit(1);
    }
    List_item<Type> *p = new List_item<Type>(firstp->obj);
    assert(p);
    push(p);
}

#endif

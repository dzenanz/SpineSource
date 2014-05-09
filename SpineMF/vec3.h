/*
	math3d++ - A C++ 3d math library
	Copyright (c) 2004, Trenkwalder Markus
	All rights reserved. 
	
	Redistribution and use in source and binary forms, with or without
	modification, are permitted provided that the following conditions
	are met:
	
	- Redistributions of source code must retain the above copyright notice,
	  this list of conditions and the following disclaimer.
	  
	- Redistributions in binary form must reproduce the above copyright
	  notice, this list of conditions and the following disclaimer in the
	  documentation and/or other materials provided with the distribution.
	  
	- Neither the name of the math3d++'s copyright owner nor the names
	  of its contributors may be used to endorse or promote products
	  derived from this software without specific prior written permission.
	  
	THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
	"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
	LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
	A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
	OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
	SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
	TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
	PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
	LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING 
	NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
	SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
	
	Contact info:
	email: trenki2@gmx.net web: trenki.50free.org
*/

#ifndef MATH3DPP_VEC3_H
#define MATH3DPP_VEC3_H

#include <cmath>
typedef float Real;

/// 3d vector class
class vec3 {
private:
	Real data[3];
public:
	vec3();
	vec3(Real x, Real y, Real z);
    vec3(Real *ptr);
	
	vec3& operator+= (const vec3& v);
	vec3& operator-= (const vec3& v);
	
	vec3& operator*= (Real s);	
	vec3& operator/= (Real s);
	
	operator Real* ();
	operator const Real* () const;

    Real length() const;
    Real length2() const;
    const Real* ptr() const;

    Real normalize();
};

// constructor

/// initialization to (0,0,0)
inline vec3::vec3() {
	for ( int i = 0; i < 3; ++i )
		data[i] = 0;
}

/// initialization to (x,y,z)
inline vec3::vec3(Real x, Real y, Real z) {
	data[0] = x;
	data[1] = y;
	data[2] = z;
}

inline vec3::vec3(Real *ptr) {
	for ( int i = 0; i < 3; ++i )
		data[i] = ptr[i];
}

// Vector arithmetic operators

inline vec3& vec3::operator+= (const vec3& v) {
	for ( int i = 0; i < 3; ++i )
		data[i] += v.data[i];
		
	return *this;
}


inline vec3& vec3::operator-= (const vec3& v) {
	for ( int i = 0; i < 3; ++i )
		data[i] -= v.data[i];
		
	return *this;
}


inline vec3& vec3::operator*= (Real s) {
	for ( int i = 0; i < 3; ++i )
		data[i] *= s;
		
	return *this;
}


inline vec3& vec3::operator/= (Real s) {
	for ( int i = 0; i < 3; ++i )
		data[i] /= s;
		
	return *this;
}

/// Conversion operator for subscripting.
/// You can write v[2] = value; or value = v[2].
inline vec3::operator Real* () {
	return data;
}

/// Conversion operator for subscripting.
/// You can write v[2] = value; or value = v[2].
inline const Real* vec3::ptr() const {
	return data;
}

/// Conversion operator for subscripting (const version).
/// You can write value = v[2].
inline vec3::operator const Real* () const {
	return data;
}

/******************************************************************************/
// non-member operators

/// @relates vec3
/// Vector addition
inline vec3 operator+ (const vec3& a, const vec3& b) {
	vec3 r = a;
	return r += b;
}

/// @relates vec3
/// Vector subtraction
inline vec3 operator- (const vec3& a, const vec3& b) {
	vec3 r = a;
	return r -= b;
}

/// @relates vec3
/// Vector multiplication by a scalar
inline vec3 operator* (Real s, const vec3& v) {
	vec3 r = v;
	return r *= s;
}

/// @relates vec3
/// Vector multiplication by a scalar
inline vec3 operator* (const vec3& v, Real s) {
	vec3 r = v;
	return r *= s;
}


/// @relates vec3
/// Divides each vector component by a scalar value
inline vec3 operator/ (const vec3& v, Real s) {
	vec3 r = v;
	return r /= s;
}

/// @relates vec3
/// Unary minus
inline vec3 operator- (const vec3& v) {
	vec3 r = v;
	for ( int i = 0; i < 3; ++i )
		r[i] = -r[i];

	return r;
}

/// @relates vec3
/// Comparison operators
inline bool operator== (const vec3& a, const vec3& b) {
	for ( int i = 0; i < 3; ++i ) {
	#ifndef USE_FEQUAL_COMPARE
		if ( a[i] != b[i] ) return false;
	#else
		if ( !fequal(a[i], b[i]) ) return false;
	#endif
	}
	
	return true;
}

/// @relates vec3
/// Vector comparison
inline bool operator!= (const vec3& a, const vec3& b) {
	return !operator==(a,b);
}

/// @relates vec3
/// Vector dot product
inline Real operator* (const vec3& a, const vec3& b) {
	Real dotprod = 0.0;
	for ( int i = 0; i < 3; ++i )
		dotprod += a[i] * b[i];
	
	return dotprod;
}

/// @relates vec3
/// Vector dot product
inline Real dot(const vec3& a, const vec3& b) {
	return a * b;
}

/// @relates vec3
/// Returns the vector length
inline Real vec3::length() const {
	return std::sqrt(length2());
}

/// @relates vec3
/// Returns the length^2
inline Real vec3::length2() const {
	Real l = 0.0;
	for ( int i = 0; i < 3; ++i )
		l += data[i] * data[i];
	
	return l;
}

/// @relates vec3
/// Vector cross product only for 3d vectors
inline vec3 cross(const vec3& a, const vec3& b) {
	vec3 r;
	r[0] = a[1] * b[2] - b[1] * a[2];
	r[1] = a[2] * b[0] - b[2] * a[0];
	r[2] = a[0] * b[1] - b[0] * a[1];
	return r;
}

/// @relates vec3
/// Vector cross product only for 3d vectors
inline vec3 operator^ (const vec3& a, const vec3& b) {
	return cross(a,b);
}

/// @relates vec3
/// Returns the vector scaled to unit length
inline vec3 normalized(const vec3& v)  {
	Real l = v.length();
	
	if ( l != 0.0 ) {
		vec3 result(v[0] / l, v[1] / l, v[2] / l);
		return result;
	}
	
	return v;
}

/// @relates vec3
/// Normalizes the vector and returns previous length
inline Real vec3::normalize()  {
    Real l = length();

	if ( l != 0.0 ) {
		data[0]/=l;
        data[1]/=l;
        data[2]/=l;
	}
	
	return l;
}

/// @relates vec3
/// Scales the passed vector to unit length and returns a reference to it
inline vec3& normalize(vec3& v) {
	Real l = v.length();
	
	if ( l != 0.0 ) {
		v[0] /= l;
		v[1] /= l;
		v[2] /= l;
	}
	
	return v;
}

#endif

#include "IntVector.hpp"
#include "Msg.hpp"



IntVector::IntVector(void) {

    v = NULL;
    allocateVector(1);
    v[0] = 0;
}

IntVector::IntVector(size_t n) {

    v = NULL;
    allocateVector(n);
    for (int i=0; i<dim; i++)
        v[i] = 0;
}

IntVector::IntVector(size_t n, int val) {

    v = NULL;
    allocateVector(n);
    for (int i=0; i<dim; i++)
        v[i] = val;
}

IntVector::IntVector(const IntVector& vec) {

    v = NULL;
    allocateVector(vec.dim);
    for (int i=0; i<dim; i++)
        v[i] = vec.v[i];
}

IntVector::IntVector(const std::vector<int>& vec) {

    v = NULL;
    allocateVector(vec.size());
    for (int i=0; i<dim; i++)
        v[i] = vec[i];
}

IntVector::~IntVector(void) {

    freeVector();
}

void IntVector::add(IntVector& iVec) {

	if (iVec.dim != dim)
        Msg::error("Vector sizes don't match in add1.");
    for (int i=0; i<dim; i++)
        v[i] += iVec.v[i];
}

void IntVector::add(std::vector<int>& iVec) {

	if (iVec.size() != dim)
        Msg::error("Vector sizes don't match in add2.");
    for (int i=0; i<dim; i++)
        v[i] += iVec[i];
}
    
void IntVector::add(int* iVec, int n) {

	if (n != dim)
        Msg::error("Vector sizes don't match in add2.");
    for (int i=0; i<dim; i++)
        v[i] += iVec[i];
}

void IntVector::addMultiple(IntVector& iVec, int iMultiple) {

	if (iVec.dim != dim)
        Msg::error("Vector sizes don't match in addMultiple1.");
    for (int i=0; i<dim; i++)
        v[i] += iVec.v[i] * iMultiple;
}

void IntVector::addMultiple(std::vector<int>& iVec, int iMultiple) {

	if (iVec.size() != dim)
        Msg::error("Vector sizes don't match in addMultiple2.");
    for (int i=0; i<dim; i++)
        v[i] += iVec[i] * iMultiple;
}

void IntVector::addMultiple(int* iVec, int n, int iMultiple) {

	if (n != dim)
        Msg::error("Vector sizes don't match in addMultiple2.");
    for (int i=0; i<dim; i++)
        v[i] += iVec[i] * iMultiple;
}

void IntVector::allocateVector(size_t n) {

    if (v != NULL)
        freeVector();
    if (n > 0)
        {
        dim = n;
        v = new int[n];
        }
    else
        Msg::error("Cannot initialize a vector of zero length");
}

void IntVector::clean(void) {

    for (int i=0; i<dim; i++)
        v[i] = 0;
}

void IntVector::clean(void) const {

    for (int i=0; i<dim; i++)
        v[i] = 0;
}

void IntVector::freeVector(void) {

    if (v != NULL)
        delete [] v;
    v = NULL;
    dim = 0;
}

int IntVector::innerProduct(IntVector& iVec) {

	if (iVec.dim != dim)
        Msg::error("Vector sizes don't match in innerProduct1.");
    int sum = 0;
    for (int i=0; i<dim; i++)
        sum += iVec.v[i] * v[i];
    return sum;
}

int IntVector::innerProduct(std::vector<int>& iVec) {

	if (iVec.size() != dim)
        Msg::error("Vector sizes don't match in innerProduct2 (" + std::to_string(iVec.size()) + " vs " + std::to_string(dim) + ").");
    int sum = 0;
    for (int i=0; i<dim; i++)
        sum += iVec[i] * v[i];
    return sum;
}

int IntVector::innerProduct(int* iVec, int n) {

	if (n != dim)
        Msg::error("Vector sizes don't match in innerProduct2 (" + std::to_string(n) + " vs " + std::to_string(dim) + ").");
    int sum = 0;
    for (int i=0; i<dim; i++)
        sum += iVec[i] * v[i];
    return sum;
}

void IntVector::subtract(IntVector& iVec) {

	if (iVec.dim != dim)
        Msg::error("Vector sizes don't match in subtract1.");
    for (int i=0; i<dim; i++)
        v[i] -= iVec.v[i];
}
    
void IntVector::subtract(std::vector<int>& iVec) {

	if (iVec.size() != dim)
        Msg::error("Vector sizes don't match in subtract2.");
    for (int i=0; i<dim; i++)
        v[i] -= iVec[i];
}

void IntVector::subtract(int* iVec, int n) {

	if (n != dim)
        Msg::error("Vector sizes don't match in subtract2.");
    for (int i=0; i<dim; i++)
        v[i] -= iVec[i];
}

bool IntVector::zeroEntry(void) {

    for (size_t i=0; i<dim; i++)
        {
        if (v[i] == 0)
            return true;
        }
    return false;
}

IntVector& IntVector::operator=(const int& a) {

    for (size_t i=0; i<dim; i++)
        v[i] = a;
    return *this;
}

IntVector& IntVector::operator=(const IntVector& vec) {

    if (dim != vec.dim)
        {
        Msg::warning("Assigning vector of size " + std::to_string(dim) + " to be of size " + std::to_string(vec.dim));
        freeVector();
        allocateVector(vec.dim);
        }
    for (size_t i=0; i<dim; i++)
        v[i] = vec.v[i];
    return *this;
}

bool IntVector::operator==(const IntVector& vec) const {

    if (dim != vec.dim)
        return false;
    for (size_t i=0; i<dim; i++)
        {
        if (v[i] != vec.v[i])
            return false;
        }
    return true;
}

std::ostream& operator<<(std::ostream& s, const IntVector& vec) {

	size_t n = vec.size();
	s << "[" << n << "] (";
	for (size_t i=0; i<n; i++)
        {
		s << vec[i];
		if (i != n-1)
			s << ",";
        }
	s << ")";
	return s;
}

IntVector operator+(const IntVector& vecA, const IntVector& vecB) {

	size_t n = vecA.size();
	if (vecB.size() != n )
		return IntVector();
    IntVector vecC(n);
    for (size_t i=0; i<n; i++)
        vecC[i] = vecA[i] + vecB[i];
    return vecC;
}

IntVector operator-(const IntVector& vecA, const IntVector& vecB) {

	size_t n = vecA.size();
	if (vecB.size() != n )
		return IntVector();
    IntVector vecC(n);
    for (size_t i=0; i<n; i++)
        vecC[i] = vecA[i] - vecB[i];
    return vecC;
}

IntVector operator/(const IntVector& vecA, const IntVector& vecB) {

	size_t n = vecA.size();
	if (vecB.size() != n )
		return IntVector();
    IntVector vecC(n);
    for (size_t i=0; i<n; i++)
        vecC[i] = vecA[i] / vecB[i];
    return vecC;
}

IntVector& operator+=(IntVector& vecA, const IntVector& vecB) {

	size_t n = vecA.size();
	if (vecB.size() == n)
        {
		for (size_t i=0; i<n; i++)
			vecA[i] += vecB[i];
        }
	return vecA;
}

IntVector& operator-=(IntVector& vecA, const IntVector& vecB) {

	size_t n = vecA.size();
	if (vecB.size() == n)
        {
		for (size_t i=0; i<n; i++)
			vecA[i] -= vecB[i];
        }
	return vecA;
}

IntVector operator*(const IntVector& vecA, const IntVector& vecB) {

	size_t n = vecA.size();
	if (vecB.size() != n )
		return IntVector();
    IntVector vecC(n);
    for (size_t i=0; i<n; i++)
        vecC[i] = vecA[i] * vecB[i];
    return vecC;
}

IntVector& operator*=(IntVector& vecA, const IntVector& vecB) {

	size_t n = vecA.size();
	if (vecB.size() == n)
        {
		for (size_t i=0; i<n; i++)
			vecA[i] *= vecB[i];
        }
	return vecA;
}

IntVector& operator*=(IntVector& vecA, const int b) {

    for (size_t i=0; i<vecA.size(); i++)
        vecA[i] *= b;
	return vecA;
}

IntVector& operator/=(IntVector& vecA, const IntVector& vecB) {

	size_t n = vecA.size();
	if (vecB.size() == n)
        {
		for (size_t i=0; i<n; i++)
			vecA[i] /= vecB[i];
        }
	return vecA;
}

IntVector& operator/=(IntVector& vecA, const int b) {

    for (size_t i=0; i<vecA.size(); i++)
        vecA[i] /= b;
	return vecA;
}

bool operator!=(const IntVector& vecA, const IntVector& vecB) {

	if (vecA == vecB)
		return false;
	else
		return true;
}

IntVector operator*(const IntVector& vecA, const int b) {

    size_t n = vecA.size();
    IntVector vecC(n);
    for (size_t i=0; i<n; i++)
        vecC[i] = vecA[i] * b;
    return vecC;
}

IntVector operator*(const int b, const IntVector& vecA) {

    size_t n = vecA.size();
    IntVector vecC(n);
    for (size_t i=0; i<n; i++)
        vecC[i] = vecA[i] * b;
    return vecC;
}



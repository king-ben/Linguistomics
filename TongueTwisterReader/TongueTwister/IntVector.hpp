#ifndef IntVector_hpp
#define IntVector_hpp

#include <iostream>
#include <vector>



class IntVector {

    public:
                                IntVector(void);
                                IntVector(size_t n);
                                IntVector(size_t n, int val);
                                IntVector(const IntVector& vec);
                                IntVector(const std::vector<int>& vec);
                               ~IntVector(void);
        IntVector&              operator=(const int& a);
        IntVector&              operator=(const IntVector& vec);
        bool                    operator==(const IntVector& vec) const;
        int&                    operator[](size_t i) { return v[i]; }
        const int&              operator[](size_t i) const { return v[i]; };
        friend std::ostream&    operator<<(std::ostream& s, const IntVector& vec);
        friend IntVector        operator+(const IntVector& vecA, const IntVector& vecB);
        friend IntVector        operator-(const IntVector& vecA, const IntVector& vecB);
        friend IntVector        operator*(const IntVector& vecA, const IntVector& vecB);
        friend IntVector        operator*(const IntVector& vecA, const int b);
        friend IntVector        operator*(const int b, const IntVector& vecA);
        friend IntVector        operator/(const IntVector& vecA, const IntVector& vecB);
        friend IntVector&       operator+=(IntVector& vecA, const IntVector& vecB);
        friend IntVector&       operator-=(IntVector& vecA, const IntVector& vecB);
        friend IntVector&       operator*=(IntVector& vecA, const IntVector& vecB);
        friend IntVector&       operator*=(IntVector& vecA, const int b);
        friend IntVector&       operator/=(IntVector& vecA, const IntVector& vecB);
        friend IntVector&       operator/=(IntVector& vecA, const int b);
        friend bool             operator!=(const IntVector& vecA, const IntVector& vecB);
        void                    add(IntVector& iVec);
        void                    add(std::vector<int>& iVec);
        void                    add(int* iVec, int n);
        void                    addMultiple(IntVector& iVec, int iMultiple);
        void                    addMultiple(std::vector<int>& iVec, int iMultiple);
        void                    addMultiple(int* iVec, int n, int iMultiple);
        void                    clean(void);
        void                    clean(void) const;
        int                     innerProduct(IntVector& iVec);
        int                     innerProduct(std::vector<int>& iVec);
        int                     innerProduct(int* iVec, int n);
        size_t                  size(void) { return dim; }
        size_t                  size(void) const { return dim; }
        void                    subtract(IntVector& iVec);
        void                    subtract(std::vector<int>& iVec);
        void                    subtract(int* iVec, int n);
        bool                    zeroEntry(void);
            
    private:
        void                    allocateVector(size_t n);
        void                    freeVector(void);
        size_t                  dim;
        int*                    v;
};

struct CompIntVector {

    bool operator()(const IntVector& v1, const IntVector& v2) const {
        
        if (v1.size() > v2.size())
            return true;
        else if ( v1.size() == v2.size())
            {
            for (int i=0; i<v1.size(); i++)
                {
                if (v1[i] > v2[i])
                    return true;
                else if (v1[i] < v2[i])
                    return false;
                }
            }
        return false;
        }

    bool operator()(const IntVector* v1, const IntVector* v2) const {

        if (v1->size() > v2->size())
            return true;
        else if ( v1->size() == v2->size())
            {
            for (int i=0; i<v1->size(); i++)
                {
                if ((*v1)[i] > (*v2)[i])
                    return true;
                else if ((*v1)[i] < (*v2)[i])
                    return false;
                }
            }
        return false;
        }
};

#endif

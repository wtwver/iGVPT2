#ifndef MATRIX
#define MATRIX
#include <vector>
#include <iostream>
#include <cmath>
#include <complex>

using namespace std;
typedef complex<double> cdouble;

template <class T>
class Matrix
{
	int numberOfRows;
	int numberOfColumns;
	vector<vector<T> > elements;
	public :
	Matrix();
	Matrix(int rows, int cols);
	Matrix(const Matrix<T> &M);
	Matrix<T> &operator= (const Matrix<T> &M);
	vector<T> &operator[] (int row);
	int getNumberOfRows() const;
	int getNumberOfColumns() const;
	vector<T> getRow(int num) const;
	vector<T> getColumn(int num) const;
	Matrix<T> transpose();
	Matrix<T> scale(T factor);
	void clear();
	Matrix<T> operator+ (Matrix<T> &M);
	Matrix<T> operator- (Matrix<T> &M);
	Matrix<T> operator* (Matrix<T> &M);
	vector<T> operator* (const vector<T> &v);
	Matrix<T> multDiag(const vector<T> &v);
	T cofactorSign(int row, int col);
	Matrix<T> minminor(int row, int col);
	T determinant();
	// Get/Set the value of one element.
	T operator()(const int, const int) const;
	// Get the value of one element.
	T& operator()(const int, const int);
	bool checkBounds(int a, int b) const;
	void add(int a, int b);
	std::vector<T> slice(int row);
	void resize(int n);
};
#endif

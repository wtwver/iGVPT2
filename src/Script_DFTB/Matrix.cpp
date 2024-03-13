#include <vector>
#include <iostream>
#include <cmath>
#include <complex>
#include "Matrix.h"

using namespace std;

template <class T>
Matrix<T>::Matrix() : numberOfRows(0), numberOfColumns(0) {}

template <class T>
Matrix<T>::Matrix(int rows, int cols)
{
	numberOfRows = rows;
  	numberOfColumns = cols;
  
	for (int i=0; i < numberOfRows; i++)
    	{
      		vector<T> col(numberOfColumns);
       		elements.push_back(col);
	}
}
template <class T>
Matrix<T>::Matrix(const Matrix<T> &M)
{
	numberOfRows = M.numberOfRows;
	numberOfColumns = M.numberOfColumns;
	elements = M.elements;
}

template <class T>
Matrix<T> &Matrix<T>::operator= (const Matrix<T> &M)
{
	numberOfRows = M.numberOfRows;
	numberOfColumns = M.numberOfColumns;
	elements = M.elements;
	return *this;
}

template <class T>
vector<T> &Matrix<T>::operator[] (int row)
{
	return elements[row];
}

template <class T>
int Matrix<T>::getNumberOfRows() const
{
	return numberOfRows;
}

template <class T>
int Matrix<T>::getNumberOfColumns() const
{
	return numberOfColumns;
}

template <class T>
vector<T> Matrix<T>::getRow(int num) const
{
	vector<T> Result(numberOfColumns);
	for(int x=0; x<numberOfColumns; x++)
		Result[x] = (elements[num][x]); 
	return Result;
}

template <class T>
vector<T> Matrix<T>::getColumn(int num) const
{
	vector<T> Result(numberOfRows);
	for(int y=0; y<numberOfRows; y++) Result[y] = (elements[y][num]);
	return Result;
}

template <class T>
Matrix<T> Matrix<T>::transpose()
{
	Matrix<T> M_Trans(numberOfColumns, numberOfRows);
	for(int row=0; row<M_Trans.numberOfRows; row++)
	for(int mrow=0; mrow<numberOfRows; mrow++)
		M_Trans[row][mrow] = elements[mrow][row];
	return M_Trans;
}

template <class T>
Matrix<T> Matrix<T>::scale(T factor)
{
	Matrix Result(numberOfRows, numberOfColumns);
	for(int row=0; row<Result.numberOfRows; row++)
	for(int col=0; col<Result.numberOfColumns; col++)
		Result[row][col] = factor*(elements[row][col]);
	return Result;
} 

template <class T>
Matrix<T> Matrix<T>::operator+ (Matrix<T> &M)
{
	Matrix<T> M2 = *this;
	Matrix<T> Result(M2.numberOfRows, M2.numberOfColumns);
  
	for(int row=0; row<M2.numberOfRows; row++)
	for(int col=0; col<M2.numberOfColumns; col++)
		Result[row][col] = M2[row][col] + M[row][col];

	return Result;
}
  
  
template <class T>
Matrix<T> Matrix<T>::operator- (Matrix<T> &M)
{
	Matrix<T> M2 = *this;
	Matrix<T> Result(M2.numberOfRows, M2.numberOfColumns);
  
	for(int row=0; row<M2.numberOfRows; row++)
		for(int col=0; col<M2.numberOfColumns; col++)
      			Result[row][col] = M2[row][col] - M[row][col];

	return Result;
}

template <class T>
Matrix<T> Matrix<T>::operator* (Matrix<T> &M)
{
	Matrix<T> M2 = *this;
	Matrix<T> Result(M2.numberOfRows, M2.numberOfColumns);
	T tmp;

	for(int row=0; row<M2.numberOfRows; row++)
		for(int col=0; col<M2.numberOfColumns; col++)
	{
		tmp = 0;
		vector<T> row_vector = M2.getRow(row);
		vector<T> col_vector = M.getColumn(col);
		for(unsigned int x=0; x<row_vector.size(); x++)
	  		tmp += (row_vector[x]*col_vector[x]);
		Result[row][col] = tmp;	
	}

	return Result;
}

template <class T>
vector<T> Matrix<T>::operator* (const vector<T> &v)
{
	vector<T> Result(numberOfRows);
	if(numberOfColumns != (int)v.size()) return Result;

	for(int row=0; row<numberOfRows; row++)
	{
		Result[row] = 0;
		for(int j=0;j<numberOfColumns;j++) Result[row] =  Result[row] + elements[row][j]*v[j];
	}
	return Result;
}
template <class T>
Matrix<T> Matrix<T>::multDiag(const vector<T> &diag)
{
	if(diag.size() != (unsigned int)numberOfColumns) return  Matrix<T> (numberOfRows, numberOfColumns);
	Matrix<T> Result = *this;

	for(int row=0; row<numberOfRows; row++)
		for(int col=0; col<numberOfColumns; col++)
	{
		Result[row][col] =  Result[row][col]*diag[col];	
	}

	return Result;
}

template <class T>
T Matrix<T>::cofactorSign(int row, int col)
{
	return (T) (pow(-1.0,double(row)) * pow(-1.0,double(col)));
}

template <class T>
Matrix<T> Matrix<T>::minminor(int row, int col)
{
	Matrix<T> result((getNumberOfRows()-1), (getNumberOfColumns()-1));
	int xpt=0, ypt=0;

	for(int y=0; y<(getNumberOfRows()); y++)
	{
		if(y != row)
		{
			for(int x=0; x<(getNumberOfColumns()); x++)
			{
	      			if(x != col)
				{
		  			result[ypt][xpt] = elements[y][x];
		  			xpt++;
				}
	    		}
	  		ypt++;
		}
		xpt=0;
	}
  
	return result;
}
	    
template <class T>
T Matrix<T>::determinant()
{
	T result = 0; 
	Matrix<T> tmp;

	if(getNumberOfRows() == 2)
	{
      		result = (elements[0][0]*elements[1][1] - elements[0][1]*elements[1][0]);
	}
	else
	{
      		for(int col=0; col < getNumberOfColumns(); col++)
		{
	  		tmp = minminor(0, col);
	  		result += (tmp.determinant() * elements[0][col] * cofactorSign(0, col));
		}
	}
	return result;  
}
template<class T>  
T Matrix<T>::operator()(const int x, const int y) const
{
	if(checkBounds(x, y))
		return elements[x][y];
}


template<class T>  
T& Matrix<T>::operator()(const int x, const int y)
{
	if(checkBounds(x, y))
		return elements[x][y];
	else
	{
		add(x+1, y+1);
		return elements[x][y];
	}
}
template<class T>  
inline bool Matrix<T>::checkBounds(int a, int b) const
{
	if(a < numberOfRows || b < numberOfColumns)
		return 1;
	else
		return 0;
}
template<class T>  
void Matrix<T>::add(int a, int b)
{
	if(a > numberOfRows)
	{
		numberOfRows = a;
		elements.resize(numberOfRows);
	}
	if(b > numberOfColumns)
	{
		numberOfColumns = b;
		for(int i = 0; i < numberOfRows; i++)
			elements[i].resize(numberOfColumns);
	}
}
template<class T>  
std::vector<T> Matrix<T>::slice(int row)
{
	std::vector<T> vec;
	if(row < numberOfRows)
		vec = elements[row];
	//cout<<"Number of rows : "<<row <<endl;
	//cout<<"Number of numberofrows : "<<numberOfRows <<endl;
	//cout<<"sizevect : "<<vec.size()<<endl;
	return vec;		
}

template<class T>
void Matrix<T>::clear()
{
	Matrix<T> M(0, 0);
    *this = M;
}
template<class T>  
void Matrix<T>::resize(int n)
{
	elements.resize(n);
	//cout <<"Rows : "<<getNumberOfRows()<<endl;
	//cout <<"Columns : "<<getNumberOfColumns()<<endl;
}
template class Matrix<int>;
template class Matrix<double>;
//template class Matrix<bool>;
template class Matrix<cdouble>;

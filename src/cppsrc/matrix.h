#include <math.h>
#include <iostream>
using namespace std;
class Matrix
{
private:
	int row,col;
	double **p;
public:

	Matrix(int a,int b)
	{ 
		this->row=a;
		this->col=b;
		this->p=new double*[row];
		for(int i=0;i<a;i++)
			this->p[i]=new double[col];
		for(int i=0;i<row;i++)
			for(int j=0;j<col;j++)
			{
				this->p[i][j]=0;
			}
	}
	Matrix(const Matrix &src )
	{
		this->row=src.row;
		this->col=src.col;
		this->p=new double*[row];
		for(int i=0;i<this->row;i++)
			this->p[i]=new double[col];
		for(int i=0;i<row;i++)
			for(int j=0;j<col;j++)
			{
				this->p[i][j]=src.p[i][j];
			}
	}

	    ~Matrix()//析构函数  
    {  
        for(int i=0;i<row;i++)  
            delete []p[i];  
        delete []p; 
		p=0;
    }  
		Matrix operator=(const Matrix& m1);
	    Matrix operator+(const Matrix& m1);
		Matrix operator-(const Matrix& m1);
		Matrix operator*(const Matrix& m1);
		double* operator[](int row);
		Matrix trans();
		Matrix inverse();
		double dist();
		void Input();
		void Display();
		void change(int row,int col,double a);
		void swap(bool a,int b1,int b2);
		double get(int row, int col);
};
double Matrix::dist()
{
	double n=0;
	for(int i=0;i<row;i++)
		for(int j=0;j<col;j++)
		{
			n=p[i][j]*(p[i][j])+n;
		}
		n=sqrt(n)/(row*col);
		return n;
}
double* Matrix::operator[](int nrow)
{
	if(nrow<this->row)
		return (double*)this->p[nrow];
}
void Matrix::Input()
{
	for(int i=0;i<row;i++)
		for(int j=0;j<col;j++)
		{
			cin>>p[i][j];
		}
}
void Matrix::Display()
{
	for(int i=0;i<row;i++)
	{	
		for(int j=0;j<col;j++)
		{
			cout<<p[i][j]<<" ";
		}
		cout<<endl;
	}
	return ;
}
void Matrix::change(int row,int col,double a)
{
	this->p[row][col]=a;
}
double Matrix::get(int row, int col)
{
	return this->p[row][col];

}
Matrix Matrix::operator=(const Matrix& m1)
{
	int i,j;
	if(p!=0)
	{
		for(i=0;i<row;i++)  
			delete []p[i];  
		delete []p;
	}
    row=m1.row;
	col=m1.col;  
	this->p=new double*[row];
	for(int i=0;i<this->row;i++)
		this->p[i]=new double[col];
	for(i=0;i<row;i++)
		for(j=0;j<col;j++)
		{
			p[i][j]=m1.p[i][j];
		}
		return *this;
}
Matrix Matrix::operator+(const Matrix& m1)
{
	if(this->row!=m1.row||this->col!=m1.col)
	{
		cout<<"行列不同的矩阵无法进行相加运算\n"<<endl;
	}
	Matrix m2(this->row,this->col);
	for(int i=0;i<this->row;i++)
		for(int j=0;j<this->col;j++)
			m2.p[i][j]=this->p[i][j]+m1.p[i][j];
	return m2;
}
Matrix Matrix::operator-(const Matrix& m1)
{
	if(this->row!=m1.row||this->col!=m1.col)
	{
		cout<<"行列不同的矩阵无法进行相减运算\n"<<endl;
	}
	Matrix m2(this->row,this->col);
	for(int i=0;i<this->row;i++)
		for(int j=0;j<this->col;j++)
			m2.p[i][j]=this->p[i][j]-m1.p[i][j];
	return m2;
}
Matrix Matrix::operator*(const Matrix& m1)
{
	if(this->col!=m1.row)
	{
		cout<<"不匹配的矩阵无法进行相乘运算\n"<<endl;
	}
	Matrix m2(this->row,m1.col);
	for(int i=0;i<m2.row;i++)
		for(int j=0;j<m2.col;j++)
		{
			for(int k=0;k<this->col;k++)
				m2.p[i][j]+=this->p[i][k]*m1.p[k][j];
		}
		return m2;
}
Matrix Matrix::trans()
{
	Matrix t(this->col,this->row);
	for(int i=0;i<this->row;i++)
		for(int j=0;j<this->col;j++)
		{
			t.p[j][i]=this->p[i][j];
		}
		return t;
}
void Matrix::swap(bool a,int b1,int b2)
{
	double p;
	if(a)
	{
		for(int i=0;i<this->col;i++)
		{
			p=this->p[b1][i];
			this->p[b1][i]=this->p[b2][i];
			this->p[b2][i]=p;
		}
	}
	else
	{
		for(int i=0;i<this->row;i++)
		{
			p=this->p[i][b1];
			this->p[i][b1]=this->p[i][b2];
			this->p[i][b2]=p;
		}
	}
}
Matrix Matrix::inverse()
{
	if(this->row!=this->col)
	{
		cout<<"奇异矩阵不可求逆"<<endl;
	}
	double d; 
	int n=this->row;
	int *row1 = new int[n];
	int *col1 = new int[n];
	Matrix A(*this);
	 for (int k=0; k<n;k++){ row1[k]=k; col1[k]=k; }
	    for (int k = 0; k <n; k++)  
    {  
		double dmax=0;
		//开始全选主元
		for(int r=k;r<n;r++)
		{
			for(int c=k;c<n;c++)
			{
				if(fabs(A.p[r][c])>dmax)
				{
					     dmax = fabs(A.p[r][c]);
						 row1[k]=r;
						 col1[k]=c;
				}
			}
		}
			  if (row1[k] !=k)    //交换行
				  A.swap(true,k, row1[k]);
			if (col1[k] !=k)    //交换列
				A.swap(false,k,col1[k]);
		d = 1.0/A.p[k][k];  
        A.p[k][k] = d;  
        for (int i = 0; i <n; i++)  
        {  
            if ( i != k )  
                A.p[k][i] *=-d;   
        }  
        for (int i = 0; i < n; i++)  
        {  
            if ( i != k)  
                A.p[i][k]*=d;  
        }  
        for (int i = 0; i < n; i++)  
        {  
            if ( i != k )  
            {  
                for (int j = 0; j < n; j++)  
                {  
                    if ( j != k )  
                        A.p[i][j] += A.p[i][k]*A.p[k][j]/d;   
                }  
            }  
        }  
    }  
 for (int k=n-1; k>=0; k--)  //恢复原来行列顺序
 {
  if (row1[k] != k )  //原来的行交换就列交换来恢复
   A.swap(false,k,row1[k]);
  if (col1[k] != k)  //原来的列交换就行交换来恢复
   A.swap(true,k, col1[k]);
 }
 free(row1); free(col1);
    return A;

}

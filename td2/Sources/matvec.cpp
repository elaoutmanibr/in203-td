// Produit matrice-vecteur
# include <cassert>
# include <mpi.h>
# include <vector>
# include <iostream>

// ---------------------------------------------------------------------
class Matrix : public std::vector<double>
{
public:
    Matrix (int dim);
    Matrix( int nrows, int ncols );
    Matrix( const Matrix& A ) = delete;
    Matrix( Matrix&& A ) = default;
    ~Matrix() = default;

    Matrix& operator = ( const Matrix& A ) = delete;
    Matrix& operator = ( Matrix&& A ) = default;
    
    double& operator () ( int i, int j ) {
        return m_arr_coefs[i + j*m_nrows];
    }
    double  operator () ( int i, int j ) const {
        return m_arr_coefs[i + j*m_nrows];
    }
    
    std::vector<double> operator * ( const std::vector<double>& u ) const;
    
    std::ostream& print( std::ostream& out ) const
    {
        const Matrix& A = *this;
        out << "[\n";
        for ( int i = 0; i < m_nrows; ++i ) {
            out << " [ ";
            for ( int j = 0; j < m_ncols; ++j ) {
                out << A(i,j) << " ";
            }
            out << " ]\n";
        }
        out << "]";
        return out;
    }
private:
    int m_nrows, m_ncols;
    std::vector<double> m_arr_coefs;
};
// ---------------------------------------------------------------------
inline std::ostream& 
operator << ( std::ostream& out, const Matrix& A )
{
    return A.print(out);
}
// ---------------------------------------------------------------------
inline std::ostream&
operator << ( std::ostream& out, const std::vector<double>& u )
{
    out << "[ ";
    for ( const auto& x : u )
        out << x << " ";
    out << " ]";
    return out;
}
// ---------------------------------------------------------------------
std::vector<double> 
Matrix::operator * ( const std::vector<double>& u ) const
{
	
    const Matrix& A = *this;
    assert( u.size() == unsigned(m_ncols) );
    std::vector<double> v(m_nrows, 0.);
    for ( int i = 0; i < m_nrows; ++i ) {
        for ( int j = 0; j < m_ncols; ++j ) {
            v[i] += A(i,j)*u[j];
        }            
    }
    return v;
}

// =====================================================================
Matrix::Matrix (int dim) : m_nrows(dim), m_ncols(dim),
                           m_arr_coefs(dim*dim)
{
    for ( int i = 0; i < dim; ++ i ) {
        for ( int j = 0; j < dim; ++j ) {
            (*this)(i,j) = (i+j)%dim;
        }
    }
}
// ---------------------------------------------------------------------
Matrix::Matrix( int nrows, int ncols ) : m_nrows(nrows), m_ncols(ncols),
                                         m_arr_coefs(nrows*ncols)
{
    int dim = (nrows > ncols ? nrows : ncols );
    for ( int i = 0; i < nrows; ++ i ) {
        for ( int j = 0; j < ncols; ++j ) {
            (*this)(i,j) = (i+j)%dim;
        }
    }    
}
// =====================================================================
std::vector<double> multip_vect_col(Matrix& A, std::vector<double>& u, int rank, int Nloc){
	
    std::vector<double> v(u.size(), 0.);
    for ( unsigned int i = 0; i < u.size(); ++i ) {
        for ( int j = rank*Nloc; j < (rank+1)*Nloc; ++j ) {
            v[i] += A(i,j)*u[j];
        }            
    }
    return v;
    
}
// =====================================================================

// =====================================================================
std::vector<double> multip_vect_lig(Matrix& A, std::vector<double>& u, int rank, int Nloc){
	
    std::vector<double> v(Nloc, 0.);
    for ( int i = rank*Nloc; i < (rank+1)*Nloc; ++i ) {
        for (unsigned int j = 0; j < u.size(); ++j ) {
            v[i] += A(i,j)*u[j];
        }            
    }
    return v;
    
}
// =====================================================================

int main( int nargs, char* argv[] )
{
	MPI_Init( &nargs, &argv );

	MPI_Comm globComm;
	MPI_Comm_dup(MPI_COMM_WORLD, &globComm);

	int nbp;
	MPI_Comm_size(globComm, &nbp);
	
	int rank;
	MPI_Comm_rank(globComm, &rank);
	
    const int N = 120;
    const int Nloc = N/nbp; // nbp divise N par hypothèse
    Matrix A(N);
    std::vector<double> u(N);
    std::vector<double> r(N);
    for ( int i = 0; i < N; ++i ) u[i] = i+1;
    
    //std::cout  << "A : " << A << std::endl;

    
    //std::cout << " u : " << u << std::endl;
    //std::vector<double> v1 = A*u;
    std::vector<double> v = multip_vect_col(A,u,rank,Nloc);
    MPI_Allreduce(&v, &r, N, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    
    //std::cout << "A.u = " << v1 << std::endl;
    //std::cout << "A.u = " << r << std::endl;
    MPI_Finalize();
    return EXIT_SUCCESS;
}

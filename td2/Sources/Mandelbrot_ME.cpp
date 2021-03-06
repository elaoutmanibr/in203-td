# include <mpi.h>
# include <iostream>
# include <cstdlib>
# include <string>
# include <chrono>
# include <cmath>
# include <vector>
# include <fstream>


struct Complex
{
    Complex() : real(0.), imag(0.)
    {}
    Complex(double r, double i) : real(r), imag(i)
    {}
    Complex operator + ( const Complex& z )
    {
        return Complex(real + z.real, imag + z.imag );
    }
    Complex operator * ( const Complex& z )
    {
        return Complex(real*z.real-imag*z.imag, real*z.imag+imag*z.real);
    }
    double sqNorm() { return real*real + imag*imag; }
    double real,imag;
};

std::ostream& operator << ( std::ostream& out, const Complex& c )
{
  out << "(" << c.real << "," << c.imag << ")" << std::endl;
  return out;
}

/** Pour un c complexe donné, calcul le nombre d'itérations de mandelbrot
 * nécessaires pour détecter une éventuelle divergence. Si la suite
 * converge, la fonction retourne la valeur maxIter
 **/
int iterMandelbrot( int maxIter, const Complex& c)
{
    Complex z{0.,0.};
    // On vérifie dans un premier temps si le complexe
    // n'appartient pas à une zone de convergence connue :
    // Appartenance aux disques  C0{(0,0),1/4} et C1{(-1,0),1/4}
    if ( c.real*c.real+c.imag*c.imag < 0.0625 )
        return maxIter;
    if ( (c.real+1)*(c.real+1)+c.imag*c.imag < 0.0625 )
        return maxIter;
    // Appartenance à la cardioïde {(1/4,0),1/2(1-cos(theta))}    
    if ((c.real > -0.75) && (c.real < 0.5) ) {
        Complex ct{c.real-0.25,c.imag};
        double ctnrm2 = sqrt(ct.sqNorm());
        if (ctnrm2 < 0.5*(1-ct.real/ctnrm2)) return maxIter;
    }
    int niter = 0;
    while ((z.sqNorm() < 4.) && (niter < maxIter))
    {
        z = z*z + c;
        ++niter;
    }
    return niter;
}

/**
 * On parcourt chaque pixel de l'espace image et on fait correspondre par
 * translation et homothétie une valeur complexe c qui servira pour
 * itérer sur la suite de Mandelbrot. Le nombre d'itérations renvoyé
 * servira pour construire l'image finale.
 
 Sortie : un vecteur de taille W*H avec pour chaque case un nombre d'étape de convergence de 0 à maxIter
 MODIFICATION DE LA FONCTION :
 j'ai supprimé le paramètre W étant donné que maintenant, cette fonction ne prendra plus que des lignes de taille W en argument.
 **/
void 
computeMandelbrotSetRow( int W, int H, int maxIter, int num_ligne, int* pixels)
{
    // Calcul le facteur d'échelle pour rester dans le disque de rayon 2
    // centré en (0,0)
    double scaleX = 3./(W-1);
    double scaleY = 2.25/(H-1.);
    //
    // On parcourt les pixels de l'espace image :
    for ( int j = 0; j < W; ++j ) {
       Complex c{-2.+j*scaleX,-1.125+ num_ligne*scaleY};
       pixels[j] = iterMandelbrot( maxIter, c );
    }
}

std::vector<int>
computeMandelbrotSet(int W, int H, int maxIter, int b, int e, int Hloc)
{
    std::chrono::time_point<std::chrono::system_clock> start, end;
    std::vector<int> pixels(W*Hloc);
    start = std::chrono::system_clock::now();
    // On parcourt les pixels de l'espace image :
    for ( int i = b; i < e; ++i ) {
      computeMandelbrotSetRow(W, H, maxIter, i, pixels.data() + W*(Hloc-i+b-1) );
    }
    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;
    std::cout << "Temps calcul ensemble mandelbrot : " << elapsed_seconds.count() 
              << std::endl;
    return pixels;
}

/** Construit et sauvegarde l'image finale **/
void savePicture( const std::string& filename, int W, int H, const std::vector<int>& nbIters, int maxIter )
{
    double scaleCol = 1./maxIter;//16777216
    std::ofstream ofs( filename.c_str(), std::ios::out | std::ios::binary );
    ofs << "P6\n"
        << W << " " << H << "\n255\n";
    for ( int i = 0; i < W * H; ++i ) {
        double iter = scaleCol*nbIters[i];
        unsigned char r = (unsigned char)(256 - (unsigned (iter*256.) & 0xFF));
        unsigned char b = (unsigned char)(256 - (unsigned (iter*65536) & 0xFF));
        unsigned char g = (unsigned char)(256 - (unsigned( iter*16777216) & 0xFF));
        ofs << r << g << b;
    }
    ofs.close();
}

const int W = 800;
const int H = 600;




int main( int nargs, char* argv[] )
{
	MPI_Init( &nargs, &argv );

	MPI_Comm globComm;
	MPI_Comm_dup(MPI_COMM_WORLD, &globComm);

	int nbp;
	MPI_Comm_size(globComm, &nbp);
	
	int rank;
	MPI_Comm_rank(globComm, &rank);
	MPI_Status status;
	int tag = 500;
	//MPI_Status status;
	
	
    const int maxIter = 8*65536;
    // std::vector<int> iters(W*Hloc); 
    
    std::vector<int> iters(W);
	
	
	//~ iters = computeMandelbrotSet( W,H, maxIter, (nbp-rank-1)*Hloc,(nbp-rank)*Hloc,Hloc  );
	//~ MPI_Gather(iters.data(), W * Hloc, MPI_INT, res.data(), W * Hloc, MPI_INT, 0, MPI_COMM_WORLD);

	//~ if (rank==0){
		//~ savePicture("mandelbrot.tga", W, H, res, maxIter);
	//~ }
	
	
	if (rank==0){
		std::vector<int> res(W*H);
		int row =0;
		
		for (int dest =1 ;dest<nbp; dest++){
			MPI_Send(&row,1,MPI_INT,dest,tag,MPI_COMM_WORLD);
			row++;
		}
		
		while (row<H){
			int ligne;
			MPI_Recv(iters.data(),W,MPI_INT,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
			
			if (status.MPI_SOURCE != 0){
				
				ligne = status.MPI_TAG;
				for (int k =0; k<W; k++){
					res[ligne*W+k]= iters[k];
				}
				MPI_Send(&row,1,MPI_INT,status.MPI_SOURCE,tag,MPI_COMM_WORLD);
				row++;
			}
		}
		

		row =-1; // termination signal
		for (int i=0;i<nbp;i++){
			MPI_Send(&row,1,MPI_INT,i,tag,MPI_COMM_WORLD);
		}
		
		savePicture("mandelbrot.tga", W, H, res, maxIter);
		
	}else{
		int command = 0;
		MPI_Recv(&command,1,MPI_INT,0,tag,MPI_COMM_WORLD,&status);
		while (command != -1){
			
			computeMandelbrotSetRow(W, H, maxIter, command, iters.data());
			std::cout << rank << "\t" << command << std::endl;
			MPI_Send(iters.data(),W,MPI_INT,0,command,MPI_COMM_WORLD);
			MPI_Recv(&command,1,MPI_INT,0,tag,MPI_COMM_WORLD,&status);
			std::cout << rank << "\t" << command << std::endl;
		}
	}
	

	// A la fin du programme, on doit synchroniser une dernière fois tous les processus
	// afin qu'aucun processus ne se termine pendant que d'autres processus continue à
	// tourner. Si on oublie cet instruction, on aura une plantage assuré des processus
	// qui ne seront pas encore terminés.
	MPI_Finalize();
	return EXIT_SUCCESS;
}


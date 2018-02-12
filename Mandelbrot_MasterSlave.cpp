# include <iostream>
# include <cstdlib>
# include <string>
# include <chrono>
# include <cmath>
# include "lodepng/lodepng.h"
# include <mpi.h>
# include <fstream>
# include <sstream>
# include <iomanip>

/** Une structure complexe est définie pour la bonne raison que la classe
 * complex proposée par g++ est très lente ! Le calcul est bien plus rapide
 * avec la petite structure donnée ci--dessous
 **/
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
 **/
std::vector<int> computeMandelbrotSet( int W, int H, int maxIter, int rank, int H_loc, std::ofstream& output)
{
    std::chrono::time_point<std::chrono::system_clock> start, end;
    // Calcul le facteur d'échelle pour rester dans le disque de rayon 2
    // centré en (0,0)
    double scaleX = 3./(W-1);
    double scaleY = 2.25/(H-1);
    //
    std::vector<int> pixels(W*H_loc);
    start = std::chrono::system_clock::now();
    // On parcourt les pixels de l'espace image :
    for ( int i_loc = 0; i_loc < H_loc; ++i_loc )
    {
        int i_glob = i_loc + rank * H_loc;
        for ( int j = 0; j < W; ++j ) {
            Complex c{-2.+j*scaleX,-1.125+i_glob * scaleY};
            pixels[i_loc*W+j] = iterMandelbrot( maxIter, c );
        }
    }
    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;
    output << "Temps calcul ensemble mandelbrot row " 
    << rank * H_loc << " to row " << (rank +1)*H_loc-1 << ": " << elapsed_seconds.count() << std::endl;
    return pixels;
}

/** Construit et sauvegarde l'image finale **/
void savePicture( const std::string& fileName, int W, int H, const std::vector<int>& nbIters, int maxIter )
{
    std::vector<unsigned char> image(4*W*H);
    double scaleCol = 1./maxIter;//16777216
    for ( int i = 0; i < H; ++i ) {
        for ( int j = 0; j < W; ++j ) {
            double iter = scaleCol*nbIters[i*W+j];
            unsigned r = unsigned (iter*256.) & 0xFF;
            unsigned b = (unsigned (iter*65536) & 0xFF);
            unsigned g = (unsigned( iter*16777216) & 0xFF);
            image[4*(i*W+j)+0] = (unsigned char)(256-r);
            image[4*(i*W+j)+1] = (unsigned char)(256-g);
            image[4*(i*W+j)+2] = (unsigned char)(256-b);
            image[4*(i*W+j)+3] = 255;
        }
    }
    unsigned error = lodepng::encode(fileName.c_str(), image, W, H);

    //if there's an error, display it
    if(error) std::cout << "encoder error " << error << ": "<< lodepng_error_text(error) << std::endl;
}

int main(int nargs, char** argv)
{
    const int W = 800;
    const int H = 600;
    // Normalement, pour un bon rendu, il faudrait le nombre d'itérations
    // ci--dessous :
    //const int maxIter = 16777216;
    const int maxIter = 8*65536;
    MPI_Init( &nargs, &argv );
    MPI_Comm globComm;
    MPI_Comm_dup(MPI_COMM_WORLD, &globComm);
    int nbp;
    MPI_Comm_size(globComm, &nbp);
    int rank;
    MPI_Comm_rank(globComm, &rank);
    std::stringstream fileName;
    fileName << "Output" << std::setfill('0') << std::setw(5) << rank << ".txt";
    std::ofstream output( fileName.str().c_str() );

    output << "I'm the processus " << rank << " on " << nbp << " processes." << std::endl;

    std::vector<int> pixels(W*H);

    if(rank == 0) // Master: patch line task to slaves
    {
        int currentRow[W]; // data
        int nbRowsSent = 0; // number of lines already sent
        int nbRowsRecv = 0; // number of rows already received
        MPI_Status currentStatus;
        for(int rk = 0; rk < nbp-1; rk++)
        {
            MPI_Send(&rk, 1, MPI_INT, rk+1, 0, globComm); // send first tasks
            nbRowsSent++;
        }
        
        while(nbRowsRecv < H)
        {
            MPI_Recv(&currentRow, W, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, globComm, &currentStatus);
            nbRowsRecv++;
            int currentRowNum = currentStatus.MPI_TAG;
            int slave_rk = currentStatus.MPI_SOURCE;
            pixels.insert(pixels.begin()+W*(currentRowNum), currentRow, currentRow+W);
            if(nbRowsSent < H)
            {
                MPI_Send(&nbRowsSent, 1, MPI_INT, slave_rk, 0, globComm); // send next line
                nbRowsSent++;
            }
            else
            {
                int finishSignal = -1;
                MPI_Send(&finishSignal, 1, MPI_INT, slave_rk, 0, globComm); // send next line
            }
        }
    }
    else // Slave: receive instruction and execute
    {
        int row_recv = 0;
        while(row_recv != -1)
        {
            MPI_Recv(&row_recv, 1, MPI_INT, 0, 0, globComm, NULL);
            if(row_recv != -1)
            {
                auto iters = computeMandelbrotSet( W, H, maxIter, row_recv, 1, output); // compute only one line
                MPI_Send(iters.data(), W, MPI_INT, 0, row_recv, globComm);
            }  
        }
    }
    if ( rank == 0 )
    {
        output << "Master finished, saving image ..." << std::endl;
        savePicture("mandelbrot_MasterSlave.png", W, H, pixels, maxIter);
    }

    output.close();

    MPI_Finalize();

    return EXIT_SUCCESS;
}
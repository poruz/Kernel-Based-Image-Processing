#include <iostream>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <cstdlib>
#include "bmplib.h"

using namespace std;

//============================Add function prototypes here======================
void dummy(unsigned char[][SIZE][RGB], unsigned char[][SIZE][RGB]);
void convolve(unsigned char[][SIZE][RGB], unsigned char[][SIZE][RGB], 
	      int, double[][11]);
void sobel(unsigned char[][SIZE][RGB], unsigned char[][SIZE][RGB]);
void gaussian(double [][11], int, double);
void gaussian_filter(unsigned char[][SIZE][RGB], unsigned char[][SIZE][RGB], int,
         double);
void unsharp(unsigned char[][SIZE][RGB], unsigned char[][SIZE][RGB], int,
         double, double);
//============================Do not change code in main()======================

#ifndef AUTOTEST

int main(int argc, char* argv[])
{
   //First check argc
  if(argc < 3)
    {
      //we need at least ./filter <input file> <filter name> to continue
      cout << "usage: ./filter <input file> <filter name> <filter parameters>";
      cout << " <output file name>" << endl;
      return -1;
    }
   //then check to see if we can open the input file
   unsigned char input[SIZE][SIZE][RGB];
   unsigned char output[SIZE][SIZE][RGB];
   char* outfile;
   int N;
   double sigma, alpha;

   // read file contents into input array
   int status = readRGBBMP(argv[1], input); 
   if(status != 0)
   {
      cout << "unable to open " << argv[1] << " for input." << endl;
      return -1;
   }
   //Input file is good, now look at next argument
   if( strcmp("sobel", argv[2]) == 0)
   {
     sobel(output, input);
     outfile = argv[3];
   }
   else if( strcmp("blur", argv[2]) == 0)
   {
     if(argc < 6)
       {
	 cout << "not enough arguments for blur." << endl;
	 return -1;
       }
     N = atoi(argv[3]);
     sigma = atof(argv[4]);
     outfile = argv[5];
     gaussian_filter(output, input, N, sigma);
   }
   else if( strcmp("unsharp", argv[2]) == 0)
   {
     if(argc < 7)
       {
	 cout << "not enough arguments for unsharp." << endl;
	 return -1;
       }
     N = atoi(argv[3]);
     sigma = atof(argv[4]);
     alpha = atof(argv[5]);
     outfile = argv[6];
     unsharp(output, input, N, sigma, alpha);

   }
   else if( strcmp("dummy", argv[2]) == 0)
   {
     //do dummy
     dummy(output, input);
     outfile = argv[3];
   }
   else
   {
      cout << "unknown filter type." << endl;
      return -1;
   }

   if(writeRGBBMP(outfile, output) != 0)
   {
      cout << "error writing file " << argv[3] << endl;
   }

}   

#endif 

//=========================End Do not change code in main()=====================


// Creates an identity kernel (dummy kernel) that will simply
// copy input to output image and applies it via convolve()
//
// ** This function is complete and need not be changed.
// Use this as an example of how to create a kernel array, fill it in
// appropriately and then use it in a call to convolve.
void dummy(unsigned char out[][SIZE][RGB], unsigned char in[][SIZE][RGB])
{
   double k[11][11];
   for (int i = 0; i < 3; i++)
   {
      for(int j = 0; j < 3; j++)
      {
         k[i][j] = 0;
      }
   }
   k[1][1] = 1;
   convolve(out, in, 3, k);
}


// Convolves an input image with an NxN kernel to produce the output kernel
// You will need to complete this function by following the 
//  instructions in the comments
void convolve(unsigned char out[][SIZE][RGB], unsigned char in[][SIZE][RGB], 
	      int N, double kernel[][11])
{
 
   int padded[SIZE+11][SIZE+11][RGB];  // Use for input image with appropriate 
                                       // padding
   int temp[SIZE][SIZE][RGB];          // Use for the unclamped output pixel 
                                       // values then copy from temp to out, 
                                       // applying clamping 

   //first set all of padded to 0 (black)
   for(int i = 0; i < SIZE+11; i++)
      for(int j = 0; j < SIZE+11; j++)
         for(int k = 0; k < RGB; k++)
          {
            padded[i][j][k] = 0;
          }


   //now copy input into padding to appropriate locations
   for(int i = 0; i < SIZE; i++)
      for(int j = 0; j < SIZE; j++)
         for(int k = 0; k < RGB; k++)
         {
            padded[i+N/2][j+N/2][k] = in[i][j][k];
         }


   //initialize temp pixels to 0 (black)
   for(int i = 0; i < SIZE; i++)
      for(int j = 0; j < SIZE; j++)
         for(int k = 0; k < RGB; k++)
         {
            temp[i][j][k] = 0;
         }


   //now perform convolve (using convolution equation on each pixel of the 
   // actual image) placing the results in temp (i.e. unclamped result)
   for(int x = 0; x < SIZE; x++)
      for(int y = 0; y < SIZE; y++)
         for(int z = 0; z < RGB; z++)
            for(int i = -N/2; i <= N/2; i++)
               for(int j = -N/2; j <= N/2; j++)
               {
                  temp[y][x][z] += padded[y+i+N/2][x+j+N/2][z]
                                       * kernel[N/2 + i][N/2 + j];
               }

   //now clamp and copy to output
   // You may need to cast to avoid warnings from the compiler:
   // (i.e. out[i][j][k] = (unsigned char) temp[i][j][k];)
   for(int i = 0; i < SIZE; i++)
      for(int j = 0; j < SIZE; j++)
         for(int k = 0; k < RGB; k++)
         {
            if(temp [i][j][k] < 0)
               temp[i][j][k] = 0;
            else if(temp[i][j][k] > 255)
               temp[i][j][k] = 255;
            out[i][j][k] = (unsigned char) temp[i][j][k];      
         }



}

// You will need to complete this function by following the 
//  instructions in the comments
void sobel(unsigned char out[][SIZE][RGB], unsigned char in[][SIZE][RGB])
{
   double k[11][11];
   double s_h1[3][3] = { {-1, 0, 1}, 
                         {-2, 0, 2}, 
                         {-1, 0, 1} };
   double s_h2[3][3] = { {1, 0, -1}, 
                         {2, 0, -2}, 
                         {1, 0, -1} };
   
   unsigned char h1_soble[SIZE][SIZE][RGB]; //hold intemediate images
   unsigned char h2_soble[SIZE][SIZE][RGB]; 

   for (int i = 0; i < 11; i++)
   {
      for(int j=0; j < 11; j++)
      {
         k[i][j] = 0;
      }
   }


   // Copy in 1st 3x3 horizontal sobel kernel (i.e. copy s_h1 into k)
   for (int i = 0; i < 3; i++)
   {
      for(int j=0; j < 3; j++)
      {
         k[i][j] = s_h1[i][j];
      }
   }


   // Call convolve to apply horizontal sobel kernel with result in h1_soble
   convolve(h1_soble, in, 3, k);


   // Copy in 2nd 3x3 horizontal sobel kernel (i.e. copy s_h2 into k)
   for (int i = 0; i < 3; i++)
   {
      for(int j=0; j < 3; j++)
      {
         k[i][j] = s_h2[i][j];
      }
   }


   // Call convolve to apply horizontal sobel kernel with result in h2_soble
   convolve(h2_soble, in, 3, k);


   // Add the two results (applying clamping) to produce the final output 
   for(int r = 0; r < SIZE; r++)
      for(int c = 0; c < SIZE; c++)
         for(int h = 0; h < RGB; h++)
         {
            int temp = (int)(h1_soble[r][c][h] + h2_soble[r][c][h]);
            if(temp < 0)
               temp = 0;
            else if(temp > 255)
               temp = 255;
            out[r][c][h] = (unsigned char)temp;
         }
}


// Add the rest of your functions here (i.e. gaussian, gaussian_filter, unsharp)

void gaussian(double kernel[][11], int N, double sigma)
{
   //Formula for raw,
   //g(x,y) = A*e^(-((x-x0)^2/(2*sigma^2) + (y-y0)^2/(2*sigma^2)))
   //A = 1, x0 = y0 = 0
   // => g(x,y) = e^(-(x^2/(2*sigma^2) + y^2/(2*sigma^2)))


   //Raw gaussian kernel
   double sum = 0;
   for(int i = -N/2; i <= N/2; i++)
      for(int j = -N/2; j <= N/2; j++)
         kernel[N/2 + i][N/2 + j] = exp(-1*(i*i/(2*sigma*sigma) 
                                             + j*j/(2*sigma*sigma)));
   
   //Sum of elements  
   for(int i = 0; i < N; i++)
      for(int j = 0; j < N; j++)
         sum += kernel[i][j];   
      
      
   //Normalized kernel and printing it in a nice format using setw
   for(int i = 0; i < N; i++)
   {
      for(int j = 0; j < N; j++)
      {
         kernel[i][j] = kernel[i][j]/sum;
         cout << setw(10) << kernel[i][j];
      }
      cout << endl;
   }  
}

void gaussian_filter(unsigned char output[][SIZE][RGB], 
         unsigned char input[][SIZE][RGB], int N, double sigma)
{
   double kernel[11][11];
   
   //Call gaussian to get the kernel
   gaussian(kernel, N, sigma);
   
   //Now we have the kernel, need to convolve now
   //Call convolve
   convolve(output, input, N, kernel);

}

void unsharp(unsigned char output[][SIZE][RGB], unsigned char input[][SIZE][RGB], 
         int N, double sigma, double alpha)
{
   unsigned char blur[SIZE][SIZE][RGB];
   int detail[SIZE][SIZE][RGB];
   double kernel[11][11];
   
   //Get the Blurred image
   gaussian(kernel,N,sigma);
   convolve(blur, input, N, kernel);
   
   //Get the details
   for(int i = 0; i < SIZE; i++)
      for(int j = 0; j < SIZE; j++)
         for(int k = 0; k < RGB; k++)
            detail[i][j][k] = (int)input[i][j][k] - (int)blur[i][j][k];

   //Get the sharpened version
   for(int i = 0; i < SIZE; i++)
      for(int j = 0; j < SIZE; j++)
         for(int k = 0; k < RGB; k++)
            {
               int temp = (int)input[i][j][k] + (int)(alpha*detail[i][j][k]);
               if (temp > 255)
                  temp = 255;
               else if(temp < 0)
                  temp = 0;
               output[i][j][k] = (unsigned char)temp; 
            }
}

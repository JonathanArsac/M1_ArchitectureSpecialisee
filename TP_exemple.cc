// Programme source exemple pour le TP d'Architectures specialis�es M1 informatique
// Universit� de Toulouse 3 - Paul Sabatier � 2017
//
// � compiler avec le Makefile fourni : make TP_exemple


#include <opencv/cv.h>
#include <opencv/highgui.h>
#include <iostream>
// #include <xmmintrin.h>	// SSE
// #include <emmintrin.h>	// SSE2
// #include <pmmintrin.h>	// SSE3
#include <tmmintrin.h>	// SSSE3

#include <iomanip>

using namespace std;

// affiche_xmm_i() -----------------------------------------------------//
// affiche le contenu du registre xmm entier sous diff�rents format :	//
//  8, 16 et 32 bits non sign�s.					//
//----------------------------------------------------------------------//

void affiche_xmm_i ( __m128i x )
{
  uint32_t a[4];
  _mm_store_si128((__m128i*)a,x);
  cout<<"Xmm = "<<a[0]<<" "<<a[1]<<" "<<a[2]<<" "<<a[3]<<" \n";

  uint16_t a2[8];
  _mm_store_si128((__m128i*)a2,x);
  cout<<"Xmm = "<<a2[0]<<" "<<a2[1]<<" "<<a2[2]<<" "<<a2[3]
  <<a2[4]<<" "<<a2[5]<<" "<<a2[6]<<" "<<a2[7]<<" \n";

  uint8_t a3[16];
  _mm_store_si128((__m128i*)a3,x);
  cout<<"Xmm = "<<a3[0]<<" "<<a3[1]<<" "<<a3[2]<<" "<<a3[3]
  <<a3[4]<<" "<<a3[5]<<" "<<a3[6]<<" "<<a3[7]
  <<a3[8]<<" "<<a3[9]<<" "<<a3[10]<<" "<<a3[11]
  <<a3[12]<<" "<<a3[13]<<" "<<a3[14]<<" "<<a3[15]<<" \n";

}


// Inverse() -----------------------------------------------------------//
// inverse les valeurs des pixels de l'image en niveaux de gris X	//
// de taille h*l et renvoie le r�sultat dans Y.				//
//----------------------------------------------------------------------//

void Inverse(unsigned char *X, long h, long l, unsigned char *Y)
{
    int nPixels = h*l;

   for (int i=0; i < nPixels; i++)
   {
	Y[i] = 255-X[i];
   }
}

void Inverse_sse(unsigned char *X, long h, long l, unsigned char *Y)
{
  int nbPixels = h*l;
  __m128i a,b,c;
  b=_mm_set1_epi8(255);
  for(int i=0;i<nbPixels;i+=16){
    a=_mm_load_si128((__m128i*)&X[i]);
    c=_mm_sub_epi8(b,a);
    _mm_store_si128((__m128i*)&Y[i],c);
  }

}


// Luminance_sse() -------------------------------------------------------------//
// r�duit ou augmente de la valeur du param�tre d_lum les valeurs des pixels de	//
// l'image en niveaux de gris X de taille h*l et renvoie le r�sultat dans Y.	//
//------------------------------------------------------------------------------//

void Luminance_sse(int d_lum, unsigned char *X, long h, long l, unsigned char *Y)
{
  int nbPixels=h*l;
  __m128i a,b,c;
  b=_mm_set1_epi8(d_lum);
  affiche_xmm_i(b);
  for(int i=0;i<nbPixels;i+=16){
    a=_mm_load_si128((__m128i*)&X[i]);
    c=_mm_adds_epu8(b,a);
    _mm_store_si128((__m128i*)&Y[i],c);
  }

}


// binarise() ----------------------------------------------------------//
// binarise l'image en niveaux de gris X selon le param�tre seuil	//
// et renvoie le r�sultat dans Y.					//
//----------------------------------------------------------------------//

void binarise (unsigned char seuil, unsigned char *X, long h, long l, unsigned char *Y)
{
    int nPixels = h*l;

    for (int i=0; i < nPixels; i++)
    {
	Y[i] = (X[i]>seuil)?255:0;
    }
}


void binarise_sse (unsigned char seuil, unsigned char *X, long h, long l, unsigned char *Y)
{
}



// min_et_max() ----------------------------------------------------------------//
// calcule les valeurs min et max des pixels de l'image en niveaux de gris X.	//
//------------------------------------------------------------------------------//

void min_et_max (unsigned char *X, long h, long l, unsigned char & min, unsigned char & max )
{
    int nPixels = h*l;
    min = 255;
    max = 0;

    for (int i=0; i < nPixels; i++)
    {
	if (X[i]<min)
	    min=X[i];
	if (X[i]>max)
	    max=X[i];
    }
}

void min_et_max_sse (unsigned char *X, long h, long l, unsigned char & min, unsigned char & max )
{
}



// recadrageDynamique() ------------------------------------------------//
// calcule le recadrage dynamique de l'image en niveaux de gris X	//
// et renvoie le r�sultat dans Y					//
//----------------------------------------------------------------------//

void recadrageDynamique (unsigned char *X, long h, long l, unsigned char min, unsigned char max, unsigned char *Y)
{
    long i;
    long size=h*l;

    //calcule des param�tres pour le recadrage dynamique
    int div=max-min;

    for(i=0;i<size;i++){
        Y[i]=((X[i]-min)*255)/div;
    }
}

void recadrageDynamique_sse(unsigned char *X, long h, long l, unsigned char min, unsigned char max, unsigned char *Y)
{
}



int main (int argc, char** argv)
{
  IplImage* img = NULL;

  if(argc!=2)
  {
    cerr << "Usage: " << argv[0] << " image" << endl;
    exit(1);
  }

  cout << "Ouverture de " << argv[1] << endl;
  img=cvLoadImage(argv[1], CV_LOAD_IMAGE_GRAYSCALE);	// chargement en m�moire de l'image depuis le fichier pass� en paramtre
  if(!img)
  {
    cerr << "Could not load image file: " << argv[1] << endl;
    exit(2);
  }
  cvShowImage( "Affiche_origine", img);
  cout << argv[1] << " : " << img->height << "x" << img->width << endl;

	// cr�ation d'une image pour stocker l'image invers�e
  IplImage* image_inverse = cvCreateImage (cvGetSize (img), IPL_DEPTH_8U, 1);
  Inverse((unsigned char*) img->imageData, img->height, img->width, (unsigned char*) image_inverse->imageData);
  cvShowImage( "Affiche_inverse", image_inverse);

	// cr�ation d'une image pour stocker l'image invers�e en sse
  IplImage* image_inverse_sse = cvCreateImage (cvGetSize (img), IPL_DEPTH_8U, 1);
  Inverse_sse((unsigned char*) img->imageData, img->height, img->width, (unsigned char*) image_inverse_sse->imageData);
  cvShowImage( "Affiche_inverse_sse", image_inverse_sse);

  // cr�ation d'une image pour stocker l'image illuminée en sse
  IplImage* image_luminance_sse = cvCreateImage (cvGetSize (img), IPL_DEPTH_8U, 1);
  Luminance_sse(20,(unsigned char*) img->imageData, img->height, img->width, (unsigned char*) image_inverse_sse->imageData);
  cvShowImage( "Affiche_Luminance_sse", image_luminance_sse);



	// cr�ation d'une image pour stocker l'image binaris�e
  IplImage* image_binaire = cvCreateImage (cvGetSize (img), IPL_DEPTH_8U, 1);
  binarise( 120, (unsigned char*) img->imageData, img->height, img->width, (unsigned char*) image_binaire->imageData);
  cvShowImage( "Affiche_binaire", image_binaire);




  unsigned char min, max;
  min_et_max( (unsigned char*) img->imageData, img->height, img->width, min, max);
  cout << "min : " << (int)min << "   max : " << (int)max << endl;

	// cr�ation d'une image pour stocker l'image de recadrage dynamique
  IplImage* image_ = cvCreateImage (cvGetSize (img), IPL_DEPTH_8U, 1);
  recadrageDynamique( (unsigned char*) img->imageData, img->height, img->width, min, max, (unsigned char*) image_->imageData);
  cvShowImage( "Affiche_", image_);



  cvWaitKey();		// attend une touche clavier
  exit (0);
}

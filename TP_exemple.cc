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
  cout<<"Xmm 32 = "<<a[0]<<" "<<a[1]<<" "<<a[2]<<" "<<a[3]<<" \n";

  uint16_t a2[8];
  _mm_store_si128((__m128i*)a2,x);
  cout<<"Xmm 16 = "<<a2[0]<<" "<<a2[1]<<" "<<a2[2]<<" "<<a2[3]<<" "
  <<a2[4]<<" "<<a2[5]<<" "<<a2[6]<<" "<<a2[7]<<" \n";

  uint8_t a3[16];
  _mm_store_si128((__m128i*)a3,x);
  cout<<"Xmm 8 = "<<(int)a3[0]<<" "<<(int)a3[1]<<" "<<(int)a3[2]<<" "<<(int)a3[3]<<" "
  <<(int)a3[4]<<" "<<(int)a3[5]<<" "<<(int)a3[6]<<" "<<(int)a3[7]<<" "
  <<(int)a3[8]<<" "<<(int)a3[9]<<" "<<(int)a3[10]<<" "<<(int)a3[11]<<" "
  <<(int)a3[12]<<" "<<(int)a3[13]<<" "<<(int)a3[14]<<" "<<(int)a3[15]<<" \n";

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
    // on soustrait à 255 la valeur donnée
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
  if(d_lum<0){
    b=_mm_set1_epi8(abs(d_lum));
    //affiche_xmm_i(b);

    for(int i=0;i<nbPixels;i+=16){
       a=_mm_load_si128((__m128i*)&X[i]);
       // on soustrait la luminance en faisant attention à la saturation
       c=_mm_subs_epu8(a,b);
       _mm_store_si128((__m128i*)&Y[i],c);
     }
  }else{

    b=_mm_set1_epi8(d_lum);
    for(int i=0;i<nbPixels;i+=16){
       a=_mm_load_si128((__m128i*)&X[i]);
       // on ajoute la luminance si elle est positive
       c=_mm_adds_epu8(a,b);
       _mm_store_si128((__m128i*)&Y[i],c);
     }
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

  int nbPixels = h*l;
  __m128i v_seuil = _mm_set1_epi8(seuil);
  __m128i v_moins1 = _mm_set1_epi8(-1);
  __m128i v_donnees,v_resultat;
  for(int i=0;i<nbPixels;i+=16){
    v_donnees=_mm_load_si128((__m128i*)&X[i]);
    // trick pour faire un cmpgt en int donc on va utiliser un ou exclusief sur un masque résultant de la comparaison entre le seuil et la valeur max entre le seuil et la donnée
    v_resultat=_mm_xor_si128(_mm_cmpeq_epi8(_mm_max_epu8(v_seuil,v_donnees),v_seuil),v_moins1);
    _mm_store_si128((__m128i*)&Y[i],v_resultat);
  }
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

  int nbPixels = h*l;
  __m128i v_min = _mm_set1_epi8(255);
  __m128i v_max = _mm_set1_epi8(0);
  __m128i v_donnees;
  min = 255;
  max = 0;

  for(int i=0;i<nbPixels;i+=16){
    v_donnees=_mm_load_si128((__m128i*)&X[i]);
    v_min = _mm_min_epu8(v_min,v_donnees) ;
    v_max= _mm_max_epu8(v_max,v_donnees) ;
  }
  unsigned char tmpMin[16];
  unsigned char tmpMax[16];
  _mm_store_si128((__m128i*)&tmpMin[0],v_min);
  _mm_store_si128((__m128i*)&tmpMax[0],v_max);

// on cherche le min entre tout le vecteur de 16 entiers
  for(int i=0;i<16;i++){
    min = std::min(tmpMin[i],min);
    max = std::max(tmpMax[i],max);
  }


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

  long nbPixels=h*l;

__m128i v_donnees;
__m128i v_min = _mm_set1_epi8(min);
__m128i v_255 = _mm_set1_epi16(255);
__m128 v_div = _mm_set1_ps(max-min);
//__m128i v_255 = _mm_set1_epi16(255);
__m128i v_donneeslo;
__m128i v_donneeshi;
__m128i v_donneeslolo;
__m128i v_donneeshihi;
__m128i v_donneeshilo;
__m128i v_donneeslohi;

__m128 v_donneesloloF;
__m128 v_donneeshihiF;
__m128 v_donneeshiloF;
__m128 v_donneeslohiF;


__m128i v_masque = _mm_set_epi8(0,-1,0,-1,0,-1,0,-1,0,-1,0,-1,0,-1,0,-1);
__m128i v_masque2 = _mm_set_epi16(0,-1,0,-1,0,-1,0,-1);
  for(long i=0;i<nbPixels;i+=16){
    v_donnees=_mm_load_si128((__m128i*)&X[i]);
    v_donnees = _mm_subs_epu8(v_donnees,v_min);

    // séparer en deux vecteurs de forme 0 X 0 X 0 X .... pour pouvoir multiplier sans overflow
    v_donneeslo = _mm_and_si128(v_donnees,v_masque);
    v_donneeshi = _mm_and_si128(_mm_srli_si128(v_donnees,1),v_masque);
    // multiplication  sur 16 bits
    v_donneeslo = _mm_mullo_epi16(v_donneeslo,v_255);
    v_donneeshi = _mm_mullo_epi16(v_donneeshi,v_255);

    // transformer le vecteur d'entiers de 16  bits en 32
    v_donneeshilo = _mm_and_si128(v_donneeshi,v_masque2);
    v_donneeshihi = _mm_and_si128(_mm_srli_si128(v_donneeshi,2),v_masque2);

    v_donneeslolo = _mm_and_si128(v_donneeslo,v_masque2);
    v_donneeslohi = _mm_and_si128(_mm_srli_si128(v_donneeslo,2),v_masque2);

// transformation en float pour pouvoir diviser
    v_donneeshiloF = _mm_cvtepi32_ps(v_donneeshilo);
    v_donneeshihiF = _mm_cvtepi32_ps(v_donneeshihi);
    v_donneeslohiF = _mm_cvtepi32_ps(v_donneeslohi);
    v_donneesloloF = _mm_cvtepi32_ps(v_donneeslolo);

// division sur 32 bits
    v_donneeshiloF = _mm_div_ps(v_donneeshiloF,v_div);
    v_donneeshihiF = _mm_div_ps(v_donneeshihiF,v_div);
    v_donneeslohiF = _mm_div_ps(v_donneeslohiF,v_div);
    v_donneesloloF = _mm_div_ps(v_donneesloloF,v_div);

// reconversion en entier 32 bits pour pouvoir le store
    v_donneeshilo = _mm_cvtps_epi32(v_donneeshiloF);
    v_donneeshihi = _mm_cvtps_epi32(v_donneeshihiF);
    v_donneeslohi = _mm_cvtps_epi32(v_donneeslohiF);
    v_donneeslolo = _mm_cvtps_epi32(v_donneesloloF);

// on récupère les données qui sont sous forme 0 0 0 X 0 0 0 Y .... donc on décale un par un d'un octet supplémentaire afin de se retrouver avec X Y Z W
    v_donnees = _mm_or_si128(_mm_or_si128(_mm_slli_si128(v_donneeshihi,3),_mm_slli_si128(v_donneeslohi,2)),_mm_or_si128(_mm_slli_si128(v_donneeshilo,1),v_donneeslolo));


    _mm_store_si128((__m128i*)&Y[i],v_donnees);


  }
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
  Luminance_sse(-20,(unsigned char*) img->imageData, img->height, img->width, (unsigned char*) image_luminance_sse->imageData);
  cvShowImage( "Affiche_Luminance_sse", image_luminance_sse);



	// cr�ation d'une image pour stocker l'image binaris�e
  IplImage* image_binaire = cvCreateImage (cvGetSize (img), IPL_DEPTH_8U, 1);
  binarise( 120, (unsigned char*) img->imageData, img->height, img->width, (unsigned char*) image_binaire->imageData);
  cvShowImage( "Affiche_binaire", image_binaire);


  	// cr�ation d'une image pour stocker l'image binaris�e
    IplImage* image_binaire_sse = cvCreateImage (cvGetSize (img), IPL_DEPTH_8U, 1);
    binarise_sse( 120, (unsigned char*) img->imageData, img->height, img->width, (unsigned char*) image_binaire_sse->imageData);
    cvShowImage( "Affiche_binaire_sse", image_binaire_sse);




  unsigned char min, max;
  min_et_max( (unsigned char*) img->imageData, img->height, img->width, min, max);
  cout << "min : " << (int)min << "   max : " << (int)max << endl;

  unsigned char min2, max2;
  min_et_max_sse( (unsigned char*) img->imageData, img->height, img->width, min2, max2);
  cout << "min sse : " << (int)min2 << "   max sse : " << (int)max2 << endl;

	// cr�ation d'une image pour stocker l'image de recadrage dynamique
  IplImage* image_ = cvCreateImage (cvGetSize (img), IPL_DEPTH_8U, 1);
  recadrageDynamique( (unsigned char*) img->imageData, img->height, img->width, min, max, (unsigned char*) image_->imageData);
  cvShowImage( "Affiche_", image_);


  	// cr�ation d'une image pour stocker l'image de recadrage dynamique
    IplImage* image_sse = cvCreateImage (cvGetSize (img), IPL_DEPTH_8U, 1);
    recadrageDynamique_sse( (unsigned char*) img->imageData, img->height, img->width, min, max, (unsigned char*) image_sse->imageData);
    cvShowImage( "Affiche_sse", image_sse);




  cvWaitKey();		// attend une touche clavier
  exit (0);
}

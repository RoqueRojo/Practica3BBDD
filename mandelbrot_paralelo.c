/****************** Grupo C2 y B2 respectivamente, Roque Rojo Bacete, Noelia María Granados Carrasco ******************/
/* Conjunto Mandelbrot, M.
 * Sea c un número complejo, n un número natural y {z(n)} una sucesión de números definida por:
 * 		z(0)= 0;
 * 		z(n+1) = z(n)² + c
 * c pertenece a M ssi, para cualquier n, |z(n)| <= 2

 * Como aproximación, tomaremos los valores de c tales que
 *  para cualquier n < MAX_ITER, |z(n)| <= 2

 * Compilación:
 *  gcc mandelbrot_secuencial.c -o mandelbrot_secuencial -lm -fopenmp -O1 -Wno-unused-result
 *
 * Ejecución: ./mandelbrot_secuencial <archivo de entrada>
 *
 * El archivo de entrada debe ajustarse al siguiente formato:
 * 	1.ª línea: <número de imágenes que queremos generar>
 * 	Una línea adicional por cada imagen que queramos generar con uno de los siguientes formatos según
 * 		la imagen que queramos sea rectangular o cuadrada, respectivamente:
 * 		1 <menor abscisa> <mayor abscisa> <menor ordenada> <mayor ordenada> <nombre archivo imagen>
 * 		2 <abscisa centro del cuadrado > <ordenada del centro> <lado del cuadrado> <nombre archivo imagen>
 * 		El nombre del archivo imagen no debe tener más de 15 caracteres ni punto.

 * Cada imagen será de una zona rectangular (o cuadrada) cuya abscisa menor será xmin,
 *  su abscisa mayor será xmax, su ordenada menor será ymin y su ordenada mayor será ymax.
 * Por tanto, la anchura de la zona es xmax-xmin y la altura es ymax-ymin
 * En el archivo de imagen (que tendrá el formato ppm) el número de puntos en las direcciones
 *  horizontal y vertical será proporcional a  xmax-xmin e ymax-ymin, respectivamente. */

#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <math.h>
#include <string.h>

#define MAX_ITER 5000	// Máximo número de términos calculados en cada sucesión
// Cuanto mayor sea MAX_ITER más calidad de imagen y mayor tiempo de ejecución
#define RK 1000			// RK*RK número aproximado de puntos de la imagen

/* La función mandel_val determina si el número x+i·y
 * pertenece (probablemente) o no (con seguridad) al conjunto de Mandelbrot */
int mandel_val(double x, double y, int max_iter){
	int j = 1;
	double temp_real, cr, p_real = x;
	double ci, p_imag = y;

	cr = x*x; ci = y*y;	// cr+ci es el cuadrado del módulo de x + i·y
	while (((cr+ci)<=4) && (j < max_iter)){
		temp_real = cr - ci + x;
		p_imag = 2 * p_real * p_imag + y;
		p_real = temp_real;
		j++;			// z(j) = p_real + i · p_imag
		cr = p_real * p_real;	ci = p_imag * p_imag;
	}

	if((cr+ci) <= 4)	return -1; else return j;
};

int main(int argc, char *argv[])	{

	int n_imag, imag, t_imag, i, j, value, valuer, valueg, valueb;
	float xmin, xmax, ymin, ymax;
	double x,y,l,t;
	char archivoppm[20];
    int nucleos = omp_get_num_procs();
	FILE* filein;
	
    omp_set_num_threads(nucleos);
	// A. Abrir archivo de entrada y leer número de imágenes
	if (argc == 2) {
	   filein = fopen(argv[1], "r");
	   if (filein == NULL) {
		  fprintf(stderr, "No puedo abrir %s\n", argv[1]);
		  fprintf(stderr, "uso: %s <archivo de entrada> \n",argv[0]);
		  exit(0);
	   }
	   fscanf(filein,"%d",&n_imag);
	}	else {
			fprintf(stderr, "uso: %s <archivo de entrada> \n",argv[0]);
			exit(0);
		}

	// B. Obtener una a una las imágenes y los tiempos respectivos
	for (imag = 0; imag <n_imag; imag ++)	{

		t = omp_get_wtime();

		// B1. Leer tipo de imagen (1 = rectangular, 2 = cuadrada)
		fscanf(filein,"%d",&t_imag);
		if (t_imag == 1)	{
			fscanf(filein,"%f %f %f %f %s", &xmin, &xmax, &ymin, &ymax, archivoppm);
		}	else if (t_imag == 2)	{
				fscanf(filein,"%lf %lf %lf %s", &x, &y, &l, archivoppm);
				l = l/2; xmin = x-l; xmax = x+l; ymin = y-l; ymax = y+l;
			}	else {fprintf(stderr, "\n Error en formato del archivo de entrada\n");exit(0);}

		// B2. Calcular, en pixels, la anchura, width, y la altura, height, de la imagen
		double k = RK/sqrt((xmax-xmin)*(ymax-ymin));
		int width = k*(xmax-xmin);	// Hay conversión de double a int
		int height = k*(ymax-ymin);	// Íd.
		// width x height = RK² aproximadamente
		// k = width/(xmax-xmin) = height/(ymax-ymin);

		// B3. Se abre y configura el archivo de imagen ppm
		strcat(archivoppm,".ppm"); //archivoppm = archivoppm + ".ppm";
		FILE *f_imag = fopen(archivoppm, "wb");
		fprintf(f_imag, "P3\n%d %d\n255\n", width, height);
		/* Formato P3: cada pixel se representa por sus valores de rojo, verde y azul (entre 0 y 255);
		 * width x height: dimensión imagen en pixels */

		// B4. Se calcula el punto (x,y) y se determina si incluir (x + i·y) en el conjunto de Mandelbrot
		int pixelsimg [height][width];

        #pragma omp parallel private(x, y, j) shared(pixelsimg, k, xmin, ymax)
        #pragma omp for schedule (dynamic,5)  
		for(i = 0; i < height; i++){
			for(j = 0; j < width; j++) {
				x = xmin + j/k;
				y = ymax - i/k;
                pixelsimg[i][j] = mandel_val(x, y, MAX_ITER);
			}
		}

		// B5. Se determina el color del pixel del número (x + y·i) y se imprime en el archivo
        for (i = 0; i < height; i++) {
            for (j = 0; j < width; j++) {
                if(pixelsimg[i][j] == -1)	// se supone que (x+i·y) pertenece a M --> pixel blanco
					fprintf(f_imag, " 255 255 255 ");
				else
				// se sabe que (x+i·y) no pertenece a M --> color del pixel según número de iteraciones
				{	value = ((float)pixelsimg[i][j]/MAX_ITER) * 16777215;
					// Hacemos que value quede entre 0 y 2²⁴ y extraemos sus componentes rgb
					valueg = value >> 8;	valueb = value % 256;
					valuer = valueg >> 8;   valueg = valueg % 256;
					fprintf(f_imag," %d %d %d ",valuer,valueg,valueb);
				}
            }
        }

		fprintf(f_imag, "\n");	fclose(f_imag);

		// B6. Informamos al usuario
		printf("\n En el archivo %s se ha creado una imagen de la zona explorada", archivoppm);
		printf("\n\n Tiempo empleado %0.2f segundos\n\n",omp_get_wtime()-t);
        
	}	// fin de: for (imag = 0; imag <n; imag ++)
	fclose(filein);
}

/*____Procesador:intel Core i7-4720HQ 2.60GHz_______
______Numero de hilos utilizados:4________________*/

/*____Para el secuencial (sin paralelizar) _________
______en el cuadrado ha tardado: 19.14 segundos_____
______en el rectangulo ha tardado: 8.81 segundos__*/

/*____Para el paralelo sin utilizar schedule________
______en el cuadrado ha tardado: 7.06 segundos______
______en el rectangulo ha tardado: 4.33 segundos__*/

/*____Para el paralelo schedule (static,1)__________
______en el cuadrado ha tardado: 4.96 segundos______
______en el rectangulo ha tardado: 2.52 segundos__*/

/*____Para el paralelo schedule (static,5)__________
______en el cuadrado ha tardado: 4.92 segundos______
______en el rectangulo ha tardado: 2.54 segundos__*/ 

/*____Para el paralelo schedule (static,10)_________
______en el cuadrado ha tardado: 4.92 segundos______
______en el rectangulo ha tardado: 2.56 segundos__*/

/*____Para el paralelo schedule (static,20)_________
______en el cuadrado ha tardado: 5.04 segundos______
______en el rectangulo ha tardado: 2.59 segundos__*/

/*____Para el paralelo schedule (static,40)_________
______en el cuadrado ha tardado: 5.14 segundos______
______en el rectangulo ha tardado: 2.57 segundos__*/

/*____Para el paralelo schedule (static,100)________
______en el cuadrado ha tardado: 4.92 segundos______
______en el rectangulo ha tardado: 2.54 segundos__*/

/*____ESTE ES EL PEOR TIEMPO QUE HEMOS OBTENIDO___*/
/*____Para el paralelo schedule (static,150)________
______en el cuadrado ha tardado: 5.90 segundos______
______en el rectangulo ha tardado: 3.26 segundos__*/

/*____Para el paralelo schedule (dynamic,1)_________
______en el cuadrado ha tardado: 4.84 segundos______
______en el rectangulo ha tardado: 2.46 segundos__*/

/*____ESTE ES EL TIEMPO MÁS OPTIMO OBTENIDO_______*/
/*____Para el paralelo schedule (dynamic,5)_________
______en el cuadrado ha tardado: 4.85 segundos______
______en el rectangulo ha tardado: 2.45 segundos__*/ 

/*____Para el paralelo schedule (dynamic,10)________
______en el cuadrado ha tardado: 4.92 segundos______
______en el rectangulo ha tardado: 2.56 segundos__*/

/*____Para el paralelo schedule (dynamic,20)________
______en el cuadrado ha tardado: 5.02 segundos______
______en el rectangulo ha tardado: 2.46 segundos__*/

/*____Para el paralelo schedule (dynamic,40)________
______en el cuadrado ha tardado: 5.03 segundos______
______en el rectangulo ha tardado: 2.56 segundos__*/

/*____Para el paralelo schedule (dynamic,100)_______
______en el cuadrado ha tardado: 5.96 segundos______
______en el rectangulo ha tardado: 3.20 segundos__*/

/*____Para el paralelo schedule (dynamic,150)_______
______en el cuadrado ha tardado: 5.75 segundos______
______en el rectangulo ha tardado: 3.26 segundos__*/

/*____Conclusión sobre el trabajo: hemos llegado a__
______la conclusión de que en los bordes de la______
______figura es lo que más le cuesta al programa____
______(más tiempo tarda), y que al ser figuras _____
______distintas es más óptimo utilizar el parámetro_
______dynamic, ya que obtenemos mejores tiempos,____
______pero en el caso de analizaramos una figura____
______unica con partes estáticas, el parámetro______
______static sería el que debemos utilizar para ____
______obtener mejores resultados.Gracias a este ____
______trabajo hemos podido aprender la gran ________
______diferencia que puede existir en los tiempos___
______de ejecución al paralelizar un código, al no__
______hacerlo, y la importancia en cuanto al________
______rendimiento.________________________________*/


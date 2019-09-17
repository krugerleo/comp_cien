/*
 * gauss.c
 *
 * Resolve sistemas de equações utilizando o método da eliminação de Gauss.
 * Feito de acordo com <http://rpanta.com/downloads/material/Gauss_01.PDF>.
 *
 * Autor: Pedro Sousa Lacerda, 2010
 *
 */

#include <stdio.h>
#include <stdlib.h>

int main() {

    int i, j, k;


    /* int ico
     *      número de icógnitas presentes no sistema
     *
     * double sis[ico][ico+1]
     *      sistema de equações:
     *          * o multiplicador da icógnita está na mesma coluna em todas as linhas.
     *          * equações com mais icógnitas vem antes das com menos icógnitas.
     *          * possivelmente há outras regras, mas são desconhecidas.
     */


	/*
		// EXEMPLO com uma icógnita
		int ico = 1;
		double sis[1][2] = {{10,5}};
	*/

	/*
		// EXEMPLO com três icógnitas
		int ico = 3;
		double sis[3][4] = {{1, 1, 1, 0},
		                    {2, 1, 1, 1},
		                    {1, 2, 1, 15}};
	*/

    // EXEMPLO com quatro igógnitas
    int ico = 4;
    double sis[4][5] = {{2, 3, 1, 1, 3},
                        {4, 7, 4, 3, 6},
                        {4, 7, 6, 4, 4},
                        {6, 9, 9, 8, 3}};

    /*
     * Eliminação de Gauss
     */
    for (i=1; i<ico; i++)
        for (k=0; k<i; k++)
            if (sis[i][k] != .0)
            {
                // razão para eliminação
                double r = sis[i][k] / sis[k][k];

                // eliminação da linha i, coluna j
                for (j=k; j<ico+1; j++)
                    sis[i][j] -= r*sis[k][j];
            }


    /*
     * Visualizar sistema já eliminado
     */
    puts("\nSistema Eliminado:\n");
    char xnome;
    for (i=0; i<ico; i++)
        for (j=0, xnome='z'-ico+1; j<ico+1; j++, xnome++)
        {
            if          (j < ico-1)     printf("%.2f%c + ", sis[i][j], xnome);
            else if     (j == ico-1)    printf("%.2f%c = ", sis[i][j], xnome);
            else                        printf("%.2f\n", sis[i][j]);
        }
    puts("\n");



    /*
     * Retrosubstituição para encontrar valor das variáveis
     */
    double* x = malloc(sizeof(double)*ico);
    for (i=(ico-1), xnome = 'z'; i>=0; i--, xnome--)
    {
        x[i] = sis[i][ico];
        for (j=ico-1; j>i; j--)
            x[i] -= sis[i][j] * x[j];

        x[i] /= sis[i][i];

        // Impressão dos Resultados
        // printf("%c = %+.2f\n", xnome, x[i]);
    }

    getchar();
    return (EXIT_SUCCESS);
}

/* 
	Trabalho de Cálculo Numérico - 2016
  Nome : Arthur Alexsander Martins Teodoro    Matricula : 0022427
         Saulo Ricardo Dias Fernandes         Matricula : 0021581
  Data : 23/06/2016       
*/

#include <stdio.h>
#include <stdlib.h>

/*Criação das Estruturas usadas*/
typedef struct ponto
{
	float x;
	float y;
} Ponto;

struct points
{
	int numPontos;
	Ponto *pontos;
	float maiorX;
	float maiorY;
	float menorX;
	float menorY;
};

typedef struct points *Points;

/*Declaração de Protótipos das funções*/
Points carregaArquivo(const char *arq);
void destroiTable(Points t);

/*Desenvolvimento das Funções*/
int main(int argc, char const *argv[])
{
	int i;
	Points table = carregaArquivo(argv[1]);
	destroiTable(table);
	return 0;
}

/*Função responsável por ler o arquivo*/
Points carregaArquivo(const char *arq)
{
	int i, quant = 0;
	float x, y, menorx, menory, maiorx, maiory;
	Points aux;
	/*Abre o arquivo para contar quantos pontos existem*/
	FILE *Arq = fopen(arq,"rt");
	if(Arq == NULL)
	{
		return NULL;
	}
	while(fscanf(Arq,"%f %f",&x,&y) != EOF)
	{
		quant++;
	}
	fclose(Arq);
	/*Aloca memória para a estrututa*/
	aux = (Points) malloc(sizeof(struct Point));
	aux->pontos = (Ponto*) malloc(sizeof(struct ponto)*quant);
	/*Insere dados na estrutura*/
	aux->numPontos = quant;
	Arq = fopen(arq,"rt");
	for(i = 0; i < quant; i++)
	{
		fscanf(Arq,"%f %f",&x, &y);
		if(i == 0)
		{
			menorx = x;
			menory = y;
			maiorx = x;
			maiory = y;
		}
		aux->pontos[i].x = x;
		aux->pontos[i].y = y;
		if(x < menorx)
		{
			menorx = x;
		};
		if(y < menory)
		{
			menory = y;
		};
		if(x > maiorx)
		{
			maiorx = x;
		}
		if(y > maiory)
		{
			maiory = y;
		}
	}
	aux->maiorY = maiory;
	aux->maiorX = maiorx;
	aux->menorY = menory;
	aux->menorX = menorx;
	return aux;
}

void destroiTable(Points t)
{
	free(t->pontos);
	free(t);
	t = NULL;
}
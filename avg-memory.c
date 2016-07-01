/*
	Trabalho de Cálculo Numérico - 2016
  Nome : Arthur Alexsander Martins Teodoro    Matricula : 0022427
         Saulo Ricardo Dias Fernandes         Matricula : 0021581
  Data : 23/06/2016
*/

#include <stdio.h>
#include <stdlib.h>
//#include <math.h>

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
void destroiPoints(Points p);
float* CalculaDerivadaSpline(Points p);
float AvaliaSpline(Points p, float* s2, float valor);

/*Desenvolvimento das Funções*/
int main(int argc, char const *argv[])
{
	int i;
	Points table = carregaArquivo(argv[1]);
	float *s2 = CalculaDerivadaSpline(table);

	free(s2);
	destroiPoints(table);
	return 0;
}

/*Requesito 2 - Ler o arquivo com os valores*/
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
	aux = (Points) malloc(sizeof(struct points));
	aux->pontos = (Ponto*) malloc(sizeof(struct ponto)*(quant+1));
	/*Insere dados na estrutura*/
	aux->numPontos = quant;
	Arq = fopen(arq,"rt");
	for(i = 1; i <= quant; i++)
	{
		fscanf(Arq,"%f %f",&x, &y);
		if(i == 1)
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

void destroiPoints(Points p)
{
	free(p->pontos);
	free(p);
	p = NULL;
}

/*Requesito 3 - Calcular as derivadas das Splines*/
float* CalculaDerivadaSpline(Points p)
{
  if(p->numPontos < 3)
  {
  	return NULL;
  }
  /*Declaração das variaveis*/
  int m, i;
  float Ha, Hb, DeltaA, DeltaB, t;
  float *e = (float*) malloc(sizeof(float)*(p->numPontos+1));
  float *d = (float*) malloc(sizeof(float)*(p->numPontos+1));
  float *s2 = (float*) malloc(sizeof(float)*(p->numPontos+1));

  /*Sistema Tridiagonal Simétrico*/
  m = p->numPontos - 2;
  Ha = p->pontos[2].x - p->pontos[1].x;
  DeltaA = (p->pontos[2].y - p->pontos[1].y)/Ha;
  for(i = 1; i <= m; i++)
  {
  	Hb = p->pontos[i+2].x - p->pontos[i+1].x;
  	DeltaB = (p->pontos[i+2].y - p->pontos[i+1].y)/Hb;
  	e[i] = Hb;
  	d[i] = 2*(Ha + Hb);
  	s2[i+1] = 6*(DeltaB - DeltaA);
  	Ha = Hb;
  	DeltaA = DeltaB;
  }

  /*Eliminação de Gauss*/
  for(i = 2; i <= m; i++)
  {
  	t = e[i-1]/d[i-1];
  	d[i] = d[i] - (t*e[i-1]);
  	s2[i+1] = s2[i+1] - (t*s2[i]);
  }

  /*Substituição Retroativa*/
  s2[m+1] = s2[m+1]/d[m];
  for(i = m; i >= 2; --i)
  {
  	s2[i] = (s2[i] - e[i-1] * s2[i+1])/d[i-1];
  }

  s2[1] = 0;
  s2[m+2] = 0;
  free(e);
  free(d);
  return s2;
}

/*Requesito 4 - Avaliar e Gerar o polinomio da Spline*/
float AvaliaSpline(Points p, float* s2, float valor)
{
	if((valor < p->pontos[1].x) || (valor > p->pontos[p->numPontos].x))
	{
		return 0;
	}

	/*Declaração de variaveis*/
	int inf, sup, indice;
	float a,b,c,d,h,resultado;

	/*Busca Binaria*/
	inf = 1;
	sup = p->numPontos;
	while((sup - inf) > 1)
	{
		indice = (inf+sup)/2;
		if(p->pontos[indice].x > valor)
		{
			sup = indice;
		}
		else
		{
			inf = indice;
		}
	}
	
	/*Avaliação de Horner*/
	h = p->pontos[sup].x - p->pontos[inf].x;
	a = (s2[sup] - s2[inf])/(6*h);
	b = s2[inf]*0.5;
	c = ((p->pontos[sup].y - p->pontos[inf].y)/h) - ((s2[sup] + 2*s2[inf])*h/6);
	d = p->pontos[inf].y;
	h = valor - p->pontos[inf].x;
	resultado = ((a*h+b)*h+c)*h+d;
	return resultado;
}

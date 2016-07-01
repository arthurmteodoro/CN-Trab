/*
	Trabalho de Cálculo Numérico - 2016
  Nome : Arthur Alexsander Martins Teodoro    Matricula : 0022427
         Saulo Ricardo Dias Fernandes         Matricula : 0021581
  Data : 23/06/2016
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*Criação das Estruturas usadas*/
typedef struct ponto
{
	double x;
	double y;
} Ponto;

struct points
{
	int numPontos;
	Ponto *pontos;
	double maiorX;
	double maiorY;
	double menorX;
	double menorY;
};

typedef struct points *Points;

/*Declaração de Protótipos das funções*/
Points carregaArquivo(const char *arq);
void destroiPoints(Points p);
double* CalculaDerivadaSpline(Points p);
double AvaliaSpline(Points p, double* s2, double valor);
double geraNum(double min, double max);
double IntegralMonteCarlo(long 	int n, Points p, double* s2);
double TVMI(double a, double b, double integral);
void SaidaTerminal(Points p, double mem, const char *str);
void SaidaR(Points p, double med, const char *str);

/*Desenvolvimento das Funções*/
int main(int argc, char const *argv[])
{
	int i;
	Points table = carregaArquivo(argv[1]);
	double *s2 = CalculaDerivadaSpline(table);
	SaidaTerminal(table,TVMI(table->menorX,table->maiorX,IntegralMonteCarlo(1000,table,s2)),argv[2]);
	free(s2);
	destroiPoints(table);
	return 0;
}

/*Requesito 2 - Ler o arquivo com os valores*/
Points carregaArquivo(const char *arq)
{
	int i, quant = 0;
	double x, y, menorx, menory, maiorx, maiory;
	Points aux;
	/*Abre o arquivo para contar quantos pontos existem*/
	FILE *Arq = fopen(arq,"rt");
	if(Arq == NULL)
	{
		return NULL;
	}
	while(fscanf(Arq,"%lf %lf",&x,&y) != EOF)
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
		fscanf(Arq,"%lf %lf",&x, &y);
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
double* CalculaDerivadaSpline(Points p)
{
  if(p->numPontos < 3)
  {
  	return NULL;
  }
  /*Declaração das variaveis*/
  int m, i;
  double Ha, Hb, DeltaA, DeltaB, t;
  double *e = (double*) malloc(sizeof(double)*(p->numPontos+1));
  double *d = (double*) malloc(sizeof(double)*(p->numPontos+1));
  double *s2 = (double*) malloc(sizeof(double)*(p->numPontos+1));

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
double AvaliaSpline(Points p, double* s2, double valor)
{
	if((valor < p->pontos[1].x) || (valor > p->pontos[p->numPontos].x))
	{
		return 0;
	}

	/*Declaração de variaveis*/
	int inf, sup, indice;
	double a,b,c,d,h,resultado;

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
	resultado = ((a*h+b)*h+c)*h + d;
	return resultado;
}

/*Requesito 5 - Gerador de Números Uniformes*/
double geraNum(double min, double max)
{
	return (rand()/(double)RAND_MAX)*(max-min)+min;
}

/*Requisito 6 - Integral por MonteCarlo*/
double IntegralMonteCarlo(long int n, Points p, double* s2)
{
	/*Declaração das variaveis*/
	long int i;
	double x, y,AreaTotal, Area, xMin, xMax, yMin, yMax, numAbaixo = 0;

	xMin = p->menorX;
	xMax = p->maiorX;	
	yMin = 0;
	yMax = p->maiorY + (p->maiorY*0.2);
	for(i = 1; i < n; i++)
	{
		x = geraNum(xMin,xMax);
		y = geraNum(yMin,yMax);
		if(y <= AvaliaSpline(p,s2,x))
		{
			numAbaixo++;
		}
	}
	AreaTotal = (xMax - xMin)*(yMax - yMin);
	Area = AreaTotal*(numAbaixo/n);
	return Area;
}

/*Requisito 07 - TVMI*/
double TVMI(double a, double b, double integral)
{
	return((1/(b-a))*integral);
}

/*Requisito 08 - Saida Terminal*/
void SaidaTerminal(Points p, double mem, const char *str)
{
	printf("Number of Samples : %d\n", p->numPontos);
	printf("Average Memory Usage : %.3lf Kb \n", mem);
	printf("\nRun 'Rscript %s.r' to generate Avarage Momory Usage Chart\n", str);
}

/*Requisito 09 - Script*/
void SaidaR(Points p, double med, const char *str)
{

}

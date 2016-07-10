/*==================================================================================================*/
/*                             Trabalho de Cálculo Numérico - 2016                                  */
/*                 Nome : Arthur Alexsander Martins Teodoro    Matricula : 0022427                  */
/*                 Nome : Saulo Ricardo Dias Fernandes         Matricula : 0021581                  */
/*                                     Data : 23/06/2016                                            */
/*==================================================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>  

/*==================================================================================================*/
/*===================================Criação das Estruturas Usadas==================================*/
/*==================================================================================================*/  
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

/*==================================================================================================*/
/*===============================Declaração de protótipo das funções================================*/
/*==================================================================================================*/
Points carregaArquivo(const char *arq);
void destroiPoints(Points p);
double* CalculaDerivadaSpline(Points p);
double AvaliaSpline(Points p, double* s2, double valor);
double geraNum(double min, double max);
double IntegralMonteCarlo(long 	int n, Points p, double* s2);
double TVMI(double a, double b, double integral);
void SaidaTerminal(Points p, double mem, const char *str);
void SaidaR(Points p, double med, const char *str, double *s2);

/*==================================================================================================*/
/*==========================================Função Principal========================================*/
/*==================================================================================================*/
int main(int argc, char const *argv[])
{
	srand(time(NULL));
	double integral, tvmi;
	Points table = carregaArquivo(argv[1]);
	double *s2 = CalculaDerivadaSpline(table);
	integral = IntegralMonteCarlo(10000000,table,s2);
	tvmi = TVMI(table->menorX,table->maiorX,integral);
	SaidaTerminal(table,tvmi,argv[2]);
	SaidaR(table,tvmi,argv[2],s2);
	free(s2);
	destroiPoints(table);
	return 0;
}

/*==================================================================================================*/
/*                         REQUESITO 2 - ENTRADA DE DADOS (ARQUIVO)                                 */
/*                    LE O ARQUIVO E COLOCA OS DADOS NA ESTRUTURA CRIADA                            */
/*IN = NOME DO ARQUIVO                                               OUT = PONTEIRO PARA A ESTRUTURA*/
/*==================================================================================================*/
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

/*==================================================================================================*/
/*                                           DESTROI PONTOS                                         */
/*                                      DESALOCA A ESTRUTURA USADA                                  */
/*IN = ESTRUTURA POINTS                                                                   OUT = VOID*/
/*==================================================================================================*/
void destroiPoints(Points p)
{
	free(p->pontos);
	free(p);
	p = NULL;
}

/*==================================================================================================*/
/*                           REQUESITO 3 - CALCULO DAS DERIVIDAS SEGUNDAS                           */
/*                              CALCULA AS DERIVADAS DA SPLINE CÚBICA                               */
/*IN = ESTRUTURA COM OS PONTOS                                   OUT = PONTEIRO PARA VETOR DE DOUBLE*/
/*==================================================================================================*/
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

/*==================================================================================================*/
/*                                    REQUESITO 4 - AVALIA SPLINE                                   */
/*                                  AVALIA A FUNÇÃO INTERPOLADA P(x)                                */
/*IN = 	ESTRUTURA P, VETOR DE DERIVADAS, VALOR A INTERPOLAR                 OUT = Y DE X INTERPOLADO*/
/*==================================================================================================*/
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

/*==================================================================================================*/
/*                         REQUESITO 5 - GERADOR DE NÚMEROS UNIFORMES                               */
/*                    CRIA UM NUMERO ALEATORIA SEGUINDO DISTRIBUIÇÃO UNIFORME                       */
/*IN = VALOR MAXIMO, VALOR MÍNIMO                                              OUT = VALOR ALEATORIO*/
/*==================================================================================================*/
double geraNum(double min, double max)
{
	return (rand()/(double)RAND_MAX)*(max-min)+min;
}

/*==================================================================================================*/
/*                           REQUESITO 6 - INTEGRAL POR MONTE CARLO                                 */
/*                      CALCULA A INTEGRAL NUMÉRICA PELO MÉTODO DE MONTE CARLO                      */
/*IN = NOME DO ARQUIVO                                               OUT = PONTEIRO PARA A ESTRUTURA*/
/*==================================================================================================*/
double IntegralMonteCarlo(long int n, Points p, double* s2)
{
	/*Declaração das variaveis*/
	long int i;
	double x, y,AreaTotal, Area, xMin, xMax, yMin, yMax, numAbaixo = 0;

	xMin = p->menorX;
	xMax = p->maiorX;	
	yMin = 0;
	yMax = p->maiorY + (p->maiorY*0.2);
	for(i = 1; i <= n; i++)
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

/*==================================================================================================*/
/*                                        REQUESITO 7 - TVMI                                        */
/*                 CALCULA O VALOR MÉDIO PELO TEOREMA DO VALOR MÉDIO DAS INTEGRAIS                  */
/*IN = VALOR A, VALOR B, INTEGRAL                                        OUT = VALOR MÉDIO CALCULADO*/
/*==================================================================================================*/
double TVMI(double a, double b, double integral)
{
	return ((1/(b-a))*integral);
}

/*==================================================================================================*/
/*                                REQUESITO 8 - SAIDA NO TERMINAL                                   */
/*                                   REALIZA A SAIDA NO TERMINAL                                    */
/*IN = PONTOS P, VALOR MÉDIO, ARQUIVO DE SAÍDA                                            OUT = VOID*/
/*==================================================================================================*/
void SaidaTerminal(Points p, double mem, const char *str)
{
	printf("Number of Samples : %d\n", p->numPontos);
	printf("Average Memory Usage : %.3lf Kb \n", mem);
	printf("\nRun 'Rscript %s.r' to generate Avarage Momory Usage Chart\n", str);
}

/*==================================================================================================*/
/*                            REQUESITO 9 - SAIDA DE DADOS R SCRIPT                                 */
/*                                    REALIZA A SAIDA EM R SCRIPT                                   */
/*IN = PONTOS P, VALOR MÉDIO, ARQUIVO DE SAIDA, PONTEIRO DERIVADAS                        OUT = VOID*/
/*==================================================================================================*/
void SaidaR(Points p, double med, const char *str, double *s2)
{
	/*Concatenação das strings*/
	char *Rarq = (char*) malloc(sizeof(char)*(strlen(str)+2));
	char *PNGarq = (char*) malloc(sizeof(char)*(strlen(str)+4));
	strcpy(Rarq,str);
	strcpy(PNGarq,str);
	strcat(Rarq,".r");
	strcat(PNGarq,".png");
	/*Criação das variáveis*/
	FILE *arq;
	int i;
	double pnt;
	char aspas = '"';
	arq = fopen(Rarq, "wt");
	fprintf(arq, "#\n");
	fprintf(arq, "# Generated automatically by %cavg-memory%c application\n",aspas,aspas);
	fprintf(arq, "#\n\n");
	fprintf(arq, "# Original points (x cordinates)\n");
	/*Escreve o vetor X original*/
	fprintf(arq, "xorig <- c(\n");
	for(i=1; i <= p->numPontos; i++)
	{
		if(i != p->numPontos)
		{
			fprintf(arq, "\t%lf,\n", p->pontos[i].x);
		}
		else
		{
			fprintf(arq, "\t%lf\n", p->pontos[i].x);
		}
	}
	fprintf(arq, ");\n\n");
	/*Escreve o Y Original*/
	fprintf(arq, "# Original points (y cordinates)\n");
	fprintf(arq, "yorig <- c(\n");
	for(i=1; i <= p->numPontos; i++)
	{
		if(i != p->numPontos)
		{
			fprintf(arq, "\t%lf,\n", p->pontos[i].y);
		}
		else
		{
			fprintf(arq, "\t%lf\n", p->pontos[i].y);
		}
	}
	fprintf(arq, ");\n\n");
	fprintf(arq, "# Spline points (x cordinates, sampling interval = 0.01)\n");
	/*Imprime o x de 0.01 em 0.01*/
	fprintf(arq, "xspl <- c(\n");
	for(pnt=p->menorX; pnt <= (p->maiorX-0.01); pnt = pnt + 0.01)
	{
		fprintf(arq, "\t%lf,\n", pnt);
	}
  fprintf(arq, "\t%lf\n", p->maiorX);
  fprintf(arq, ");\n\n");
  fprintf(arq, "# Splines points (y cordinates, sampling interval = 0.01)\n");
  /*Imprime o vetor de Y Interpolado*/
  fprintf(arq, "yspl <- c(\n");
  for(pnt=p->menorX; pnt <= (p->maiorX-0.01); pnt = pnt + 0.01)
	{
		fprintf(arq, "\t%lf,\n", AvaliaSpline(p,s2,pnt));
	}
	fprintf(arq, "\t%lf\n", AvaliaSpline(p,s2,p->maiorX));
	fprintf(arq, ");\n\n");
	fprintf(arq, "# Average Memory Usage\n");
	fprintf(arq, "AvgMemory <- %lf;\n\n", med);
	fprintf(arq, "# Plot the values in .png file\n");
	fprintf(arq, "png(file=%c%s%c, width=1200);\n",aspas,PNGarq,aspas);
	fprintf(arq, "title <- paste(%cAVG Memory Usage: %lf Kb (%d Samples)%c);\n", aspas,med,p->numPontos,aspas);
	fprintf(arq, "plot(xspl, yspl, type=%cl%c, col=%cblue%c, main=title, xlab=%cSamples%c, ylab=%cMem. Usage%c, lwd=3);\n", aspas,aspas,aspas,aspas,aspas,aspas,aspas,aspas);
	fprintf(arq, "points(xorig, yorig, pch=19, col=%cred%c);\n", aspas,aspas);
	fprintf(arq, "lines( c(min(xorig), max(xorig)), c(AvgMemory, AvgMemory), col=%cblack%c, lty=12, lwd=3);\n", aspas,aspas);
	fprintf(arq, "dev.off();");
	fclose(arq);
	free(Rarq);
	free(PNGarq);
}
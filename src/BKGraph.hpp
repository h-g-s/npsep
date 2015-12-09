#ifndef BKGRAPH_HPP_INCLUDED
#define BKGRAPH_HPP_INCLUDED

#include "BKVertex.hpp"
#include <string>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <set>
#include <ctime>


typedef struct lista{
    int vertice;
    struct lista* prox;
}Lista;

typedef struct clique{
    unsigned long int *vetorVertices;
    Lista* listaVerticesClique;
    int peso;
}Clique;



extern "C"
{
    #include "cgraph.h"
    #include "clique.h"
    #include "memory.h"
}

class BKGraph
{
public:
    BKGraph(const char*);
    BKGraph(const CGraph*);
    virtual ~BKGraph();
    void insereOrdenado (Clique* P, int vertice);
    Clique* criarClique(int tamanho);
    void removerVertice(Clique* P, int vertice);
    void subtrairPeso(Clique* c, int peso);
    void adicionarPeso(Clique* c, int peso);
    void liberarClique(Clique* c);
    void liberaRec (Lista* l);
    void copiaClique1(Clique* C, Clique* Caux);
    int escolherVerticeMaiorGrauModificado(Clique* P);
    void adicionarVerticeClique1(Clique* P, int vertice, unsigned mask[]);
    int BronKerbosch(Clique *C, Clique *P, Clique *S, int, double, clock_t, unsigned mask[], int **bit);
    void intersecaoOrdenado(Clique* P, int v, Clique* vetorAux, int cont, int P_sem_vizinho[]);
    int execute(int, double);
    void intersecao1(Clique* S, Clique* Saux, int quantidadeVertices, int** bit, int vertice, unsigned mask[]);
    void excluirVizinhos(int P_sem_vizinhos_U[], Clique* P, int u);
    int busca(int cont, int P_sem_vizinho[], int vertice);
    CliqueSet* getCliqueSet();

private:
    const CGraph *cgraph;
    std::vector<BKVertex> vertices;
    CliqueSet *clqSet;
};

struct sortByMDegree
{
    bool operator()(const BKVertex &x, const BKVertex &y)
    {
        if(x.getMDegree() != y.getMDegree())
            return x.getMDegree() > y.getMDegree();

        else if(x.getDegree() != y.getDegree())
            return x.getDegree() > y.getDegree();

      return x.getNumber() < y.getNumber();
    }
};



#endif // BKGRAPH_HPP_INCLUDED

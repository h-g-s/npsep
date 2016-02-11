#include <cassert>
#include "BKGraph.hpp"
#define INT_SIZE (8*sizeof(int))


using namespace std;

BKGraph::BKGraph(const char* fileName)
{
    cgraph = cgraph_load(fileName);
    const int* pesos;
    BKVertex aux;

    int nVertices = cgraph_size(cgraph);

    pesos = cgraph_get_node_weights(cgraph);

    aux.setNumber(0);
    aux.setWeight(0);
    aux.setDegree(0);
    aux.setMDegree(0);

    vertices.push_back(aux);

    for(int i = 0; i < nVertices; i++)
    {
        aux.setNumber(i + 1);

        aux.setWeight(pesos[i]);

        int neighs[10 * cgraph_degree(cgraph, i)];

        int realDegree = cgraph_get_all_conflicting( cgraph, i, neighs, 10 * nVertices );

        aux.setDegree(realDegree);

        aux.setMDegree(realDegree);

        for(int j = 0; j < aux.getDegree(); j++)
            aux.insertConflict(neighs[j] + 1);

        vertices.push_back(aux);

        aux.clearAllConflicts();
    }

    for(int i = 1; i <= nVertices; i++)
    {
        int mdegree = vertices[i].getDegree();

        for(int j = 0; j < vertices[i].getNumberOfConflicts(); j++)
        {
            int position = vertices[i].getNodeConflict(j);
            assert( position != i && position >= 0 && position < (int)vertices.size() );
            mdegree += vertices[position].getDegree();
        }
        vertices[i].setMDegree(mdegree);
    }

    clqSet = clq_set_create();
}

BKGraph::BKGraph(const CGraph *cgraph)
{
    this->cgraph = cgraph;
    const int* pesos;
    BKVertex aux;

    int nVertices = cgraph_size(cgraph);

    pesos = cgraph_get_node_weights(cgraph);

    aux.setNumber(0);
    aux.setWeight(0);
    aux.setDegree(0);
    aux.setMDegree(0);

    vertices.push_back(aux);

    for(int i = 0; i < nVertices; i++)
    {
        aux.setNumber(i + 1);

        aux.setWeight(pesos[i]);

        int neighs[10 * cgraph_degree(cgraph, i)];

        int realDegree = cgraph_get_all_conflicting( cgraph, i, neighs, 10 * nVertices );

        aux.setDegree(realDegree);

        aux.setMDegree(realDegree);

        for(int j = 0; j < aux.getDegree(); j++)
            aux.insertConflict(neighs[j] + 1);

        vertices.push_back(aux);

        aux.clearAllConflicts();
    }

    for(int i = 1; i <= nVertices; i++)
    {
        int mdegree = vertices[i].getDegree();

        for(int j = 0; j < vertices[i].getNumberOfConflicts(); j++)
        {
            int position = vertices[i].getNodeConflict(j);
            assert( position != i && position >= 0 && position < (int)vertices.size() );
            mdegree += vertices[position].getDegree();
        }
        vertices[i].setMDegree(mdegree);
    }

    clqSet = clq_set_create();
}

BKGraph::~BKGraph()
{
    clq_set_free( &(clqSet) );
}

int retornarPeso(Clique* c)
{
    return c->peso;
}

int vazio1(Clique* S, int quantidadeVertices)
{
    unsigned int i;
    for(i = 0; i < quantidadeVertices/INT_SIZE+1; ++i)
    {
        if(S->vetorVertices[i] != 0)
            return 1;
    }
    return 0;
}

int BKGraph::escolherVerticeMaiorGrauModificado(Clique* P)
{
    Lista* p = P->listaVerticesClique;
    int valorMaiorGrauModificado = 0;
    int posicaoMaiorGrauModificado = 0;
    while(p != NULL)
    {
        if(vertices[(p->vertice)+1].getMDegree() > valorMaiorGrauModificado)
        {

            posicaoMaiorGrauModificado = p->vertice;
            valorMaiorGrauModificado = vertices[(p->vertice)+1].getMDegree();
        }
        p = p->prox;
    }
    return posicaoMaiorGrauModificado;
}

void BKGraph::liberarClique(Clique* c)
{
    free(c->vetorVertices);
    liberaRec(c->listaVerticesClique);
}
void BKGraph::liberaRec (Lista* l)
{
     if (l != NULL)
     {
        liberaRec (l->prox);
        free (l);
     }
}

void BKGraph::excluirVizinhos( int P_sem_vizinhos_U[], Clique* P, int u)
{
    Lista* p = P->listaVerticesClique;
    int contador = 1;
    while(p != NULL)
    {
        if(cgraph_conflicting_nodes(cgraph, u, p->vertice) == 0)
        {
            P_sem_vizinhos_U[contador] = p->vertice;
            contador++;
        }
        p = p->prox;
    }
    P_sem_vizinhos_U[0] = contador;
}

void BKGraph::copiaClique1(Clique* C, Clique* Caux)
{
    unsigned int i;
    for(i = 0; i < (vertices.size()-1)/INT_SIZE+1; ++i)
         Caux->vetorVertices[i] = C->vetorVertices[i];
    adicionarPeso(Caux, retornarPeso(C));
}

void BKGraph::adicionarVerticeClique1(Clique* P, int vertice, unsigned mask[])
{
    P->vetorVertices[vertice/INT_SIZE] |= mask[vertice%INT_SIZE];
}

void BKGraph::removerVertice(Clique* P, int vertice)
{
   P->vetorVertices[vertice] = 0;
   Lista* ant = NULL;
   Lista* p = P->listaVerticesClique;
   while (p != NULL && p->vertice != vertice)
   {
         ant = p;
         p = p->prox;
   }
   if(p == NULL)
      return;
   if (ant == NULL)
      P->listaVerticesClique = p->prox;      // Elemento foi encontrado na 1a posição da lista
   else
      ant->prox = p->prox; // Elemento encontrado no meio da lista
   free (p);
}

int BKGraph::busca(int cont, int P_sem_vizinho[], int vertice)
{
    int i;
    for(i = 1; i < cont; ++i)
    {
        if(P_sem_vizinho[i] == vertice)
            return 1;
    }
    return 0;
}

void BKGraph::intersecaoOrdenado( Clique* P, int v, Clique* vetorAux, int cont, int P_sem_vizinho[])
{

     Lista* p;
     for (p = P->listaVerticesClique; p != NULL; p = p->prox)
     {
         if(cgraph_conflicting_nodes(cgraph, v, p->vertice) == 1 && (busca(cont, P_sem_vizinho, v) == 0))
         {
             insereOrdenado(vetorAux, p->vertice);
             vetorAux->vetorVertices[p->vertice] = 1;
             adicionarPeso(vetorAux, vertices[p->vertice+1].getWeight());
         }
     }
}

void BKGraph::intersecao1(Clique* S, Clique* Saux, int quantidadeVertices, int** bit, int vertice, unsigned mask[])
{
    unsigned int i;
    for(i = 0; i < quantidadeVertices/INT_SIZE+1; ++i)
        Saux->vetorVertices[i] = S->vetorVertices[i] & bit[vertice][i];

}

void BKGraph::subtrairPeso(Clique* c, int peso)
{
    c->peso -= peso;
}

int BKGraph::BronKerbosch(Clique *C, Clique*P, Clique *S, int minWeight, double timeLimit, clock_t init, unsigned mask[], int **bit)
{
    clock_t end = clock();
    double sec = ((double)(end - init)) / ((double)CLOCKS_PER_SEC);
    if(sec >= timeLimit)
        return 1;
    if((retornarPeso(P) == 0) && (vazio1(S, (vertices.size()-1)) == 0))
    {
        if(retornarPeso(C) >= minWeight)
        {
          //  printf("S = %d\n", S->vetorVertices[1]);

            //printf("%d\n", cliquesa);
             // printf("\n[%d] ", retornarPeso(C));
              int contador;
              unsigned int valor, t;
              int nodes[vertices.size()-1];
              int cont = 0;
             // printf("\n");
              for(t = 0; t < ((vertices.size()-1)/INT_SIZE + 1); ++t)
              {
                  contador = (INT_SIZE * t)+1;
                  valor =  C->vetorVertices[t];
                  while(valor > 1)
                  {
                      if(valor % 2 == 1)
                      {
                       //  printf("%d ", contador);
                          nodes[cont] = contador - 1;
                          cont++;
                      }
                      valor = valor/2;
                      contador++;
                  }
                  if(valor == 1)
                  {
                    //  printf("%d ", contador);
                      nodes[cont] = contador - 1;
                      cont++;
                  }
              }
         //     printf("\n");
              clq_set_add(clqSet, cont, nodes, retornarPeso(C));
        }
    }

    if(retornarPeso(C) + retornarPeso(P) >= minWeight)
    {
                int u = escolherVerticeMaiorGrauModificado(P);
                int P_sem_vizinhos_U[vertices.size()];
                excluirVizinhos(P_sem_vizinhos_U, P, u);
                int cont;
                Clique* Paux;
                Clique* Saux;
                Clique* Caux;
               for(cont = 1; cont < P_sem_vizinhos_U[0]; ++cont)
               {
                    int v = P_sem_vizinhos_U[P_sem_vizinhos_U[0]-cont];
                    Paux = criarClique(vertices.size());
                    Saux = criarClique((vertices.size()-1)/INT_SIZE + 1);
                    Caux = criarClique((vertices.size()-1)/INT_SIZE + 1);
                    intersecaoOrdenado(P, v, Paux, (P_sem_vizinhos_U[0]-cont), P_sem_vizinhos_U);
                    //intersecao(P, matriz, v, Paux);
                    intersecao1(S, Saux, (vertices.size()-1), bit, v, mask);
                    copiaClique1(C, Caux);
                    adicionarVerticeClique1(Caux, v, mask);
                    adicionarPeso(Caux, vertices[v+1].getWeight());
                    BronKerbosch(Caux, Paux, Saux, minWeight, timeLimit, init, mask, bit);
                    liberarClique(Paux);
                    liberarClique(Saux);
                    liberarClique(Caux);
                    free(Paux);
                    free(Saux);
                    free(Caux);
                    subtrairPeso(P, vertices[v+1].getWeight());
                    removerVertice(P, v);
                    adicionarVerticeClique1(S, v, mask);
                }
    }

    return 0;
}


Clique* BKGraph::criarClique(int tamanho)
{
    Clique* c = (Clique*) malloc(sizeof(Clique));
    c->vetorVertices = (unsigned long int*) calloc (tamanho, sizeof(unsigned long int));
    c->listaVerticesClique = NULL;
    c->peso = 0;
    return c;
}

void BKGraph::adicionarPeso(Clique* c, int peso)
{
    c->peso += peso;
}

void BKGraph::insereOrdenado (Clique* P, int vertice)
{
       Lista* novo = (Lista*) malloc(sizeof(Lista));
       if(P->listaVerticesClique == NULL)
       {
           novo->vertice = vertice;
           novo->prox = P->listaVerticesClique;
           P->listaVerticesClique = novo;
           return;
       }
       Lista* ant = NULL;
       Lista* p = P->listaVerticesClique;
       while(p!= NULL && vertices[p->vertice+1].getMDegree() > vertices[vertice+1].getMDegree())
       {
            ant = p;
            p = p->prox;
       }
       if(p == NULL)
       {
           novo->vertice = vertice;
           novo->prox = NULL;
           ant->prox = novo;
       }
       if(ant == NULL)
       {
           novo->vertice = vertice;
           novo->prox = P->listaVerticesClique;
           P->listaVerticesClique = novo;
           return;
       }
       else
       {
           novo->vertice = vertice;
           novo->prox = p;
           ant->prox = novo;
       }
}


int BKGraph::execute(int minWeight, double timeLimit)
{
    clock_t init;
    unsigned int i;
    Clique* C = criarClique((vertices.size()-1)/INT_SIZE + 1);
    Clique* P = criarClique(vertices.size()-1);
    Clique* S = criarClique((vertices.size()-1)/INT_SIZE + 1);
    unsigned mask[INT_SIZE];
    mask[0] = 1;
    for(unsigned h=1;h<INT_SIZE;h++)
        mask[h] = mask[h-1]<<1;
    int **bit;
    bit = (int**) calloc ((vertices.size()-1), sizeof(int*));
    for(i = 0; i < (vertices.size()-1); ++i)
        bit[i] = (int*) calloc((vertices.size()-1)/INT_SIZE + 1, sizeof(int));
    for(unsigned int i = 0; i < (vertices.size()-1); i++)
    {
        insereOrdenado(P, i);
        P->vetorVertices[i] = 1;
        adicionarPeso(P, vertices[i+1].getWeight());
    }
    unsigned int v,y;
    for(v = 0; v < (vertices.size()-1); ++v)
    {
        for(y = (v+1); y < (vertices.size()-1); ++y)
        {
            if(cgraph_conflicting_nodes(cgraph, v, y) == 1)
            {
                 bit[y][v/INT_SIZE] |= mask[v%INT_SIZE];
            bit[v][y/INT_SIZE] |= mask[y%INT_SIZE];

            }
        }
    }

    init = clock();

    int stat = BronKerbosch(C, P, S, minWeight, timeLimit, init, mask, bit);

    return stat;


}

CliqueSet* BKGraph::getCliqueSet() { return clqSet; }

#include <cassert>
#include <climits>
#include "BKGraph.hpp"
#define INT_SIZE (8*sizeof(int))


using namespace std;

BKGraph::BKGraph(const CGraph *cgraph)
{
    this->cgraph = cgraph;
    BKVertex aux;
    const int* pesos = cgraph_get_node_weights(cgraph);
    int nVertices = cgraph_size(cgraph);
    int *neighs = new int[nVertices*2];

    for(int i = 0; i < nVertices; i++)
    {
        aux.setNumber(i);
        aux.setWeight(pesos[i]);
        int realDegree = cgraph_degree(cgraph, i);
        int check = cgraph_get_all_conflicting( cgraph, i, neighs, nVertices*2 );
        assert(check == realDegree);
        int mdegree = realDegree;
        for(int j = 0; j < realDegree; j++)
        {
            assert(neighs[j] != i && neighs[j] >= 0 && neighs[j] < nVertices);
            mdegree += cgraph_degree(cgraph, neighs[j]);
        }

        aux.setDegree(realDegree);
        aux.setMDegree(mdegree);
        vertices.push_back(aux);
    }

    clqSet = clq_set_create();
    status = 0;
    maxWeight = 0;
    it = 0;
    maxIt = (INT_MAX/100);
    minWeight = 0;

    delete[] neighs;
}

BKGraph::~BKGraph()
{
    clq_set_free( &(clqSet) );
}

int retornarPeso(Clique* c)
{
    return c->peso;
}

int vazio(Clique* P)
{
    if(P->listaVerticesClique == NULL)
        return 1;
    else
        return 0;
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
        if(vertices[p->vertice].getMDegree() > valorMaiorGrauModificado)
        {

            posicaoMaiorGrauModificado = p->vertice;
            valorMaiorGrauModificado = vertices[p->vertice].getMDegree();
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
    for(i = 0; i < vertices.size()/INT_SIZE+1; ++i)
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
            adicionarPeso(vetorAux, vertices[p->vertice].getWeight());
//             adicionarVerticeCliqueOrdenado(vetorAux, p->vertice, grauModificado);
        }
    }
}

void BKGraph::intersecao1(Clique* S, Clique* Saux, int** bit, int vertice, unsigned mask[])
{
    unsigned int i;
    for(i = 0; i < vertices.size()/INT_SIZE+1; ++i)
        Saux->vetorVertices[i] = S->vetorVertices[i] & bit[vertice][i];

}

void BKGraph::subtrairPeso(Clique* c, int peso)
{
    c->peso -= peso;
}

void BKGraph::algoritmoBronKerbosch(Clique *C, Clique*P, Clique *S, unsigned mask[], int **bit)
{
    if(it > maxIt)
        return;
    it++;

    if((vazio(P) == 1) && (vazio1(S, (int)vertices.size()) == 0))
    {
        if(retornarPeso(C) >= minWeight)
        {
            int contador;
            unsigned int valor, t;
            int nodes[vertices.size()];
            int cont = 0;

            for(t = 0; t < vertices.size()/INT_SIZE + 1; ++t)
            {
                contador = INT_SIZE * t;
                valor =  C->vetorVertices[t];
                while(valor > 1)
                {
                    if(valor % 2 == 1)
                    {
                        nodes[cont] = contador;
//                        printf("%d ", contador);
                        cont++;
                    }
                    valor = valor/2;
                    contador++;
                }
                if(valor == 1)
                {
                    nodes[cont] = contador;
//                    printf("%d ", contador);
                    cont++;
                }
            }
//            printf("\n");
            clq_set_add(clqSet, cont, nodes, retornarPeso(C));

            if(retornarPeso(C) > maxWeight)
                maxWeight = retornarPeso(C);
        }
    }

    if(retornarPeso(C) + retornarPeso(P) >= minWeight)
    {
        int u = escolherVerticeMaiorGrauModificado(P);
        int P_sem_vizinhos_U[vertices.size()+1];
        excluirVizinhos(P_sem_vizinhos_U, P, u);
        int cont;
        Clique* Paux;
        Clique* Saux;
        Clique* Caux;
        for(cont = 1; cont < P_sem_vizinhos_U[0]; ++cont)
        {
            int v = P_sem_vizinhos_U[P_sem_vizinhos_U[0]-cont];
            Paux = criarClique((int)vertices.size()+1);
            Saux = criarClique((int)vertices.size()/INT_SIZE + 1);
            Caux = criarClique((int)vertices.size()/INT_SIZE + 1);
            intersecaoOrdenado(P, v, Paux, (P_sem_vizinhos_U[0]-cont), P_sem_vizinhos_U);
            //intersecao(P, matriz, v, Paux);
            intersecao1(S, Saux, bit, v, mask);
            copiaClique1(C, Caux);
            adicionarVerticeClique1(Caux, v, mask);
            adicionarPeso(Caux, vertices[v].getWeight());
            algoritmoBronKerbosch(Caux, Paux, Saux, mask, bit);
            liberarClique(Paux);
            liberarClique(Saux);
            liberarClique(Caux);
            free(Paux);
            free(Saux);
            free(Caux);
            subtrairPeso(P, vertices[v].getWeight());
            removerVertice(P, v);
            adicionarVerticeClique1(S, v, mask);
        }
    }
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
    while(p!= NULL && vertices[p->vertice].getMDegree() > vertices[vertice].getMDegree())
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


int BKGraph::execute()
{
    unsigned int i;
    Clique* C = criarClique((int)vertices.size()/INT_SIZE + 1);
    Clique* P = criarClique((int)vertices.size());
    Clique* S = criarClique((int)vertices.size()/INT_SIZE + 1);
    unsigned mask[INT_SIZE];
    mask[0] = 1;
    for(unsigned h=1; h<INT_SIZE; h++)
        mask[h] = mask[h-1]<<1;
    int **bit;
    bit = (int**) calloc (vertices.size(), sizeof(int*));
    for(i = 0; i < vertices.size(); ++i)
        bit[i] = (int*) calloc(vertices.size()/INT_SIZE + 1, sizeof(int));
    for(unsigned int i = 0; i < vertices.size(); i++)
    {
        insereOrdenado(P, i);
        P->vetorVertices[i] = 1;
        adicionarPeso(P, vertices[i].getWeight());
    }
    unsigned int v,y;
    for(v = 0; v < vertices.size(); ++v)
    {
        for(y = (v+1); y < vertices.size(); ++y)
        {
            if(cgraph_conflicting_nodes(cgraph, v, y) == 1)
            {
                bit[y][v/INT_SIZE] |= mask[v%INT_SIZE];
                bit[v][y/INT_SIZE] |= mask[y%INT_SIZE];

            }
        }
    }
    
    algoritmoBronKerbosch(C, P, S, mask, bit);

    for(i = 0; i < vertices.size(); ++i)
        free(bit[i]);
    free(bit);
    liberarClique(P);
    liberarClique(S);
    liberarClique(C);
    free(P);
    free(S);
    free(C);

    return (it > maxIt);
}

CliqueSet* BKGraph::getCliqueSet()
{
    return clqSet;
}

int BKGraph::getMaxWeight()
{
    return maxWeight;
}

void BKGraph::setMinWeight(int _minWeight)
{
    assert(_minWeight >= 0);
    minWeight = _minWeight;
}

void BKGraph::setMaxIt(int _maxIt)
{
    assert(_maxIt > 0);
    maxIt = _maxIt;
}

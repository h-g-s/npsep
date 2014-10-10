#include <cassert>
#include "BKGraph.hpp"

using namespace std;

BKGraph::BKGraph(const string& fileName)
{
    char ch,SMP,edge;
    char line[100];
    int vert,edges,weight, max = 0;
    FILE* file;

    file = fopen(fileName.c_str(), "r");

    instance = fileName;

    if(file == NULL)
    {
        fprintf( stderr, "Could not open file " );
        exit( EXIT_FAILURE );
    }

    int status;
    status = fscanf(file,"%c",&ch);
    assert( status>0 );


    while(ch!='e')
    {
        if(ch == 'p')
        {
            char* strRead = fgets(line, 100, file) ;
            assert( strRead );
            int statusRead  = sscanf( line + 1, "%c %s %d %d", &SMP, &edge, &nVertices, &nEdges);
            assert( statusRead == 4 );
        }
        else if(ch == 'c')
        {
            char* strRead = fgets(line, 100, file);
            assert( strRead );
        }
        int chRead = fscanf(file,"%c",&ch);
        assert( chRead );
    }

    vertices.resize(nVertices + 1);

    for(int i = 0; i <= nVertices; i++)
    {
        vertices[i].setDegree(0);
        vertices[i].setNumber(i);
    }

    for(int i = 0; i < nEdges; i++)
    {
        char *strRead = fgets(line, 100, file);
        assert( strRead );
        sscanf( line+1, "%d %d", &vert, &edges);

        assert( vert <= nVertices );
        assert( edges <= nVertices );

        vertices[vert].insertConflict(edges);
        vertices[edges].insertConflict(vert);
        vertices[vert].setDegree(vertices[vert].getDegree() + 1);
        vertices[edges].setDegree(vertices[edges].getDegree() + 1);
    }

    for(int i = 0; i < nVertices; i++)
    {
        char *strRead = fgets(line, 100, file);
        assert( strRead );
        sscanf( line + 1, "%d %d", &vert, &weight);
        vertices[vert].setWeight(weight);
        vertices[vert].setMDegree(vertices[vert].getDegree());
    }

    for(int i = 1; i <= nVertices; i++)
    {
        if( max < vertices[i].getDegree() )
            max = vertices[i].getDegree();
    }

    maxDegree = max;

    for(int i = 1; i <= nVertices;i++)
    {
        for(int j = 0; j < vertices[i].getNumberOfConflicts(); j++)
        {
            int position = vertices[i].getNodeConflict(j);
            vertices[i].setMDegree( vertices[i].getMDegree() + vertices[position].getDegree() );
        }
    }

    fclose(file);
}

BKGraph::BKGraph(const CGraph* cgraph)
{
    const int* pesos;
    BKVertex aux;

    nVertices = cgraph_size(cgraph);

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

        aux.setMDegree(aux.getDegree());

        for(int j = 0; j < aux.getDegree(); j++)
            aux.insertConflict(neighs[j] + 1);

        vertices.push_back(aux);

        aux.clearAllConflicts();
    }

    for(int i = 1; i <= nVertices; i++)
    {
        int mdegree = vertices[i].getMDegree();

        for(int j = 0; j < vertices[i].getNumberOfConflicts(); j++)
        {
            int position = vertices[i].getNodeConflict(j);
            assert( position < (int)vertices.size() );
            mdegree += vertices[position].getDegree();
            vertices[i].setMDegree(vertices[i].getMDegree() + vertices[position].getDegree());
        }
    }

    //free((int*)pesos);
}

BKGraph::~BKGraph(){}

vector<BKVertex> BKGraph::getVertices() const { return vertices; }

int BKGraph::getNVertices() const { return nVertices; }

int BKGraph::getNEdges() const { return nEdges; }

int BKGraph::getMaxDegree() const { return maxDegree; }

string BKGraph::getInstance() const { return instance; }

BKVertex BKGraph::getVertex(int v) const { return vertices[v]; }

void BKGraph::setNVertices(int nv){ nVertices = nv; }

void BKGraph::setNEdges(int ne){ nEdges = ne; }

void BKGraph::setMaxDegree(int md) { maxDegree = md; }

int BKGraph::weightCompute(const set<int>& clique)
{
    int weight = 0;

    for(set<int>::iterator it = clique.begin(); it != clique.end(); ++it)
        weight += vertices[*it].getWeight();

    return weight;
}

int BKGraph::weightEstimate(const set<int>& P)
{
    int weight = 0;

    for(set<int>::iterator it = P.begin(); it != P.end(); ++it)
        weight += vertices[*it].getWeight();

    return weight;
}

int BKGraph::writeSolutions(const char* filename)
{
    set<set<int> >::iterator extIterator;
    set<int>::iterator intIterator;
    FILE* arq;
    int totalWeight = 0;

    if((arq = fopen(filename, "w")) == NULL)
        perror("Couldn't create file!\n");

    for(extIterator = cliques.begin(); extIterator != cliques.end(); ++extIterator)
    {
        fprintf(arq, "[%d] ", weightCompute(*extIterator));
        totalWeight += weightCompute(*extIterator);
        for(intIterator = extIterator->begin(); intIterator != extIterator->end(); ++intIterator)
            fprintf(arq, "%d ", *intIterator);
        fprintf(arq,"\n");
    }

    return totalWeight;
}

int BKGraph::BronKerbosch(set<int> C, set<int> P, set<int> S, int minWeight, double timeLimit, clock_t init)
{
    clock_t end = clock();
    double sec = ((double)(end - init)) / ((double)CLOCKS_PER_SEC);
    if(sec >= timeLimit)
        return 1;

    if( P.size() == 0 && S.size() == 0 )
    {
        if(weightCompute(C) >= minWeight)
            cliques.insert(C);
    }

    if(weightCompute(C) + weightEstimate(P) >= minWeight)
    {
        set<BKVertex, sortByMDegree> pivot;

        for(set<int>::iterator a = P.begin(); a != P.end();++a)
            pivot.insert(vertices[*a]);

        if(pivot.size() != 0)
        {
            BKVertex u;

            u = *pivot.begin();

            for(int i = 0; i < u.getDegree(); i++)
                pivot.erase(vertices[u.getNodeConflict(i)]);
        }

        for(set<BKVertex, sortByMDegree>::const_iterator it = pivot.begin(); it != pivot.end(); ++it)
        {
            set<int> C2 = C;
            C2.insert(it->getNumber());
            set<int> P2;

            for(set<int>::const_iterator it2 = P.begin(); it2 != P.end(); ++it2)
            {
                const vector<int> &conflicts = it->getConflicts();

                if(binary_search(conflicts.begin(), conflicts.end(), *it2))
                    P2.insert(*it2);
            }

            set<int> S2;

            for(set<int>::iterator it3 = S.begin(); it3 != S.end(); ++it3)
            {
                const vector<int> &conflicts = it->getConflicts();

                if(binary_search(conflicts.begin(), conflicts.end(), *it3))
                    S2.insert(*it3);
            }

            BronKerbosch(C2, P2, S2, minWeight, timeLimit, init);

            end = clock();
            sec = ((double)(end - init)) / ((double)CLOCKS_PER_SEC);

            if(sec >= timeLimit)
                return 1;

            S.insert(it->getNumber());
            P.erase(it->getNumber());
        }
    }

    return 0;
}

int BKGraph::execute(int minWeight, double timeLimit)
{
    set<int> C, P, S;
    clock_t init;

    for(unsigned int i = 1; i < vertices.size(); i++)
        P.insert(vertices[i].getNumber());

    init = clock();

    int stat = BronKerbosch(C, P, S, minWeight, timeLimit, init);

    return stat;
}

CliqueSet* BKGraph::convertToClqSet()
{
    CliqueSet *clqSet;
    set< set<int> >::iterator extIt;
    set<int>::iterator intIt;

    clqSet = clq_set_create();

    for(extIt = cliques.begin(); extIt != cliques.end(); ++extIt)
    {
        int nodes[extIt->size()];
        int count = 0;

        for(intIt = extIt->begin(); intIt != extIt->end(); ++intIt)
        {
            nodes[count] = *intIt -1 ;
            count++;
        }

        int stat = clq_set_add(clqSet, extIt->size(), nodes, weightCompute(*extIt));

        if(!stat) cout<<"Couldn't insert clique!"<<endl;
    }

    return clqSet;
}

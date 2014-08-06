#ifndef BKGRAPH_HPP_INCLUDED
#define BKGRAPH_HPP_INCLUDED

#include "BKVertex.hpp"
#include <string>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <set>
#include <ctime>

extern "C"
{
    #include "cgraph.h"
    #include "clique.h"
    #include "memory.h"
}

class BKGraph
{
public:
    BKGraph(const std::string&);
    BKGraph(const CGraph*);
    virtual ~BKGraph();
    std::vector<BKVertex> getVertices() const;
    int getNVertices() const;
    int getNEdges() const;
    int getMaxDegree() const;
    std::string getInstance() const;
    BKVertex getVertex(int) const;
    void setNVertices(int);
    void setNEdges(int);
    void setMaxDegree(int);
    int weightCompute(const std::set<int>&);
    int weightEstimate(const std::set<int>&);
    int writeSolutions(const char*);
    int BronKerbosch(std::set<int>, std::set<int>, std::set<int>, int, double, clock_t);
    int execute(int, double);
    CliqueSet* convertToClqSet();

private:
    std::vector<BKVertex> vertices;
    std::set< std::set<int> > cliques;
    int nVertices;
    int nEdges;
    int maxDegree;
    std::string instance;

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

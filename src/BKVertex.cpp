#include "BKVertex.hpp"

using namespace std;

BKVertex::BKVertex(){}

BKVertex::~BKVertex(){}

int BKVertex::getWeight() const { return weight; }

int BKVertex::getDegree() const { return degree; }

int BKVertex::getNumber() const { return number; }

int BKVertex::getMDegree() const { return mdegree; }

int BKVertex::getNumberOfConflicts() const { return conflicts.size(); }

int BKVertex::getNodeConflict(int n) const { return conflicts[n]; }

vector<int> BKVertex::getConflicts() const { return conflicts; }

void BKVertex::setWeight(int w){ weight = w; }

void BKVertex::setDegree(int dg){ degree = dg; }

void BKVertex::setNumber(int num){ number = num; }

void BKVertex::setMDegree(int md){ mdegree = md; }

void BKVertex::setConflicts(const vector<int>& c){ conflicts = c; }

int BKVertex::insertConflict(int c) { conflicts.push_back(c); return 0; }

int BKVertex::removeConflict(int c)
{
    vector<int>::iterator status = find(conflicts.begin(), conflicts.end(), c);
    if(status != conflicts.end())
    {
        conflicts.erase(status);
        return 0;
    }
    else return 1;
}

int BKVertex::clearAllConflicts() { conflicts.clear(); return 0; }

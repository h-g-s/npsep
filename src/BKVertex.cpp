#include "BKVertex.hpp"

using namespace std;

BKVertex::BKVertex(){}

BKVertex::~BKVertex(){}

int BKVertex::getWeight() const { return weight; }

int BKVertex::getDegree() const { return degree; }

int BKVertex::getNumber() const { return number; }

int BKVertex::getMDegree() const { return mdegree; }

void BKVertex::setWeight(int w){ weight = w; }

void BKVertex::setDegree(int dg){ degree = dg; }

void BKVertex::setNumber(int num){ number = num; }

void BKVertex::setMDegree(int md){ mdegree = md; }


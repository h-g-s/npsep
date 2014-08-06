#ifndef VERTEX_H_INCLUDED
#define VERTEX_H_INCLUDED

#include <vector>
#include <algorithm>

class BKVertex
{
private:

    int number;
    int weight;
    int degree;
    int mdegree;
    std::vector<int> conflicts;

    bool operator < ( const BKVertex &other ) const
    {
        return this->number < other.number;
    }

public:

    BKVertex();
    virtual ~BKVertex();

    int getWeight() const;
    int getDegree() const;
    int getNumber() const;
    int getMDegree() const;
    int getNumberOfConflicts() const;
    int getNodeConflict(int) const;
    std::vector<int> getConflicts() const;

    void setWeight(int);
    void setDegree(int);
    void setNumber(int);
    void setMDegree(int);
    void setConflicts(const std::vector<int>&);
    int insertConflict(int);
    int removeConflict(int);
    int clearAllConflicts();

};

#endif // VERTEX_H_INCLUDED

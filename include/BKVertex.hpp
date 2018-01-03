#ifndef VERTEX_H_INCLUDED
#define VERTEX_H_INCLUDED

class BKVertex
{
private:

    int number;
    int weight;
    int degree;
    int mdegree;

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

    void setWeight(int);
    void setDegree(int);
    void setNumber(int);
    void setMDegree(int);

};

#endif // VERTEX_H_INCLUDED

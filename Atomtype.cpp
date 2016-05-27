
#include "Atomtype.h"

Atomtype::Atomtype()
{

}

Atomtype::Atomtype(string name, double c6, double c12)
{
    this->name = name;
    this->c6 = c6;
    this->c12 = c12;
}

string Atomtype::GetName()
{
    return this->name;
}

double Atomtype::GetC6()
{
    return this->c6;
}

double Atomtype::GetC12()
{
    return this->c12;
}

#include "Edge.h"

const int Edge::GetDestID() const
{
    return dest_id;
}

const int Edge::GetStartID() const
{
    return start_id;
}

const double Edge::GetCost() const
{
    return cost;
}

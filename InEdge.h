#ifndef INEDGE_HPP_INCLUDED
#define INEDGE_HPP_INCLUDED
template <typename eD>
class InEdge
{
    public:
        InEdge(int o_id, double _cost, eD _data): origin_id(o_id), cost(_cost), data(_data) {}
        const int GetOriginID() const;
        const double GetCost() const;
        eD & GetData();
    private:
        int origin_id;
        double cost;
        eD data;
};
template <typename eD>
const int InEdge<eD>::GetOriginID() const
{
    return origin_id;
}
template <typename eD>
const double InEdge<eD>::GetCost() const
{
    return cost;
}
template <typename eD>
eD & InEdge<eD>::GetData()
{
    return data;
}

#endif // INEDGE_HPP_INCLUDED

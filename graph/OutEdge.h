#ifndef OUTEDGE_HPP_INCLUDED
#define OUTEDGE_HPP_INCLUDED
template <typename eD>
class OutEdge
{
    public:
        OutEdge(int end_id, double cost, eD _data): dest_id(end_id), data(_data) {}
        const int GetDestID() const;
        const double GetCost() const;
        eD & GetData();
    private:
        int dest_id;
        double cost;
        eD data;
};

template <typename eD>
const int OutEdge<eD>::GetDestID() const
{
    return dest_id;
}

template <typename eD>
eD & OutEdge<eD>::GetData()
{
    return data;
}


#endif // OUTEDGE_HPP_INCLUDED

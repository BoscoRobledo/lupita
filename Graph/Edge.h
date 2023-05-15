#ifndef EDGE_HPP_INCLUDED
#define EDGE_HPP_INCLUDED
template <typename eD>
class Edge
{
    public:
        Edge(int o_id, int d_id, eD _data): origin_id(o_id), destination_id(d_id), data(_data) {}
        const int GetOriginID() const;
        const int GetDestinationID() const;
        eD & GetData();
    private:
        int origin_id;
        int destination_id;
        eD data;
};

template <typename eD>
const int Edge<eD>::GetOriginID() const
{
    return origin_id;
}

template <typename eD>
const int Edge<eD>::GetDestinationID() const
{
    return destination_id;
}

template <typename eD>
eD & Edge<eD>::GetData()
{
    return data;
}

template <typename eD>
bool operator<(Edge<eD>& leftValue,Edge<eD>& rightValue)
{
    return leftValue.GetData()<rightValue.GetData();
}
#endif // EDGE_HPP_INCLUDED

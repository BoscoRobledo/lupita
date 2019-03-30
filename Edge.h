#ifndef EDGE_HPP_INCLUDED
#define EDGE_HPP_INCLUDED
class Edge
{
    public:
        Edge() {}
        Edge(int start_id, int end_id, double cost): dest_id(end_id), start_id(start_id), cost(cost) {}
        const int GetDestID() const;
        const int GetStartID() const;
        const double GetCost() const;

        bool operator < (const Edge& str) const
        {
            return (cost > str.GetCost());
        }
    private:
        int dest_id;
        int start_id;
        double cost;
};
#endif // EDGELIST_HPP_INCLUDED

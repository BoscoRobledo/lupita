#ifndef EDGE_HPP_INCLUDED
#define EDGE_HPP_INCLUDED
class Edge
{
    public:
        Edge() {}
        Edge(int start_id, int end_id, double cost): dest_id(end_id), start_id(start_id) {}
        const int GetDestID() const;
        const int GetStartID() const;
    private:
        int dest_id;
        int start_id;
};
#endif // EDGELIST_HPP_INCLUDED

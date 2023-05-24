#ifndef EDGE_HPP_INCLUDED
#define EDGE_HPP_INCLUDED

template <typename eD>
/** \brief Graph directed edge, containing object of template type eD
 * \param type eD type for object stored in Edge.
 * \author Juan Bosco Robledo Muñoz, 2017
 */
class Edge
{
    public:

        /** \brief Builds Edge and sets its data
         *
         * \param o_id int Origin Vertex id
         * \param d_id int Destination vertex id
         * \param _data eD Data to store in the Edge (typically weight data)
         *
         */
        Edge(int o_id, int d_id, eD _data): origin_id(o_id), destination_id(d_id), data(_data) {}

        /** \brief Gets origin vertex id
         *
         * \return const int Origin Vertex id
         *
         */
        const int GetOriginID() const { return origin_id; }

        /** \brief Gets destination vertex id
         *
         * \return const int destination Vertex id
         *
         */
        const int GetDestinationID() const {return destination_id; }

        /** \brief Gets data object (template type eD) in edge
         *
         * \return eD& reference to data.
         *
         */
        eD & GetData() {return data; }

        /** \brief
         *
         * \param leftValue Edge<eD>&
         * \param rightValue Edge<eD>&
         * \return bool
         *
         */
        friend bool operator<(Edge<eD>& leftValue,Edge<eD>& rightValue)
        {
            return leftValue.GetData()<rightValue.GetData();
        }

    private:

        /** \brief Origin vertex id
         *
         */
        int origin_id;

        /** \brief Destination vertex id
         *
         */
        int destination_id;

        /** \brief Data object (template type eD) in edge.
         *
         */
        eD data;
};

#endif // EDGE_HPP_INCLUDED

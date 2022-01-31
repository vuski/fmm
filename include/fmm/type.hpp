/**
 * Fast map matching.
 *
 * Definition of Data types used in the FMM algorithm
 *
 * @author: Can Yang
 * @version: 2017.11.11
 * 
 * @revision : vuski
 * @version: 2021.11.29
 * 
 */

#ifndef FMM_TYPES_HPP
#define FMM_TYPES_HPP

#include <vector>
#include <list>
#include <string>
#include <set>
#include <unordered_map>


#include "fmm/geometry.hpp"




namespace FMM {






    typedef long long NodeID; /**< Node ID in the network, can be discontinuous
                              int */
    typedef long long EdgeID; /**< Edge ID in the network, can be negative to
                              distinguish edges in two directions */
    typedef unsigned int NodeIndex; /**< Node Index in the network, range
                                     from [0,num_vertices-1 ]*/
    typedef unsigned int EdgeIndex; /**< Edge Index in the network, range
                                     from [0,num_edges-1 ]*/

    /**
     * Vector of node id
     */
    typedef std::vector<NodeID> NodeIDVec;
    /**
     * Map of node index
     */
    typedef std::unordered_map<NodeID,NodeIndex> NodeIndexMap;
    /**
     * Map of edge index
     */
    typedef std::unordered_map<EdgeID,EdgeIndex> EdgeIndexMap;

    /**
     * Road edge class
     */
    struct Edge
    {
      EdgeIndex index; /**< Index of an edge, which is continuous [0,N-1] */
      EdgeID id; /**< Edge ID, can be discontinuous integers */
      NodeIndex source; /**< source node index */
      NodeIndex target; /**< target node index */
      double cost; //¨öA¡Æ¡Ì¨¬n¢¯e, ¡ÆA¢¬¢ç¨¬n¢¯e ¥ìi**< length of the edge polyline */
      FMM::LineString geom; /**< the edge geometry */
    };



  

} // FMM
#endif /* MM_TYPES_HPP */
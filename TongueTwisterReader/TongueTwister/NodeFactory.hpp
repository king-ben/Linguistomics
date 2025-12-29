#ifndef NodeFactory_hpp
#define NodeFactory_hpp

#include <set>
#include <vector>
class Node;



class NodeFactory {

    public:
        static NodeFactory& nodeFactory(void) {
                                static NodeFactory tp;
                                return tp;
                            }
        void                drainPool(void);
        Node*               getNode(void);
        int                 getNumAllocated(void) { return (int) allocated.size(); }
        int                 getNumInPool(void) { return (int)pool.size(); }
        void                returnToPool(Node* n);
    
    private:
                            NodeFactory(void);
                           ~NodeFactory(void);
                            NodeFactory(const NodeFactory& tp) = delete;
        std::vector<Node*>  pool;
        std::set<Node*>     allocated;
};

#endif

#ifndef NodeFactory_hpp
#define NodeFactory_hpp

#include <vector>
class Node;



class NodeFactory {

    public:
        static NodeFactory& nodeFactory(void) 
                                {
                                static NodeFactory tp;
                                return tp;
                                }
        Node*               getNode(void);
        int                 getNumAllocated(void) { return static_cast<int>(allocated.size()); }
        int                 getNumInPool(void) { return static_cast<int>(pool.size()); }
        void                returnToPool(Node* n);
    
    private:
                            NodeFactory(void);
                           ~NodeFactory(void);
                            NodeFactory(const NodeFactory& tp) = delete;
        size_t              initialSize;
        std::vector<Node*>  pool;
        std::vector<Node*>  allocated;
};

#endif

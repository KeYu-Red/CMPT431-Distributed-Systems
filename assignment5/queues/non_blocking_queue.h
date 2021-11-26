#include "../common/allocator.h"
#include <cinttypes>

#define LFENCE asm volatile("lfence" : : : "memory")
#define SFENCE asm volatile("sfence" : : : "memory")
#define LSB 0x0000FFFFFFFFFFFF
#define MSB 0xFFFF000000000000

template <class P>
class pointer_t
{
public:
    P * ptr;

    P* address()
    {
        //Get the address by getting the 48 least significant bits of ptr
        //return ptr&0x0000 ;
        return (P*) ( (uintptr_t)ptr & (uintptr_t)LSB );
    }

    uint count(){
        //Get the count from the 16 most significant bits of ptr
        return (uint)( ((uintptr_t)ptr & (uintptr_t)MSB) >> 48);
    }

    void concate(P* node, uint count)
    {
        ptr = (P*)( ((uintptr_t)count) <<48 | (uintptr_t) node );
    }
};



template <class T>
class Node
{
public:
    T value;
    pointer_t<Node<T>> next;
};

template <class T>
class NonBlockingQueue
{
private:
    pointer_t<Node<T>> q_head;
    pointer_t<Node<T>> q_tail;
    CustomAllocator my_allocator_;
public:
    NonBlockingQueue() : my_allocator_()
    {
        std::cout << "Using NonBlockingQueue\n";
    }

    void initQueue(long t_my_allocator_size){
        std::cout << "Using Allocator\n";
        my_allocator_.initialize(t_my_allocator_size, sizeof(Node<T>));
        // Initialize the queue head or tail here
        Node<T>* newNode = (Node<T>*)my_allocator_.newNode();
        newNode->next.ptr = nullptr;
        q_head.ptr = newNode;
        q_tail.ptr = newNode;
        //my_allocator_.freeNode(newNode);
    }

    void enqueue(T value)
    {
        // Use LFENCE and SFENCE as mentioned in pseudocode
        Node<T>* node = (Node<T>*)my_allocator_.newNode();
        node->value = value;
        node->next.ptr = nullptr;
        pointer_t<Node<T>> tail;
        pointer_t<Node<T>> swap_pointer;
        SFENCE;
        while(true)
        {
            tail = q_tail;
            LFENCE;
            pointer_t<Node<T>> next = tail.address()->next;
            LFENCE;
            if(tail.ptr==q_tail.ptr)
            {
                if(next.address() == nullptr)
                {
                    swap_pointer.concate(node, next.count()+1);
                    if(CAS(&tail.address()->next, next, swap_pointer))
                        break;
                }
                else
                {
                    swap_pointer.concate(next.address(), tail.count()+1);
                    CAS(&q_tail, tail, swap_pointer);
                }
            }
        }
        SFENCE;
        swap_pointer.concate(node, tail.count()+1);
        CAS(&q_tail, tail,swap_pointer);

    }

    bool dequeue(T *value)
    {
        // Use LFENCE and SFENCE as mentioned in pseudocode
        pointer_t<Node<T>> head;
        pointer_t<Node<T>> swap_pointer;
        while(true)
        {
            head = q_head;
            LFENCE;
            pointer_t<Node<T>> tail = q_tail;
            LFENCE;
            pointer_t<Node<T>> next = head.address()->next;
            LFENCE;
            if(head.ptr==q_head.ptr)
            {
                if(head.address() == tail.address())
                {
                    if(next.address() == nullptr)
                        return false;
                    swap_pointer.concate(next.address(), tail.count()+1);
                    CAS(&q_tail, tail, swap_pointer);
                }
                else{
                    *value = next.address()->value;
                    swap_pointer.concate(next.address(), head.count()+1);
                    if(CAS(&q_head, head, swap_pointer))
                        break;

                }
            }
        }
        my_allocator_.freeNode(head.address());
        return true;
    }

    void cleanup()
    {
        my_allocator_.cleanup();
    }

};


#include "../common/allocator.h"
#include <mutex>

template <class T>
class Node
{
public:
    T value;
    Node<T> * next;
};

template <class T>
class OneLockQueue
{
private:
    Node<T>* q_head;
    Node<T>* q_tail;
    std::mutex one_lock;
    CustomAllocator my_allocator_;
    
public:
    OneLockQueue() : my_allocator_()
    {
        std::cout << "Using OneLockQueue\n";
    }

    void initQueue(long t_my_allocator_size){
        std::cout << "Using Allocator\n";
        my_allocator_.initialize(t_my_allocator_size, sizeof(Node<T>));
        // Initialize the queue head or tail here
        Node<T>* newNode = (Node<T>*)my_allocator_.newNode();
        newNode->next = NULL;
        q_head = newNode;
        q_tail = newNode;
    }

    void enqueue(T value)
    {
        one_lock.lock();
        Node<T>* newNode = (Node<T>*)my_allocator_.newNode();
        newNode->next = NULL;
        newNode->value = value;
        q_tail->next = newNode;
        q_tail = newNode;
        one_lock.unlock();
    }

    bool dequeue(T *value)
    {
        one_lock.lock();
        Node<T>* newNode = q_head;
        Node<T>* newHead = q_head->next;
        if(newHead == NULL)
        {
            one_lock.unlock();
            return false;
        }
        q_head = newHead;
        *value = q_head->value;
        my_allocator_.freeNode(newNode);
        one_lock.unlock();
        return true;
    }

    void cleanup()
    {
        my_allocator_.cleanup();
    }
};
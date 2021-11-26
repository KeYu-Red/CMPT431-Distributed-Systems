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
class TwoLockQueue
{
private:
    Node<T>* q_head;
    Node<T>* q_tail;
    std::mutex enqueue_lock, dequeue_lock;
    CustomAllocator my_allocator_;

public:
    TwoLockQueue() : my_allocator_()
    {
        std::cout << "Using TwoLockQueue\n";
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
        enqueue_lock.lock();
        Node<T>* newNode = (Node<T>*)my_allocator_.newNode();
        newNode->next = NULL;
        newNode->value = value;
        q_tail->next = newNode;
        q_tail = newNode;
        enqueue_lock.unlock();
    }

    bool dequeue(T *value)
    {
        dequeue_lock.lock();
        Node<T>* newNode = q_head;
        Node<T>* newHead = q_head->next;
        if(newHead == NULL)
        {
            dequeue_lock.unlock();
            return false;
        }
        q_head = newHead;
        *value = q_head->value;
        my_allocator_.freeNode(newNode);
        dequeue_lock.unlock();
        return true;
    }

    void cleanup()
    {
        my_allocator_.cleanup();
    }
};
#pragma once
#include <atomic>
#include <iostream>
#include <memory>


namespace NParallel {

template <std::size_t N>
class Plane;

template <std::size_t N>
class Optimizer;

template <std::size_t N>
class PointTree;

template <class T>
class AtomicList;


template <class T>
class Node {
protected:
    std::shared_ptr<T> value;
    std::shared_ptr<Node<T>> next;

public:
    template <std::size_t N>
    friend class Plane;

    template <std::size_t N>
    friend class Optimizer;

    template <std::size_t N>
    friend class PointTree;

    template <class M>
    friend class AtomicList;

    Node() { }
    Node(T&& value_) 
        : value(std::make_shared<T>(std::move(value_))) { }

    void GetNodeValue(std::shared_ptr<T>& to) {
        to = value;
    }

    // Not thread safe! You should ensure on your own that only one process is modifying it
    void UpdateValue(T&& value_) {
        value = std::make_shared<T>(std::move(value_));
    }

    void Next(std::shared_ptr<Node<T>>& to) const {
        to = next;
    }

    // Warning: not thread safe: you should not iterate this list
    void InsertAfter(T&& value_) {
        std::shared_ptr<Node<T>> newNode = std::make_shared<Node<T>>(std::move(value_));
        std::shared_ptr<Node<T>> nextNode(next);
        newNode->next = nextNode;
        next = newNode;
    }
};


template <class T>
class AtomicList {
private:
    enum ListStatus {
        FREE,
        WORKING
    };

    std::atomic<std::size_t> size;
    std::atomic<ListStatus> listStatus;
    std::shared_ptr<Node<T>> begin;
    std::shared_ptr<Node<T>> end;

public:
    AtomicList()
        : size(0)
        , listStatus(ListStatus::FREE) { }

    std::size_t Size() const {
        return size;
    }

    // Warning: not thread-safe
    void InsertAfter(std::shared_ptr<Node<T>>& node, T&& value) {
        node->InsertAfter(std::move(value));
        ++size;
    }

    void InsertBack(T&& value) {
        std::shared_ptr<Node<T>> newNode{std::make_shared<Node<T>>(std::move(value))};
        ListStatus curStatus{ListStatus::FREE};
        if (!listStatus.compare_exchange_strong(curStatus, ListStatus::WORKING)) {
            curStatus = ListStatus::FREE;
            while (!listStatus.compare_exchange_weak(curStatus, ListStatus::WORKING)) {
                curStatus = ListStatus::FREE;
            }
        }
        if (size == 0) {
            begin = end = newNode;
        } else {
            end->next = newNode;
            end = newNode;
        }
        ++size;
        curStatus = ListStatus::WORKING;
        assert(listStatus.compare_exchange_strong(curStatus, ListStatus::FREE));
    }

    void GetBegin(std::shared_ptr<Node<T>>& to) {
        to = begin;
    }

    ~AtomicList() {
        if (size > 2) {
            --size;
            while (--size) {
                begin->next = std::move(begin->next->next);
            }
        }
    }
};

} // end namespace NParallel

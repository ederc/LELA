/*
 * Allocator.h
 * Copyright 2012 Martani Fayssal (UPMC University Paris 06 / INRIA)
 *
 * Created on: 1 august 2012
 *      Author: martani (UPMC University Paris 06 / INRIA)
 * 
 */

#ifndef ALLOCATOR_H_
#define ALLOCATOR_H_


#include <stdlib.h>
#include <stddef.h>

template <class T , int Alignment=16>
class aligned_allocator
{
public:

    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;
    typedef T*        pointer;
    typedef const T*  const_pointer;
    typedef T&        reference;
    typedef const T&  const_reference;
    typedef T         value_type;


    template <class U>
    struct rebind {
        typedef aligned_allocator<U> other;
    };


    pointer address ( reference value ) const {
        return &value;
    };

    const_pointer address ( const_reference value ) const {
        return &value;
    };


    aligned_allocator() throw() {
    };

    aligned_allocator ( const aligned_allocator& ) throw() {
    };

    template <class U>
    aligned_allocator ( const aligned_allocator<U>& ) throw() {
    };

    ~aligned_allocator() throw() {
    };

    //max capacity
    size_type max_size () const throw() {
        return (size_type)-1;
    };


    pointer allocate ( size_type size, const_pointer *hint = 0 ) {
	    pointer p;
	    posix_memalign((void**)&p, Alignment, size * sizeof (T));
	    
	    return p;
        //return ( pointer ) _aligned_malloc ( num*sizeof ( T ),Alignment );
        //return ( pointer ) memalign(Alignment, num*sizeof ( T ));
    };


    void construct ( pointer p, const T& value ) {

        // memcpy( p, &value, sizeof T );
        *p=value;
        //  new ( (void *) p ) T ( value );
    };


    void destroy ( pointer p ) {

        p->~T();
    };


    void deallocate ( pointer p, size_type num ) {

        //_aligned_free ( p );
        free(p);
    };
};


template <class T, class U, int Alignment>
bool operator==(const aligned_allocator<T, Alignment>& a,
           const aligned_allocator<U, Alignment>& b)
{
	return true;
}

template <class T, class U, int Alignment>
bool operator!=(const aligned_allocator<T, Alignment>& a,
           const aligned_allocator<U, Alignment>& b)
{
	return !(a == b);
}


#endif 	//ALLOCATOR_H

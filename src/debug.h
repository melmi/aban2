/**
 * debug.h
 *
 * by Mohammad Elmi
 * Created on 16 Sep 2014
 */

#ifndef _DEBUG_H_
#define _DEBUG_H_

#include <iostream>
#include <typeinfo>

#define watch(x) print(#x,x)

template<class T>
void print(std::string s, T o)
{
    std::cout << s << ": " << o << std::endl << std::flush;
}

#endif

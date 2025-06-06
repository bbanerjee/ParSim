// Copyright 2015, Tobias Hermann and the FunctionalPlus contributors.
// https://github.com/Dobiasd/FunctionalPlus
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include <submodules/fplus/internal/function_traits_asserts.hpp>

namespace fplus
{
namespace internal
{

struct bind_1st_of_2_tag
{
};

struct bind_2nd_of_2_tag
{
};

struct bind_1st_of_3_tag
{
};

struct bind_1st_and_2nd_of_3_tag
{
};

struct bind_2nd_and_3rd_of_3_tag
{
};

template <typename F, typename X, typename Y>
struct function_traits_asserts<bind_1st_of_2_tag, F, X, Y>
{
    static_assert(utils::function_traits<F>::arity == 2,
                  "Function must take two parameters.");
    typedef typename utils::function_traits<F>::template arg<0>::type FIn0;
    typedef typename utils::function_traits<F>::template arg<1>::type FIn1;
    static_assert(std::is_convertible<X, FIn0>::value,
                  "Function can not take bound parameter type");
    static_assert(std::is_convertible<Y, FIn1>::value,
                  "Function can not take provided parameter type");
};

template <typename F, typename X, typename Y>
struct function_traits_asserts<bind_2nd_of_2_tag, F, X, Y>
{
    static_assert(utils::function_traits<F>::arity == 2,
                  "Function must take two parameters.");
    typedef typename utils::function_traits<F>::template arg<0>::type FIn0;
    typedef typename utils::function_traits<F>::template arg<1>::type FIn1;
    static_assert(std::is_convertible<X, FIn0>::value,
                  "Function can not take provided parameter type");
    static_assert(std::is_convertible<Y, FIn1>::value,
                  "Function can not take bound parameter type");
};

template <typename F, typename X, typename Y, typename Z>
struct function_traits_asserts<bind_1st_of_3_tag, F, X, Y, Z>
{
    static_assert(utils::function_traits<F>::arity == 3,
                  "Function must take three parameters.");
    typedef typename utils::function_traits<F>::template arg<0>::type FIn0;
    typedef typename utils::function_traits<F>::template arg<1>::type FIn1;
    typedef typename utils::function_traits<F>::template arg<2>::type FIn2;
    static_assert(std::is_convertible<X, FIn0>::value,
                  "Function can not take bound parameter type");
    static_assert(std::is_convertible<Y, FIn1>::value,
                  "Function can not take provided first parameter type");
    static_assert(std::is_convertible<Z, FIn2>::value,
                  "Function can not take provided second parameter type");
};

template <typename F, typename X, typename Y, typename Z>
struct function_traits_asserts<bind_1st_and_2nd_of_3_tag, F, X, Y, Z>
{
    static_assert(utils::function_traits<F>::arity == 3,
                  "Function must take three parameters.");
    typedef typename utils::function_traits<F>::template arg<0>::type FIn0;
    typedef typename utils::function_traits<F>::template arg<1>::type FIn1;
    typedef typename utils::function_traits<F>::template arg<2>::type FIn2;
    static_assert(std::is_convertible<X, FIn0>::value,
                  "Function can not take first bound parameter type");
    static_assert(std::is_convertible<Y, FIn1>::value,
                  "Function can not take second bound parameter type");
    static_assert(std::is_convertible<Z, FIn2>::value,
                  "Function can not take provided parameter type");
};

template <typename F, typename X, typename Y, typename Z>
struct function_traits_asserts<bind_2nd_and_3rd_of_3_tag, F, X, Y, Z>
{
    static_assert(utils::function_traits<F>::arity == 3,
                  "Function must take three parameters.");
    typedef typename utils::function_traits<F>::template arg<0>::type FIn0;
    typedef typename utils::function_traits<F>::template arg<1>::type FIn1;
    typedef typename utils::function_traits<F>::template arg<2>::type FIn2;
    static_assert(std::is_convertible<X, FIn0>::value,
                  "Function can not take provided parameter type");
    static_assert(std::is_convertible<Y, FIn1>::value,
                  "Function can not take second bound parameter type");
    static_assert(std::is_convertible<Z, FIn2>::value,
                  "Function can not take first bound parameter type");
};

}
}

#ifndef UTILS_TUPLE_HPP
#define UTILS_TUPLE_HPP

namespace oopse {

template <class T1, class T2, class T3>
struct tuple3 {
    typedef T1 first_type;
    typedef T2 second_type;
    typedef T3 third_type;

    T1 first;
    T2 second;
    T3 third;
    tuple3() {}
    tuple3(const T1& a, const T2& b, const T3& c) : first(a), second(b), third(c) {}
};

template <class T1, class T2, class T3>
tuple3<T1,T2,T3> make_tuple3( T1 t1, T2 t2, T3 t3 ) {
    return tuple3<T1,T2,T3>( t1, t2, t3 ); 
}


template <class T1, class T2, class T3, class T4>
struct tuple4 {
    typedef T1 first_type;
    typedef T2 second_type;
    typedef T3 third_type;
    typedef T4 fourth_type;

    T1 first;
    T2 second;
    T3 third;
    T4 fourth;
    tuple4() {}
    tuple4(const T1& a, const T2& b, const T3& c, const T4& d)
    : first(a), second(b), third(c), fourth(d) {}
};

template <class T1, class T2, class T3, class T4>
tuple4<T1,T2,T3,T4> make_tuple4( T1 t1, T2 t2, T3 t3, T4 t4 ) {
    return tuple4<T1,T2,T3,T4>( t1, t2, t3, t4 );
}

}
#endif //UTILS_TUPLE_HPP


//  Boost string_algo library classification.hpp header file  ---------------------------//

//  Copyright Pavol Droba 2002-2003. Use, modification and
//  distribution is subject to the Boost Software License, Version
//  1.0. (See accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt)

//  See http://www.boost.org for updates, documentation, and revision history.
 #ifndef UTILS_PREDICATE_HPP
 #define UTILS_PREDICATE_HPP
#include <locale>



#if defined(__SUNPRO_CC)
#  define USE_FACET(Type, loc) std::use_facet(loc, static_cast<Type*>(0))
#else
#  define USE_FACET(Type, loc) std::use_facet< Type >(loc)
#endif

namespace OpenMD {

template<typename Derived>
struct PredFacade{};

struct CharClassification: public PredFacade<CharClassification> {

    CharClassification(std::ctype_base::mask type, std::locale const & loc = std::locale()) :
        type_(type), loc_(loc) {}

    template<typename CharT>
    bool operator()( CharT c ) const {
        return USE_FACET(std::ctype<CharT>, loc_).is( type_, c );
    }

    private:
    const std::ctype_base::mask type_;
    const std::locale loc_;
};

template<typename CharT>
struct FromRangeFunctor : public PredFacade< FromRangeFunctor<CharT> > {

    FromRangeFunctor( CharT from, CharT to ) : from_(from), to_(to) {}
    
    template<typename Char2T>
    bool operator()( Char2T c ) const {
        return ( from_ <= c ) && ( c <= to_ ); 
    }

private:
    CharT from_;
    CharT to_;
};

template<typename Pred1T, typename Pred2T>
struct PredAndFunctor : public PredFacade< PredAndFunctor<Pred1T,Pred2T> > {
public:

    PredAndFunctor( Pred1T pred1, Pred2T pred2 ) : pred1_(pred1), pred2_(pred2) {}

    template<typename CharT>
    bool operator()( CharT c ) const {
        return pred1_(c) && pred2_(c);
    }

private:
    Pred1T pred1_;
    Pred2T pred2_;
};

template<typename Pred1T, typename Pred2T>
struct PredOrFunctor : public PredFacade< PredOrFunctor<Pred1T,Pred2T> > {
public:

    PredOrFunctor( Pred1T pred1, Pred2T pred2 ) : pred1_(pred1), pred2_(pred2) {}

    template<typename CharT>
    bool operator()( CharT c ) const {
        return pred1_(c) || pred2_(c);
    }

private:
    Pred1T pred1_;
    Pred2T pred2_;
};

template< typename PredT >
struct PredNotFunctor : public PredFacade< PredNotFunctor<PredT> > {
public:

    PredNotFunctor( PredT pred ) : pred_(pred) {}

    template<typename CharT>
    bool operator()( CharT c ) const {
        return !pred_(c);
    }

private:
    PredT pred_;
};
       
  inline CharClassification isSpace(const std::locale& loc=std::locale()){
    return CharClassification(std::ctype_base::space, loc);
 }

 inline CharClassification isAlnum(const std::locale& loc=std::locale()){
    return CharClassification(std::ctype_base::alnum, loc);
 }

 inline CharClassification isAlpha(const std::locale& loc=std::locale()){
    return CharClassification(std::ctype_base::alpha, loc);
 }

inline CharClassification isCntrl(const std::locale& loc=std::locale()) {
    return CharClassification(std::ctype_base::cntrl, loc);
}

inline CharClassification isDigit(const std::locale& loc=std::locale()) {
    return CharClassification(std::ctype_base::digit, loc);
}

inline CharClassification isGraph(const std::locale& loc=std::locale()) {
    return CharClassification(std::ctype_base::graph, loc);
}

inline CharClassification isLower(const std::locale& loc=std::locale()) {
    return CharClassification(std::ctype_base::lower, loc);
}

inline CharClassification isPrint(const std::locale& loc=std::locale()) {
    return CharClassification(std::ctype_base::print, loc);
}

inline CharClassification isPunct(const std::locale& loc=std::locale()) {
    return CharClassification(std::ctype_base::punct, loc);
}      

inline CharClassification isUpper(const std::locale& loc=std::locale()) {
    return CharClassification(std::ctype_base::upper, loc);
}    

inline CharClassification isXDigit(const std::locale& loc=std::locale()) {
    return CharClassification(std::ctype_base::xdigit, loc);
} 


template<typename CharT>
inline FromRangeFunctor<CharT> isFromRange(CharT from, CharT to) {
    return FromRangeFunctor<CharT>(from, to); 
}

template<typename Pred1T, typename Pred2T>
inline PredAndFunctor<Pred1T, Pred2T>
operator&&(const PredFacade<Pred1T>& Pred1, const PredFacade<Pred2T>& Pred2) {    
    // Doing the static_cast with the pointer instead of the reference
    // is a workaround for some compilers which have problems with
    // static_cast's of template references, i.e. CW8. /grafik/
    return PredAndFunctor<Pred1T,Pred2T>(
        *static_cast<const Pred1T*>(&Pred1), 
        *static_cast<const Pred2T*>(&Pred2) );
}


template<typename Pred1T, typename Pred2T>
inline PredOrFunctor<Pred1T, Pred2T> 
operator||( const PredFacade<Pred1T>& Pred1, const PredFacade<Pred2T>& Pred2 ) {    
    // Doing the static_cast with the pointer instead of the reference
    // is a workaround for some compilers which have problems with
    // static_cast's of template references, i.e. CW8. /grafik/
    return PredOrFunctor<Pred1T,Pred2T>(
        *static_cast<const Pred1T*>(&Pred1), 
        *static_cast<const Pred2T*>(&Pred2));
}

template<typename PredT>
inline PredNotFunctor<PredT>
operator!( const PredFacade<PredT>& Pred ) {
    // Doing the static_cast with the pointer instead of the reference
    // is a workaround for some compilers which have problems with
    // static_cast's of template references, i.e. CW8. /grafik/
    return PredNotFunctor<PredT>(*static_cast<const PredT*>(&Pred)); 
}

}

 #endif 

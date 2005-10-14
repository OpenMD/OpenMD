/*
 * Copyright (c) 2005 The University of Notre Dame. All Rights Reserved.
 *
 * The University of Notre Dame grants you ("Licensee") a
 * non-exclusive, royalty free, license to use, modify and
 * redistribute this software in source and binary code form, provided
 * that the following conditions are met:
 *
 * 1. Acknowledgement of the program authors must be made in any
 *    publication of scientific results based in part on use of the
 *    program.  An acceptable form of acknowledgement is citation of
 *    the article in which the program was described (Matthew
 *    A. Meineke, Charles F. Vardeman II, Teng Lin, Christopher
 *    J. Fennell and J. Daniel Gezelter, "OOPSE: An Object-Oriented
 *    Parallel Simulation Engine for Molecular Dynamics,"
 *    J. Comput. Chem. 26, pp. 252-271 (2005))
 *
 * 2. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 3. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the
 *    distribution.
 *
 * This software is provided "AS IS," without a warranty of any
 * kind. All express or implied conditions, representations and
 * warranties, including any implied warranty of merchantability,
 * fitness for a particular purpose or non-infringement, are hereby
 * excluded.  The University of Notre Dame and its licensors shall not
 * be liable for any damages suffered by licensee as a result of
 * using, modifying or distributing the software or its
 * derivatives. In no event will the University of Notre Dame or its
 * licensors be liable for any lost revenue, profit or data, or for
 * direct, indirect, special, consequential, incidental or punitive
 * damages, however caused and regardless of the theory of liability,
 * arising out of the use of or inability to use software, even if the
 * University of Notre Dame has been advised of the possibility of
 * such damages.
 */
 
#ifndef IO_PARAMCONSTRAINT_HPP
#define IO_PARAMCONSTRAINT_HPP

/**
  * This class allows to recognize constraint predicates, so that they can be combined using
  * composition operators. Every constraint predicate must be derived from this class
  */
template<typename Derived>
struct ParamConstraintFacade {
};


struct NonConstraint : public ParamConstraintFacade<NonConstraint>{
    template<typename DataType>
    bool operator()( DataType data ) const {
        return true;
    }
};

struct NotEmptyConstraint : public ParamConstraintFacade<NotEmptyConstraint>{
    bool operator()( const std::string data ) const {
        return !data.empty();
    }
};

struct ZeroConstraint : public ParamConstraintFacade<ZeroConstraint>{
    template<typename DataType>
    bool operator()( DataType data ) const
    {
        return data == 0;
    }
};

struct NotZeroConstraint : public ParamConstraintFacade<NotZeroConstraint>{
    template<typename DataType>
    bool operator()( DataType data ) const
    {
        return data != 0;
    }
};

struct PositiveConstraint : public ParamConstraintFacade<PositiveConstraint>{
    template<typename DataType>
    bool operator()( DataType data ) const
    {
        return data > 0;
    }
};

struct NotPositiveConstraint : public ParamConstraintFacade<NotPositiveConstraint>{
    template<typename DataType>
    bool operator()( DataType data ) const
    {
        return data <= 0;
    }
};
struct NotNegativeConstraint : public ParamConstraintFacade<NotNegativeConstraint>{
    template<typename DataType>
    bool operator()( DataType data ) const
    {
        return data >= 0;
    }
};

struct NegativeConstraint : public ParamConstraintFacade<NegativeConstraint>{
    template<typename DataType>
    bool operator()( DataType data ) const
    {
        return data < 0;
    }
};

struct NotNegativeConstraint : public ParamConstraintFacade<NotNegativeConstraint>{
    template<typename DataType>
    bool operator()( DataType data ) const
    {
        return true;
    }
};

template<typename T>
struct LessThanConstraint : public ParamConstraintFacade<LessThanConstraint> {
    
    LessThanConstraint(T rhs) : rhs_(rhs){}
    template<typename DataType>
    bool operator()( DataType data ) const {
        return data < rhs_; 
    }
    private:
        T rhs_;        
};

template<typename T>
struct LessEqualConstraint : public ParamConstraintFacade<LessEqualConstraint> {
    
    LessEqualConstraint(T rhs) : rhs_(rhs){}
    template<typename DataType>
    bool operator()( DataType data ) const {
        return data <= rhs_; 
    }
    private:
        T rhs_;        
};

template<typename T>
struct EqualConstraint : public ParamConstraintFacade<EqualConstraint> {
    
    EqualConstraint(T rhs) : rhs_(rhs){}
    template<typename DataType>
    bool operator()( DataType data ) const {
        return data == rhs_; 
    }
    private:
        T rhs_;        
};

// class_and composition predicate
template<typename Cons1T, typename Cons2T>
struct AndParamConstraint:
    public ParamConstraintFacade< AndParamConstraint<Cons1T,Cons2T> > {
public:

    AndParamConstraint( Cons1T cons1, Cons2T cons2 ) :
        cons1_(cons1, cons2_(cons2) {}

    template<typename DataType>
    bool operator()( DataType data ) const {
        return cons1_(data) && cons2_(data);
    }

private:
    Cons1T cons1_;
    Cons2T cons2_;
};



template<typename Cons1T, typename Cons2T>
struct OrParamConstraint:
    public ParamConstraintFacade< OrParamConstraint<Cons1T,Cons2T> > {
public:


    OrParamConstraint( Cons1T cons1, Cons2T cons2 ) :
        cons1_(cons1, cons2_(cons2) {}

    template<typename DataType>
    bool operator()( DataType data ) const {
        return cons1_(data) || cons2_(data);
    }

private:
    Cons1T cons1_;
    Cons2T cons2_;
};

template<typename ConsT>
struct NotParamConstraint:
    public ParamConstraintFacade< NotParamConstraint<ConsT> >
{
public:


    NotParamConstraint( ConsT cons) :
        cons_(cons) {}

    template<typename DataType>
    bool operator()( DataType data ) const {
        return !cons_(data);
    }

private:
    ConsT cons_;
};    


template<typename Cons1T, typename Cons2T>
inline AndParamConstraint<Cons1T, Cons2T>
operator &&(const ParamConstraintFacade<Cons1T>& cons1,  const ParamConstraintFacade<Cons2T>& cons2 ) {    

    return AndParamConstraint<Cons1T,Cons2T>(
        *static_cast<const Cons1T*>(&cons1), 
        *static_cast<const Cons2T*>(&cons2) );
}

template<typename Cons1T, typename Cons2T>
inline OrParamConstraint<Cons1T, Cons2T>
operator ||( const ParamConstraintFacade<Cons1T>& cons1, const ParamConstraintFacade<Cons2T>& cons2 ) {    

    return OrParamConstraint<Cons1T,Cons2T>(
        *static_cast<const Cons1T*>(&cons1), 
        *static_cast<const Cons2T*>(&cons2) );
}


template<typename ConsT>
inline NotParamConstraint<ConsT>
operator !( const ParamConstraintFacade<ConsT>& cons ) {

    return NotParamConstraint<ConsT>(*static_cast<const ConsT*>(&cons)); 
}



#endif

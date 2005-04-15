////////////////////////////////////////////////////////////////////////////////
// The Loki Library
// Copyright (c) 2001 by Andrei Alexandrescu
// This code accompanies the book:
// Alexandrescu, Andrei. "Modern C++ Design: Generic Programming and Design 
//     Patterns Applied". Copyright (c) 2001. Addison-Wesley.
// Permission to use, copy, modify, distribute and sell this software for any 
//     purpose is hereby granted without fee, provided that the above copyright 
//     notice appear in all copies and that both that copyright notice and this 
//     permission notice appear in supporting documentation.
// The author or Addison-Wesley Longman make no representations about the 
//     suitability of this software for any purpose. It is provided "as is" 
//     without express or implied warranty.
////////////////////////////////////////////////////////////////////////////////
#ifndef _TYPE_INFO_H_
#define _TYPE_INFO_H_
#include <typeinfo>


class TypeInfo
{
 public:
  // Constructors
  TypeInfo(); // needed for containers
  TypeInfo(const std::type_info&); // non-explicit

  // Access for the wrapped std::type_info
  const std::type_info& Get() const;
  // Compatibility functions
  bool before(const TypeInfo& rhs) const;
  const char* name() const;

 private:
  const std::type_info* pInfo_;
};

// Implementation

inline TypeInfo::TypeInfo(){
  class Nil {};
  pInfo_ = &typeid(Nil);
}

inline TypeInfo::TypeInfo(const std::type_info& ti): pInfo_(&ti)
{ }

inline bool TypeInfo::before(const TypeInfo& rhs) const{
  return pInfo_->before(*rhs.pInfo_);
}

inline const std::type_info& TypeInfo::Get() const{
  return *pInfo_;
}

inline const char* TypeInfo::name() const{
  return pInfo_->name();
}

// Comparison operators

inline bool operator==(const TypeInfo& lhs, const TypeInfo& rhs){
  return lhs.Get() == rhs.Get(); 
}

inline bool operator<(const TypeInfo& lhs, const TypeInfo& rhs){
  return lhs.before(rhs); 
}

inline bool operator!=(const TypeInfo& lhs, const TypeInfo& rhs){
  return !(lhs == rhs); 
}    

inline bool operator>(const TypeInfo& lhs, const TypeInfo& rhs){ 
  return rhs < lhs;
}

inline bool operator<=(const TypeInfo& lhs, const TypeInfo& rhs){
  return !(lhs > rhs);
}
 
inline bool operator>=(const TypeInfo& lhs, const TypeInfo& rhs){
  return !(lhs < rhs); 
}



#endif //iindef _TYPE_INFO_H_

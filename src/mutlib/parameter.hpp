/*
 * Copyright (c) Medical Research Council 2001. All rights reserved.
 *
 * Permission to use, copy, modify and distribute this software and its
 * documentation for any purpose is hereby granted without fee, provided that
 * this copyright and notice appears in all copies.
 *
 * This file was written as part of the Staden Package at the MRC Laboratory
 * of Molecular Biology, Hills Road, Cambridge, CB2 2QH, United Kingdom.
 *
 * MRC disclaims all warranties with regard to this software.
 *
 */


#ifndef _NUMERIC_PARAMETER_HPP_
#define _NUMERIC_PARAMETER_HPP_


#include <cassert>
#include <cstring>      // For strlen(), strcpy()



template <typename T>
class NumericParameter
{
 public:
    // Constructors/Destructor
    NumericParameter();
    NumericParameter( T nDefaultValue, T nMinValue, T nMaxValue, const char* pName );
   ~NumericParameter()                       { delete [] m_pName; m_pName=0; }
 
    

 public:
    // Services
    bool        IsValid() const;
    const char* Name() const                 { return m_pName; }
    T           Value() const                { return m_nValue; }
    void        Value( T v )                 { m_nValue = v; }
    T           Minimum() const              { return m_nMinimumValue; }
    T           Maximum() const              { return m_nMaximumValue; }
    T           Default() const              { return m_nDefaultValue; }
    void        RevertToDefault() const      { m_nValue=m_nDefaultValue; }
    void        RevertToMaximum() const      { m_nValue=m_nMaximumValue; }
    void        RevertToMinimum() const      { m_nValue=m_nMinimumValue; }



 private:
    // Data
    char* m_pName;
    T     m_nValue;
    T     m_nDefaultValue;
    T     m_nMinimumValue;
    T     m_nMaximumValue;
};



//--------------
// Constructors
//--------------

template <typename T>
NumericParameter<T>::NumericParameter()
{
    m_nValue        = -1.0;
    m_nDefaultValue = -1.0;
    m_nMinimumValue = -1.0;
    m_nMaximumValue = -1.0;
    m_pName         = 0;
}

template <typename T>
NumericParameter<T>::NumericParameter( T nDefaultValue, T nMinValue, T nMaxValue, const char* pName )
{
    assert(pName != NULL);
    m_nValue        = nDefaultValue;
    m_nDefaultValue = nDefaultValue;
    m_nMinimumValue = nMinValue;
    m_nMaximumValue = nMaxValue;
    if( !pName )
        pName = 0;
    else
    {
        m_pName = new char[ std::strlen(pName)+1 ];
        std::strcpy( m_pName, pName );
    }
}



//---------
// IsValid
//---------

template <typename T>
bool NumericParameter<T>::IsValid() const
{
    if( (m_nValue>m_nMaximumValue) || (m_nValue<m_nMinimumValue) )
        return false;
    return true;
}



#endif

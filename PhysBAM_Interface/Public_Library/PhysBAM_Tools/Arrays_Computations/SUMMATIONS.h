//#####################################################################
// Copyright 2009, Michael Lentine.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SUMMATIONS
//#####################################################################
#ifndef __SUMMATIONS__
#define __SUMMATIONS__

#include <PhysBAM_Tools/Utilities/TYPE_UTILITIES.h>
namespace PhysBAM{
template<class T,class T_ARRAY,class ID>  class ARRAY_BASE;

namespace ARRAYS_COMPUTATIONS
{
    template<class T,class T_ARRAY,class ID>
    T Sum(const ARRAY_BASE<T,T_ARRAY,ID>& a)
    {const T_ARRAY& self=a.Derived();T result=T();ID m=self.Size();for(ID i(1);i<=m;i++) result+=self(i);return result;}

    template<class T,class T_ARRAY,class ID>
    double Sum_Double_Precision(const ARRAY_BASE<T,T_ARRAY,ID>& a)
    {const T_ARRAY& self=a.Derived();double result=0;ID m=self.Size();for(ID i(1);i<=m;i++) result+=self(i);return result;}

    template<class T,class T_ARRAY,class ID>
    T Sumabs(const ARRAY_BASE<T,T_ARRAY,ID>& a)
    {const T_ARRAY& self=a.Derived();T result=T();ID m=self.Size();for(ID i(1);i<=m;i++) result+=abs(self(i));return result;}
    
    template<class T,class T_ARRAY,class ID>
    T Average(const ARRAY_BASE<T,T_ARRAY,ID>& a)
    {const T_ARRAY& self=a.Derived();return self.Size()?Sum(a)/typename ARRAY_BASE<T,T_ARRAY,ID>::SCALAR(self.Size()):T();}
    
    template<class T,class T_ARRAY1,class ID,class T2,class T_ARRAY2>
    T2 Weighted_Sum(const ARRAY_BASE<T,T_ARRAY1,ID>& weights,const ARRAY_BASE<T2,T_ARRAY2,ID>& array)
    {STATIC_ASSERT_SAME(T,typename T_ARRAY2::SCALAR);assert(weights.Size()==array.Size());
    T2 result((T2()));ID m=array.Size();for(ID i(1);i<=m;i++) result+=weights(i)*array(i);return result;}
}
}
#endif

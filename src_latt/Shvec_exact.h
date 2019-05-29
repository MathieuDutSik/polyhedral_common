#ifndef INCLUDE_TEMPLATE_SHVEC_H
#define INCLUDE_TEMPLATE_SHVEC_H


#include "LatticeDefinitions.h"
#include "NumberTheory.h"
#include "MAT_Matrix.h"

#define CHECK_BASIC_CONSISTENCY



namespace TempShvec_globals {
  const int TEMP_SHVEC_MODE_UNDEF=-1;
  const int TEMP_SHVEC_MODE_BOUND=0;
  const int TEMP_SHVEC_MODE_SHORTEST_VECTORS=1;
  const int TEMP_SHVEC_MODE_MINIMUM=2;
  const int TEMP_SHVEC_MODE_THETA_SERIES=3;
  const int TEMP_SHVEC_MODE_VINBERG=4;
  const int STOP_COMPUTATION=666;
  const int NORMAL_TERMINATION_COMPUTATION=555;
}

template<typename T>
struct T_shvec_request {
  int dim;
  int mode;
  T bound;
  MyVector<T> coset;
  MyMatrix<T> gram_matrix;
};

template<typename T, typename Tint>
struct T_shvec_info {
  T_shvec_request<T> request;
  std::vector<MyVector<Tint>> short_vectors;
  T minimum;
};



// We return
// floor(sqrt(A) + epsilon + B) 
template<typename T>
int Infinitesimal_Floor_V1(T const& a, T const& b)
{
  double epsilon=0.000000001;
#ifdef CHECK_BASIC_CONSISTENCY
  if (a < 0) {
    std::cerr << "Error in Infinitesimal_Floor_V1\n";
    std::cerr << "calling with a<0 which gives NAN with sqrt\n";
    std::cerr << "a=" << a << "\n";
    throw TerminalException{1};
  }
#endif
  double a_doubl = UniversalTypeConversion<double,T>(a);
  double b_doubl = UniversalTypeConversion<double,T>(b);
  //  std::cerr << "a_doubl=" << a_doubl << "\n";
  //  std::cerr << "b_doubl=" << b_doubl << "\n";
  double alpha=sqrt(a_doubl) + epsilon + b_doubl;
  //  std::cerr << "alpha=" << alpha << "\n";
  double eD1=floor(alpha);
  //  std::cerr << "eD1=" << eD1 << "\n";
  long int eD2=lround(eD1);
  //  std::cerr << "eD2=" << eD2 << "\n";
  int eD3=eD2;
  //  std::cerr << "eD3=" << eD3 << "\n";
  return eD3;
}

template<typename T>
int Infinitesimal_Ceil_V1(T const& a, T const& b)
{
  double epsilon=0.000000001;
#ifdef CHECK_BASIC_CONSISTENCY
  if (a < 0) {
    std::cerr << "Error in Infinitesimal_Ceil_V1\n";
    std::cerr << "calling with a<0 which gives NAN with sqrt\n";
    std::cerr << "a=" << a << "\n";
    throw TerminalException{1};
  }
#endif
  double a_doubl = UniversalTypeConversion<double,T>(a);
  double b_doubl = UniversalTypeConversion<double,T>(b);
  //  std::cerr << "a_doubl=" << a_doubl << "\n";
  //  std::cerr << "b_doubl=" << b_doubl << "\n";
  double alpha=-sqrt(a_doubl) - epsilon + b_doubl;
  //  std::cerr << "alpha=" << alpha << "\n";
  double eD1=ceil(alpha);
  //  std::cerr << "eD1=" << eD1 << "\n";
  long int eD2=lround(eD1);
  //  std::cerr << "eD2=" << eD2 << "\n";
  int eD3=eD2 - 1;
  //  std::cerr << "eD3=" << eD3 << "\n";
  return eD3;
}




// We return floor(sqrt(a) + b)
// n=floor(sqrt(a) + b) is equivalent to
// n<= sqrt(a) + b < n+1
// And so to n - b <= sqrt(a) (and opposite for n+1)
// And so to (n-b)^2 <= a
template<typename T, typename Tint>
Tint Infinitesimal_Floor(T const& a, T const& b)
{
  double epsilon=0.000000001;
#ifdef CHECK_BASIC_CONSISTENCY
  if (a < 0) {
    std::cerr << "Error in Infinitesimal_Floor\n";
    std::cerr << "calling with a<0 which gives NAN with sqrt\n";
    std::cerr << "a=" << a << "\n";
    throw TerminalException{1};
  }
#endif
  //  std::cerr << "a=" << a << " b=" << b << "\n";
  double a_doubl = UniversalTypeConversion<double,T>(a);
  double b_doubl = UniversalTypeConversion<double,T>(b);
  double alpha=sqrt(a_doubl) + epsilon + b_doubl;
  //  std::cerr << "alpha=" << alpha << "\n";
  double eD1=floor(alpha);
  //  std::cerr << "eD1=" << eD1 << "\n";
  long int eD2=lround(eD1);
  //  std::cerr << "eD2=" << eD2 << "\n";
  Tint eReturn=eD2;
  //  std::cerr << "initial value eReturn=" << eReturn << "\n";
  auto f=[&](Tint const& x) -> bool {
    T eDiff=x - b;
    if (eDiff <= 0)
      return true;
    if (eDiff*eDiff <= a)
      return true;
    return false;
  };
  //  std::cerr << "Infinitesimal_floor, before while loop\n";
  while (true) {
    //    std::cerr << "eReturn=" << eReturn << "\n";
    bool test1=f(eReturn);
    bool test2=f(eReturn+1);
    //    std::cerr << "test1=" << test1 << " test2=" << test2 << "\n";
    if (test1 && !test2)
      break;
    if (!test1)
      eReturn--;
    if (test2)
      eReturn++;
  }
  //  std::cerr << "Infinitesimal_floor, after while loop\n";
  return eReturn;
}

// We return floor(sqrt(a) + b)
// n=ceil(-sqrt(a) + b) is equivalent to
// n-1 < -sqrt(a) + b <= n
// And so to -sqrt(a) <= n - b  (and opposite for n-1)
// And so to (n-b)^2 <= a (and opposite for n-1)
template<typename T, typename Tint>
Tint Infinitesimal_Ceil(T const& a, T const& b)
{
  double epsilon=0.000000001;
#ifdef CHECK_BASIC_CONSISTENCY
  if (a < 0) {
    std::cerr << "Error in Infinitesimal_Ceil\n";
    std::cerr << "calling with a<0 which gives NAN with sqrt\n";
    std::cerr << "a=" << a << "\n";
    throw TerminalException{1};
  }
#endif
  double a_doubl = UniversalTypeConversion<double,T>(a);
  double b_doubl = UniversalTypeConversion<double,T>(b);
  double alpha=-sqrt(a_doubl) - epsilon + b_doubl;
  double eD1=ceil(alpha);
  long int eD2=lround(eD1);
  Tint eReturn=eD2;
  auto f=[&](Tint const& x) -> bool {
    T eDiff=x - b;
    if (eDiff >= 0)
      return true;
    if (eDiff*eDiff <= a)
      return true;
    return false;
  };
  //  std::cerr << "Infinitesimal_ceil, before while loop\n";
  while (true) {
    bool test1=f(eReturn -1);
    bool test2=f(eReturn);
    if (!test1 && test2)
      break;
    if (test1)
      eReturn--;
    if (!test2)
      eReturn++;
  }
  //  std::cerr << "Infinitesimal_ceil, after while loop\n";
  return eReturn - 1;
}






template<typename T, typename Tint>
int computeIt_Kernel(T_shvec_info<T,Tint> & info)
{
  static_assert(is_ring_field<T>::value, "Requires T to be a field");
  int i, j;
  int result = 0;
  int dim = info.request.dim;
  T bound = info.request.bound;
  // The value of bound is assumed to be correct.
  // Thus the Trem values should be strictly positive.
  //  std::cerr << "computeIt : bound=" << bound << "\n";
  MyVector<Tint> Lower(dim);
  MyVector<Tint> Upper(dim);
  MyVector<T> Trem(dim);
  MyVector<T> U(dim);
  MyVector<T> C(dim);
  MyVector<Tint> x(dim);
  MyMatrix<T> g = info.request.gram_matrix;
  /*
  std::cerr << "g=\n";
  for (i=0; i<dim; i++) {
    for (j=0; j<dim; j++)
      std::cerr << " " << g(i,j);
    std::cerr << "\n";
    }*/
  MyMatrix<T> q = info.request.gram_matrix;
  for (i=0; i<dim; i++) {
    for (j=i+1; j<dim; j++) {
      q(i,j) = q(i,j)/q(i,i);
      q(j,i) = q(j,i)/q(i,i);
    }
    for (int i2=i+1; i2<dim; i2++)
      for (int j2=i+1; j2<dim; j2++)
	q(i2,j2)=q(i2,j2) - q(i,i)*q(i,i2)*q(i,j2);
    //    std::cerr << "diag q=" << q(i,i) << "\n";
    /*for (int j=i+1; j<dim; j++)
      std::cerr << "   j=" << j << " q=" << q(i,j) << "\n";
    */
  }
  for (int ik=0; ik<dim; ik++) {
    x(ik)=0;
    Lower(ik)=0;
    Upper(ik)=0;
  }
  bool coset = false;
  i = 0;
  while (i < dim && !coset) {
    coset = (info.request.coset(i) != 0);
    i++;
  }
  bool central=!coset;
  //  std::cerr << "computeIt, coset=" << coset << "\n";
  for (i = 0; i < dim; i++)
    C(i) = info.request.coset(i);
  bool not_finished = true;
  bool needs_new_bound = true;
  i = dim - 1;
  Trem(i) = bound;
  U(i) = 0;
  while (not_finished) {
    //    std::cerr << "Case 1 i=" << i << "\n";
    if (needs_new_bound) {
      T eQuot = Trem(i) / q(i,i);
      T eSum = - U(i) - C(i);
      //      std::cerr << "i=" << i << " eQuot=" << eQuot << "\n";
      Upper(i) = Infinitesimal_Floor<T,Tint>(eQuot, eSum);
      Lower(i) = Infinitesimal_Ceil<T,Tint>(eQuot, eSum);
      x(i) = Lower(i);
      needs_new_bound = false;
    }
    x(i) += 1;
    //    std::cerr << "dim=" << dim << "\n";
    //    std::cerr << "needs_new_bound=" << needs_new_bound << "\n";
    //    std::cerr << "Case 2 i=" << i << "\n";
    //    std::cerr << "vect=";
    //    for (int j=0; j<dim; j++)
    //      std::cerr << " " << x(j);
    //    std::cerr << "\n";
    //    for (int j=0; j<dim; j++) {
    //      std::cerr << "j=" << j << "  L/U=" << Lower(j) << "," << Upper(j) << "\n";
    //    }
    //    std::cerr << "\n";
    //    for (int j=0; j<dim; j++) {
    //      std::cerr << "j=" << j << "  T=" << Trem(j) << "\n";
    //    }
    //    std::cerr << "\n";
    //    std::cerr << "Case 3 i=" << i << "\n";
    if (x(i) <= Upper(i)) {
      //      std::cerr << "Case 4 i=" << i << "\n";
      if (i == 0) {
	//	std::cerr << "Case 5 i=" << i << "\n";
	//	std::cerr << "Insert vector x=";
	//	for (int j3=0; j3<dim; j3++)
	//	  std::cerr << " " << x(j3);
	//	std::cerr << "\n";
	if (central) {
	  j = dim - 1;
	  not_finished = false;
	  while (j >= 0 && !not_finished) {
	    not_finished = (x(j) != 0);
	    j--;
	  }
	  if (!not_finished) return result;
	}
	T hVal=x(0) + C(0) + U(0);
	T eNorm=bound - Trem(0) + q(0,0) * hVal * hVal;
#ifdef CHECK_BASIC_CONSISTENCY
	T norm=0;
	for (int i2=0; i2<dim; i2++)
	  for (int j2=0; j2<dim; j2++)
	    norm += g(i2,j2)*(x(i2) + C(i2)) * (x(j2) + C(j2));
	if (norm != eNorm) {
	  std::cerr << "Norm inconsistency\n";
	  std::cerr << "norm=" << norm << "\n";
	  std::cerr << "eNorm=" << eNorm << "\n";
	  throw TerminalException{1};
	}
	if (eNorm > bound) {
	  std::cerr << "eNorm is too large\n";
	  //	  double bound_doubl, eNorm_doubl, eDiff_doubl;
	  T eDiff=eNorm - bound;
	  double bound_doubl = UniversalTypeConversion<double,T>(bound);
	  double eNorm_doubl = UniversalTypeConversion<double,T>(eNorm);
	  double eDiff_doubl = UniversalTypeConversion<double,T>(eDiff);
	  
	  std::cerr << "bound_doubl=" << bound_doubl << "\n";
	  std::cerr << "eNorm_doubl=" << eNorm_doubl << "\n";
	  std::cerr << "eDiff_doubl=" << eDiff_doubl << "\n";
	  std::cerr << "bound=" << bound << "\n";
	  std::cerr << "eNorm=" << eNorm << "\n";
	  throw TerminalException{1};
	}
#endif

	if (eNorm < info.minimum || info.minimum == -1) {
	  info.short_vectors.clear();
	  info.minimum = eNorm;
	  return TempShvec_globals::STOP_COMPUTATION;
	}
	/*
	std::cerr << "x=";
	for (int u=0; u<dim; u++)
	  std::cerr << " " << x(u);
	  std::cerr << "\n";*/
	info.short_vectors.push_back(x);
	if (central)
	  info.short_vectors.push_back(-x);
      }
      else {
	i--;
	T sum=0;
	for (j = i + 1; j < dim; j++)
	  sum += q(i,j) * (x(j) + C(j));
	U(i) = sum;
	T hVal=x(i+1) + C(i+1) + U(i+1);
	Trem(i) = Trem(i+1) - q(i+1,i+1) * hVal * hVal;
	needs_new_bound = true;
      }
      //      std::cerr << "Case 8 i=" << i << "\n";
    }
    else {
      //      std::cerr << "Case 9 i=" << i << "\n";
      i++;
      if (i >= dim) not_finished = false;
      //      std::cerr << "Case 10 i=" << i << "\n";
    }
    //    std::cerr << "Case 11 i=" << i << "\n";
  }
  std::cerr << "Normal return of computeIt\n";
  return TempShvec_globals::NORMAL_TERMINATION_COMPUTATION;
}


template<typename T, typename Tint>
inline typename std::enable_if<is_ring_field<T>::value,int>::type computeIt(T_shvec_info<T,Tint> & info)
{
  return computeIt_Kernel(info);
}

template<typename T, typename Tint>
inline typename std::enable_if<(not is_ring_field<T>::value),int>::type computeIt(T_shvec_info<T,Tint> & info)
{
  using Tfield=typename overlying_field<T>::field_type;
  //
  T_shvec_request<Tfield> request_field{info.request.dim, info.request.mode,
      UniversalTypeConversion<Tfield,T>(info.request.bound),
      ConvertVectorUniversal<Tfield,T>(info.request.coset),
      ConvertMatrixUniversal<Tfield,T>(info.request.gram_matrix)};
  //
  T_shvec_info<Tfield,Tint> info_field{request_field, info.short_vectors, UniversalTypeConversion<Tfield,T>(info.minimum)};
  int retVal = computeIt_Kernel(info_field);
  info.short_vectors = info_field.short_vectors;
  info.minimum = UniversalTypeConversion<T,Tfield>(info_field.minimum);
  return retVal;
}







template<typename T, typename Tint>
int computeMinimum(T_shvec_info<T,Tint> &info)
{
  int result, coset, i, j;
  int dim = info.request.dim;
  MyVector<T> C(dim);
  //  std::cerr << "Passing by computeMinimum\n";
  coset = 0;
  i = 0;
  while (i < dim && !coset) {
    coset = (info.request.coset(i) != 0);
    i++;
  }
  for (i = 0; i < dim; i++)
    C(i) = info.request.coset(i);
  if (coset) {
    T eNorm=0;
    for (i=0; i<dim; i++)
      for (j=0; j<dim; j++)
	eNorm += info.request.gram_matrix(i,j)*C(i)*C(j);
    info.minimum=eNorm;
  }
  else {
    info.minimum = info.request.gram_matrix(0,0);
    for (i = 1; i < dim; i++)
      if (info.minimum > info.request.gram_matrix(i,i))
	info.minimum = info.request.gram_matrix(i,i);
  }
  while (true) {
    info.request.bound = info.minimum;
    //    double step_size_doubl = UniversalTypeConversion<double,T>(step_size);
    //    double min_doubl = UniversalTypeConversion<double,T>(info.minimum);
    //    std::cerr << "min            =" << info.minimum << "\n";
    //    std::cerr << "min_doubl      =" << min_doubl << "\n";
    //    std::cerr << "step_size      =" << step_size << "\n";
    //    std::cerr << "step_size_doubl=" << step_size_doubl << "\n";
    //    std::cerr << "min_step_size  =" << min_step_size << "\n";
    //    std::cerr << "Before computeIt, case 1\n";
    result = computeIt(info);
    //    std::cerr << "result=" << result << "\n";
    if (result == TempShvec_globals::NORMAL_TERMINATION_COMPUTATION) {
      break;
    }
  }  
  //  std::cerr << "After while loop\n";
  //  std::cerr << "info.minimum=" << info.minimum << "\n";
  //  info.minimum = info.request.bound + step_size;
  //  std::cerr << "Exiting from computeMinimum\n";
  return TempShvec_globals::NORMAL_TERMINATION_COMPUTATION;
}



template<typename T, typename Tint>
void initShvecReq(int dim,
		  MyMatrix<T> const& gram_matrix,
		  T_shvec_info<T,Tint> &info)
{
  if (dim < 2) {
    std::cerr << "shvec: (initShvecReq) wrong input!\n";
    std::cerr << "dimension too low or unassigned matrix\n";
    throw TerminalException{1};
  }
  info.request.dim = dim;
  info.request.coset = MyVector<T>(dim);
  info.request.gram_matrix = gram_matrix;
  info.request.mode = 0;
  info.request.bound = 0;
  info.minimum = -1;
}

template<typename T, typename Tint>
int T_computeShvec(T_shvec_info<T,Tint> &info)
{
  //  int dim = info.request.dim;
  //  std::cerr << "computeShvec mode=" << info.request.mode << "\n";
  switch (info.request.mode)
    {
    case TempShvec_globals::TEMP_SHVEC_MODE_UNDEF:
      {
	std::cerr << "globals::TEMP_SHVEC_MODE is undefined\n";
	throw TerminalException{1};
      }
    case TempShvec_globals::TEMP_SHVEC_MODE_BOUND:
      if (info.request.bound <= 0) {
	std::cerr << "bound=" << info.request.bound << "\n";
	std::cerr << "shvec.c (computeShvec): MODE_BOUND info.request.bound !\n";
	throw TerminalException{1};
      }
      break;
    case TempShvec_globals::TEMP_SHVEC_MODE_SHORTEST_VECTORS:
      if (info.request.bound != 0) {
	std::cerr << "shvec.c (computeShvec): wrong options MODE_SHORTEST_VECTORS!\n";
	throw TerminalException{1};
      }
      break;
    case TempShvec_globals::TEMP_SHVEC_MODE_MINIMUM:
      if (info.request.bound != 0) {
	std::cerr << "shvec.c (computeShvec): wrong options MODE_MINIMUM!\n";
	throw TerminalException{1};
      }
      break;
    case TempShvec_globals::TEMP_SHVEC_MODE_VINBERG:
      break;
    default:
      {
	std::cerr << "shvec.c (computeShvec): wrong options (default)!\n";
	throw TerminalException{1};
      }
    }
  //  std::cerr << "After switch loop\n";
  int result = -99;
  if (info.request.mode == TempShvec_globals::TEMP_SHVEC_MODE_BOUND) {
    //    std::cerr << "Before computeIt, case 2\n";
    result = computeIt<T,Tint>(info);
  }
  else if (info.request.mode == TempShvec_globals::TEMP_SHVEC_MODE_SHORTEST_VECTORS) {
    //    std::cerr << "Before computeMinimum\n";
    result = computeMinimum<T,Tint>(info);
    //    info.request.bound = info.minimum;
    //    std::cerr << "Assign info.request.bound\n";
    //    std::cerr << "Before computeIt, case 3\n";
    //    result = computeIt<T,Tint>(info);
  }
  else if (info.request.mode == TempShvec_globals::TEMP_SHVEC_MODE_MINIMUM) {
    //    std::cerr << "Before computeIt, case 4\n";
    result = computeMinimum<T,Tint>(info);
  }
  else if (info.request.mode == TempShvec_globals::TEMP_SHVEC_MODE_VINBERG) {
    //    std::cerr << "Before computeIt, case 5\n";
    result = computeIt<T,Tint>(info);
  }
  return result;
  //  std::cerr << "result=" << result << "\n";
}


template<typename T,typename Tint>
resultCVP<T,Tint> CVPVallentinProgram_exact(MyMatrix<T> const& GramMat, MyVector<T> const& eV)
{
  int dim=GramMat.rows();
  MyVector<T> cosetVect(dim);
  //  std::cerr << "coset=";
  for (int i=0; i<dim; i++) {
    cosetVect(i)=-eV(i);
    //    std::cerr << " " << cosetVect(i);
  }
  //  std::cerr << "\n";
  if (IsIntegralVector(eV)) {
    T TheNorm=0;
    MyMatrix<Tint> ListVect(1,dim);
    for (int i=0; i<dim; i++)
      ListVect(0,i)=UniversalTypeConversion<Tint,T>(eV(i));
    return {TheNorm, ListVect};
  }
  T bound=0; // should not be used
  int mode = TempShvec_globals::TEMP_SHVEC_MODE_SHORTEST_VECTORS;
  T_shvec_info<T,Tint> info;
  initShvecReq<T>(dim, GramMat, info);
  info.request.bound = bound;
  info.request.mode = mode;
  info.request.coset = cosetVect;
  info.minimum = -44;
  //  std::cerr << "Before T_computeShvec\n";
  int result=T_computeShvec(info);
  //  std::cerr << "After T_computeShvec\n";
  if (result == -2444) {
    std::cerr << "Probably wrong\n";
    throw TerminalException{1};
  }
  int nbVect=info.short_vectors.size();
  MyMatrix<Tint> ListClos(nbVect, dim);
  for (int iVect=0; iVect<nbVect; iVect++)
    for (int i=0; i<dim; i++)
      ListClos(iVect,i)=info.short_vectors[iVect](i);
  MyVector<T> eDiff(dim);
  for (int i=0; i<dim; i++)
    eDiff(i)=ListClos(0,i) - eV(i);
  T TheNorm=EvaluationQuadForm<T,T>(GramMat, eDiff);
  return {TheNorm, ListClos};
}



template<typename T, typename Tint>
MyMatrix<Tint> T_ShortVector_exact(MyMatrix<T> const& GramMat, T const&MaxNorm)
{
  int dim=GramMat.rows();
  if (dim == 1) {
    std::vector<MyVector<Tint>> ListVect;
    MyVector<Tint> eVect(1);
    eVect(0)=0;
    ListVect.push_back(eVect);
    int idx=1;
    while (true) {
      T norm = idx*idx * GramMat(0,0);
      if (norm > MaxNorm)
	break;
      MyVector<Tint> eVect1(1), eVect2(1);
      eVect1(0)=-idx;
      eVect2(0)= idx;
      ListVect.push_back(eVect1);
      ListVect.push_back(eVect2);
      idx++;
    }
    return MatrixFromVectorFamily(ListVect);
  }
  T bound=MaxNorm;
  int mode = TempShvec_globals::TEMP_SHVEC_MODE_BOUND;
  MyVector<T> cosetVect=ZeroVector<T>(dim);
  T_shvec_info<T,Tint> info;
  initShvecReq<T>(dim, GramMat, info);
  info.request.bound = bound;
  info.request.mode = mode;
  info.request.coset = cosetVect;
  info.minimum = -44;
  //
  int result=T_computeShvec(info);
  if (result == -2444) {
    std::cerr << "Probably wrong\n";
    throw TerminalException{1};
  }
  //
  return MatrixFromVectorFamily(info.short_vectors);
}






#endif

#ifndef INCLUDE_TEMPLATE_SHVEC_H
#define INCLUDE_TEMPLATE_SHVEC_H


#include "LatticeDefinitions.h"
#include "NumberTheory.h"
#include "MAT_Matrix.h"


namespace TempShvec_globals {
  const int TEMP_SHVEC_MODE_UNDEF=-1;
  const int TEMP_SHVEC_MODE_BOUND=0;
  const int TEMP_SHVEC_MODE_SHORTEST_VECTORS=1;
  const int TEMP_SHVEC_MODE_MINIMUM=2;
  const int TEMP_SHVEC_MODE_THETA_SERIES=3;
  const int TEMP_SHVEC_MODE_VINBERG=4;
  const int STOP_COMPUTATION=666;
}

template<typename T>
struct T_shvec_request {
  int dim;
  int mode;
  T bound;
  int number;
  MyVector<T> coset;
  MyMatrix<T> gram_matrix;
};

template<typename T>
struct T_shvec_info {
  T_shvec_request<T> request;
  int short_vectors_count;
  int short_vectors_number;
  std::vector<MyVector<int>> short_vectors;
  T minimum;
};


template<typename T, typename Tint>
int insertBound(T_shvec_info<T> &info,
		MyVector<Tint> const& vector, 
		int coset,
		T norm)
{
  info.short_vectors.push_back(vector);
  //  std::cerr << "FUNCTION insertBound\n";
  info.short_vectors_count++;
  info.short_vectors_number++;
  //  std::cerr << "short_vectors_number=" << info.short_vectors_number << "\n";
  
  if (norm < info.minimum || info.minimum == -1)
    info.minimum = norm;

  if (info.short_vectors_number == info.request.number)
    return TempShvec_globals::STOP_COMPUTATION;
  else
    return 0;
}

// We return
// floor(sqrt(A) + epsilon + B) 
template<typename T>
int Infinitesimal_Floor_V1(T const& a, T const& b)
{
  double a_doubl, b_doubl;
  double epsilon=0.000000001;
  GET_DOUBLE(a, a_doubl);
  GET_DOUBLE(b, b_doubl);
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
  double a_doubl, b_doubl;
  double epsilon=0.000000001;
  GET_DOUBLE(a, a_doubl);
  GET_DOUBLE(b, b_doubl);
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
template<typename T>
Tint Infinitesimal_Floor(T const& a, T const& b)
{
  double a_doubl, b_doubl;
  double epsilon=0.000000001;
  GET_DOUBLE(a, a_doubl);
  GET_DOUBLE(b, b_doubl);
  double alpha=sqrt(a_doubl) + epsilon + b_doubl;
  double eD1=floor(alpha);
  long int eD2=lround(eD1);
  Tint eReturn=eD2;
  auto f=[&](Tint const& x) -> bool {
    T eDiff=x - b;
    if (eDiff <= 0)
      return true;
    if (eDiff*eDiff <= a)
      return true;
    return false;
  };
  //  std::cerr << "a=" << a << "\n";
  //  std::cerr << "b=" << b << "\n";
  //  std::cerr << "Infinitesimal_floor, before while loop\n";
  while(true) {
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
template<typename T>
Tint Infinitesimal_Ceil(T const& a, T const& b)
{
  double a_doubl, b_doubl;
  double epsilon=0.000000001;
  GET_DOUBLE(a, a_doubl);
  GET_DOUBLE(b, b_doubl);
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
  while(true) {
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
int insertStop(T_shvec_info<T> &info,
	       MyVector<Tint> const& vector, 
	       int coset,
	       T norm)
{
  //  std::cerr << "Assignation of info.minimum\n";
  //  std::cerr << "norm=" << norm << "\n";
  info.minimum = norm;
  return TempShvec_globals::STOP_COMPUTATION;
}

template<typename T, typename Tint>
int computeIt(T_shvec_info<T> &info,
	      int (*insert)(T_shvec_info<T> &info,
			    MyVector<Tint> const& vector, 
			    int coset,
			    T norm))
{
  int coset, result, not_finished, i, j, dim;
  bool needs_new_bound;
  T sum, Z, bound;
  //  double eQuot_doubl;
  result = 0;
  dim = info.request.dim;
  bound = info.request.bound;
  //  std::cerr << "bound=" << bound << "\n";
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
  coset = 0;
  i = 0;
  while (i < dim && !coset) {
    coset = (info.request.coset(i) != 0.0);
    i++;
  }
  for (i = 0; i < dim; i++)
    C(i) = info.request.coset(i);
  not_finished = 1;
  needs_new_bound = true;
  i = dim - 1;
  Trem(i) = bound;
  U(i) = 0;
  while (not_finished) {
    //    std::cerr << "Case 1 i=" << i << "\n";
    if (needs_new_bound) {
      T eQuot = Trem(i) / q(i,i);
      T eSum = - U(i) - C(i);
      //      GET_DOUBLE(eQuot, eQuot_doubl);
      //      std::cerr << "eQuot_doubl=" << eQuot_doubl << "\n";
      Upper(i) = Infinitesimal_Floor(eQuot, eSum);
      Lower(i) = Infinitesimal_Ceil(eQuot, eSum);
      x(i) = Lower(i);
      needs_new_bound = false;
    }
    x(i)=x(i)+1;
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
	if (!coset) {
	  j = dim - 1;
	  not_finished = 0;
	  while (j >= 0 && !not_finished) {
	    not_finished = (x(j) != 0);
	    j--;
	  }
	  if (!not_finished) return result;
	}
	T hVal=x(0) + C(0) + U(0);
	T eNorm=bound - Trem(0) + q(0,0) * hVal * hVal;
	T norm=0;
	for (int i2=0; i2<dim; i2++)
	  for (int j2=0; j2<dim; j2++)
	    norm=norm + g(i2,j2)*(x(i2) + C(i2)) * (x(j2) + C(j2));
	if (norm != eNorm) {
	  std::cerr << "Norm inconsistency\n";
	  std::cerr << "norm=" << norm << "\n";
	  std::cerr << "eNorm=" << eNorm << "\n";
	  throw TerminalException{1};
	}
	if (eNorm > bound) {
	  std::cerr << "eNorm is too large\n";
	  double bound_doubl, eNorm_doubl, eDiff_doubl;
	  T eDiff=eNorm - bound;
	  GET_DOUBLE(bound, bound_doubl);
	  GET_DOUBLE(eNorm, eNorm_doubl);
	  GET_DOUBLE(eDiff, eDiff_doubl);
	  
	  std::cerr << "bound_doubl=" << bound_doubl << "\n";
	  std::cerr << "eNorm_doubl=" << eNorm_doubl << "\n";
	  std::cerr << "eDiff_doubl=" << eDiff_doubl << "\n";
	  std::cerr << "bound=" << bound << "\n";
	  std::cerr << "eNorm=" << eNorm << "\n";
	  throw TerminalException{1};
	}

	//	std::cerr << "Case 6 i=" << i << "\n";
	result = insert(info, x, coset, eNorm);
	//	std::cerr << "Case 7 i=" << i << "\n";
	//	std::cerr << "result=" << result << "\n";
	if (result == TempShvec_globals::STOP_COMPUTATION) {
	  //	  std::cerr << "Before return statement\n";
	  return result;
	}
	//	std::cerr << "Case 7.1 i=" << i << "\n";
      }
      else {
	i--;
	sum=0;
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
      if (i >= dim) not_finished = 0;
      //      std::cerr << "Case 10 i=" << i << "\n";
    }
    //    std::cerr << "Case 11 i=" << i << "\n";
  }
  return 0;
}

template<typename T>
int computeMinimum(T_shvec_info<T> &info)
{
  int result, dim, coset, i, j;
  dim = info.request.dim;
  MyVector<T> C(dim);
  //  std::cerr << "Passing by computeMinimum\n";
  coset = 0;
  i = 0;
  while (i < dim && !coset) {
    coset = (info.request.coset(i) != 0.0);
    i++;
  }
  for (i = 0; i < dim; i++)
    C(i) = info.request.coset(i);
  if (coset) {
    T eNorm=0;
    for (i=0; i<dim; i++)
      for (j=0; j<dim; j++)
	eNorm = eNorm + info.request.gram_matrix(i,j)*C(i)*C(j);
    info.minimum=eNorm;
  }
  else {
    info.minimum = info.request.gram_matrix(0,0);
    for (i = 1; i < dim; i++)
      if (info.minimum > info.request.gram_matrix(i,i))
	info.minimum = info.request.gram_matrix(i,i);
  }
  T eNum=1;
  T eDen=10000;
  T min_step_size=eNum/eDen;
  T step_size=info.minimum;
  while(true) {
    info.request.bound = info.minimum - step_size;
    double step_size_doubl;
    double min_doubl;
    GET_DOUBLE(info.minimum, min_doubl);
    GET_DOUBLE(step_size, step_size_doubl);
    //    std::cerr << "min            =" << info.minimum << "\n";
    //    std::cerr << "min_doubl      =" << min_doubl << "\n";
    //    std::cerr << "step_size      =" << step_size << "\n";
    //    std::cerr << "step_size_doubl=" << step_size_doubl << "\n";
    //    std::cerr << "min_step_size  =" << min_step_size << "\n";
    //    std::cerr << "Before computeIt, case 1\n";
    result = computeIt(info, &insertStop);
    //    std::cerr << "result=" << result << "\n";
    if (result != TempShvec_globals::STOP_COMPUTATION) {
      if (step_size < min_step_size)
	break;
      step_size=step_size/2;
    }
  }  
  //  std::cerr << "After while loop\n";
  //  std::cerr << "info.minimum=" << info.minimum << "\n";
  info.minimum = info.request.bound + step_size;
  return 0;
}



template<typename T>
void initShvecReq(int dim,
		  MyMatrix<T> const& gram_matrix,
		  T_shvec_info<T> &info)
{
  if (dim < 2) {
    std::cerr << "shvec: (initShvecReq) wrong input!\n";
    std::cerr << "dimension too low or unassigned matrix\n";
    throw TerminalException{1};
  }
  info.request.dim = dim;
  info.request.coset = MyVector<T>(dim);
  info.request.gram_matrix=gram_matrix;
  info.request.mode = 0;
  info.request.bound = 0.0;
  info.request.number = 0;
  info.short_vectors_count = 0;
  info.short_vectors_number = 0;
  info.minimum = -1;
}

template<typename T>
int T_computeShvec(T_shvec_info<T> &info)
{
  int dim, coset, i;
  dim = info.request.dim;
  coset = 0;
  i = 0;
  while (i < dim && !coset) {
    coset = (info.request.coset(i) != 0.0);
    i++;
  }
  //  std::cerr << "computeShvec mode=" << info.request.mode << "\n";
  switch (info.request.mode)
    {
    case TempShvec_globals::TEMP_SHVEC_MODE_UNDEF:
      {
	std::cerr << "globals::TEMP_SHVEC_MODE is undefined\n";
	throw TerminalException{1};
      }
    case TempShvec_globals::TEMP_SHVEC_MODE_BOUND:
      if (info.request.bound <= 0.0)
	{
	  std::cerr << "bound=" << info.request.bound << "\n";
	  std::cerr << "shvec.c (computeShvec): MODE_BOUND info.request.bound !\n";
	  throw TerminalException{1};
	}
      if (info.request.number < 0) {
	std::cerr << "shvec.c (computeShvec): MODE_BOUND info.request.number!\n";
	throw TerminalException{1};
      }
      break;
    case TempShvec_globals::TEMP_SHVEC_MODE_SHORTEST_VECTORS:
      if (info.request.bound != 0.0 || info.request.number < 0) {
	std::cerr << "shvec.c (computeShvec): wrong options MODE_SHORTEST_VECTORS!\n";
	throw TerminalException{1};
      }
      break;
    case TempShvec_globals::TEMP_SHVEC_MODE_MINIMUM:
      if (info.request.bound != 0.0 || info.request.number != 0) {
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
  int result = -99;
  if (info.request.mode == TempShvec_globals::TEMP_SHVEC_MODE_BOUND) {
    std::cerr << "Before computeIt, case 2\n";
    result = computeIt(info, &insertBound);
  }
  else if (info.request.mode == TempShvec_globals::TEMP_SHVEC_MODE_SHORTEST_VECTORS) {
    result = computeMinimum(info);
    info.request.bound = info.minimum;
    //    std::cerr << "Assign info.request.bound\n";
    //    std::cerr << "Before computeIt, case 3\n";
    result = computeIt(info, &insertBound);
  }
  else if (info.request.mode == TempShvec_globals::TEMP_SHVEC_MODE_MINIMUM) {
    result = computeMinimum(info);
  }
  else if (info.request.mode == TempShvec_globals::TEMP_SHVEC_MODE_VINBERG) {
    result = computeIt(info, &insertBound);
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
  int number=0;
  T_shvec_info<T> info;
  initShvecReq<T>(dim, GramMat, info);
  info.request.bound = bound;
  info.request.mode = mode;
  info.request.number = number;
  info.request.coset = cosetVect;
  info.minimum = -44;
  //  std::cerr << "Before T_computeShvec\n";
  int result=T_computeShvec(info);
  //  std::cerr << "After T_computeShvec\n";
  if (result == -2444) {
    std::cerr << "Probably wrong\n";
    throw TerminalException{1};
  }
  int nbVect=info.short_vectors_number;
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
    while(1) {
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
  int number=0;
  MyVector<T> cosetVect=ZeroVector<T>(dim);
  T_shvec_info<T> info;
  initShvecReq<T>(dim, GramMat, info);
  info.request.bound = bound;
  info.request.mode = mode;
  info.request.number = number;
  info.request.coset = cosetVect;
  info.minimum = -44;
  //
  int result=T_computeShvec(info);
  if (result == -2444) {
    std::cerr << "Probably wrong\n";
    throw TerminalException{1};
  }
  //
  std::vector<MyVector<Tint>> ListVect;
  for (auto & eVect : info.short_vectors) {
    MyVector<Tint> eVectImg = ConvertVectorUniversal<Tint,int>(eVect);
    ListVect.push_back(eVectImg);
  }
  return MatrixFromVectorFamily(ListVect);
}






#endif

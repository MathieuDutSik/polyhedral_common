#ifndef CTYPE_MPI_TYPES
#define CTYPE_MPI_TYPES


#include "MAT_Matrix.h"

struct TypeIndex {
  int iProc;
  int idxMatrix;
  int iAdj;
};


template<typename T>
struct TypeCtypeExch {
  MyMatrix<T> eMat;
};


template<typename T>
struct PairExch {
  TypeCtypeExch<T> eCtype;
  TypeIndex eIndex;
};


template<typename T>
std::ostream& operator<<(std::ostream& os, TypeCtypeExch<T> const& obj)
{
  int nbRow=obj.eMat.rows();
  int nbCol=obj.eMat.cols();
  os << " " << nbRow << " " << nbCol;
  for (int iRow=0; iRow<nbRow; iRow++)
    for (int iCol=0; iCol<nbCol; iCol++)
      os << " " << obj.eMat(iRow, iCol);
  return os;
}


template<typename T>
bool operator==(TypeCtypeExch<T> const& obj1, TypeCtypeExch<T> const& obj2)
{
  int nbRow1=obj1.eMat.rows();
  int nbRow2=obj2.eMat.rows();
  if (nbRow1 != nbRow2)
    return false;
  int nbCol1=obj1.eMat.cols();
  int nbCol2=obj2.eMat.cols();
  if (nbCol1 != nbCol2)
    return false;
  for (int iRow=0; iRow<nbRow1; iRow++)
    for (int iCol=0; iCol<nbCol1; iCol++)
      if (obj1.eMat(iRow, iCol) != obj2.eMat(iRow, iCol))
        return false;
  return true;
}



std::ostream& operator<<(std::ostream& os, TypeIndex const& obj)
{
  os << obj.iProc << " " << obj.idxMatrix << " " << obj.iAdj;
  return os;
}






namespace std {
  template<typename T>
  struct less<TypeCtypeExch<T>> {
    bool operator()(TypeCtypeExch<T> const& eTPE1, TypeCtypeExch<T> const& eTPE2) const
    {
      int nbRow=eTPE1.eMat.rows();
      int nbCol=eTPE1.eMat.cols();
      for (int iRow=0; iRow<nbRow; iRow++)
        for (int iCol=0; iCol<nbCol; iCol++) {
          if (eTPE1.eMat(iRow,iCol) < eTPE2.eMat(iRow,iCol))
            return true;
          if (eTPE1.eMat(iRow,iCol) > eTPE2.eMat(iRow,iCol))
            return false;
        }
      return false;
    }
  };
}




namespace boost { namespace serialization {
    // TypeCtypeExch
    template<class Archive, typename T>
      inline void serialize(Archive & ar,
                            TypeCtypeExch<T> & eRecMat,
                            const unsigned int version)
      {
        int rows = eRecMat.eMat.rows();
        int cols = eRecMat.eMat.cols();
        ar & make_nvp("rows", rows);
        ar & make_nvp("cols", cols);
        eRecMat.eMat.resize(rows, cols);
        for (int r = 0; r < rows; ++r)
          for (int c = 0; c < cols; ++c)
            ar & make_nvp("val", eRecMat.eMat(r,c));
      }

    // TypeCtypeExch
    template<class Archive>
      inline void serialize(Archive & ar,
                            TypeIndex & eTypIdx,
                            const unsigned int version)
      {
        ar & make_nvp("iProc", eTypIdx.iProc);
        ar & make_nvp("idxMatrix", eTypIdx.idxMatrix);
        ar & make_nvp("iAdj", eTypIdx.iAdj);
      }

    // PairExch
    template<class Archive, typename T>
      inline void serialize(Archive & ar,
                            PairExch<T> & ePair,
                            const unsigned int version)
      {
        ar & make_nvp("perfect", ePair.eCtype);
        ar & make_nvp("index"  , ePair.eIndex);
      }

  }}

namespace std {
  template <typename Tint>
  struct hash<TypeCtypeExch<Tint>>
  {
    std::size_t operator()(const TypeCtypeExch<Tint>& k) const
    {
      std::size_t h1 = 0;
      int nbRow=k.eMat.rows();
      int nbCol=k.eMat.cols();
      for (int iRow=0; iRow<nbRow; iRow++)
        for (int iCol=0; iCol<nbCol; iCol++) {
          Tint eVal = k.eMat(iRow,iCol);
          std::size_t h2 = std::hash<Tint>()(eVal);
          h1 = h2 ^ ( h1 << 1);
        }
      return h1;
    }
  };
}




template<typename T>
TypeCtypeExch<T> ParseStringToCtypeExch(std::string const& str)
{
  std::vector<std::string> LStr = STRING_Split(str, " ");
  int nRow;
  std::istringstream(LStr[0]) >> nRow;
  int n;
  std::istringstream(LStr[1]) >> n;
  MyMatrix<T> eMat(nRow,n);
  int idx=2;
  for (int iRow=0; iRow<nRow; iRow++) {
    for (int iCol=0; iCol<n; iCol++) {
      T eVal;
      std::istringstream(LStr[idx]) >> eVal;
      eMat(iRow,iCol) = eVal;
      idx++;
    }
  }
  return {eMat};
}

TypeIndex ParseStringToTypeIndex(std::string const& str)
{
  std::vector<std::string> LStr = STRING_Split(str, " ");
  int iProc;
  std::istringstream(LStr[0]) >> iProc;
  int idxMatrixF;
  std::istringstream(LStr[1]) >> idxMatrixF;
  int iAdj;
  std::istringstream(LStr[1]) >> iAdj;
  return {iProc, idxMatrixF, iAdj};
}


#endif

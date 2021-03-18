#ifndef INCLUDE_POLY_RECURSIVE_DUAL_DESC_H
#define INCLUDE_POLY_RECURSIVE_DUAL_DESC_H


namespace datagap {

static const int int_scalar = 1;
static const int int_string = 2;
static const int int_permutation = 3;
static const int int_group = 4;
static const int int_list  = 5;
static const int int_record = 6;



template<typename T, typename Telt>
struct DataGAP {
  // 1: for scalar value
  // 2: for string
  // 3: for permutation
  // 4: for Group
  // 5: for List
  // 6: for record
  int Nature;
  // 1: value
  T scalar;
  // 2: for a string
  std::string string;
  // 3: permutation
  Telt permutation;
  // 4, 5: List of entries, which applies for group as well.
  std::vector<DataGAP> ListEnt:
  // 5: The entries of record
  std::vector<std::pair<std::string, DataGAP>> ListRec;
};




std::vector<std::string_view> ParseStringByComma(std::string_view const& estr)
{
  int LevelParenthesis=0;
  int LevelBracket=0;
  size_t n_char=estr.size();
  size_t pos_start = 0;
  std::vector<std::string_view> LStr;
  auto insert=[&](size_t const& pos1, size_t const& pos2) -> void {
    size_t len = pos2 - pos1;
    LStr.push_back(estr.substr(pos_start, len));
    pos_start = pos2 + 1;
  };
  for (size_t i_char=0; i_char<n_char; i_char++) {
    std::string echar = estr.substr(i_char, 1);
    if (echar == "(")
      LevelParenthesis++;
    if (echar == ")")
      LevelParenthesis--;
    if (echar == "[")
      LevelBracket++;
    if (echar == "]")
      LevelBracket--;
    if (LevelParenthesis == 0 && LevelBracket == 0 && echar == ",")
      insert(pos_start, i_char);
  }
  insert(pos_start, n_char);
  return LStr;
}


// The recursive function for parsing the entries of the text file.
// It should not have any space in it

template<typename T, typename Telt>
DataGAP ParseGAPString(std::string_view const& full_str)
{
  // Case 2: a string
  if (full_str.substr(0,1) == "\"") {
    size_t n_char = full_str.size();
    if (full_str.substr(n_char-1 , 1) != "\"") {
      std::cerr << "Parsing error for the string\n";
      throw TerminalException{1};
    }
    std::string str = full_str.substr(1, n_char-2);
    return {int_string, {}, str, {}, {}, {}};
  }
  // Case 5: a list
  if (full_str.substr(0,1) == "[") {
    size_t n_char = full_str.size();
    if (full_str.substr(n_char-1 , 1) != "]") {
      std::cerr << "Parsing error for the string\n";
      throw TerminalException{1};
    }
    std::vector<std::string_view> LStr = ParseStringByComma(full_str.substr(1,n_char-2));
    std::vector<DataGAP> LVal;
    for (auto & estr : LStr) {
      LVal.push_back(ParseGAPString(estr));
    }
    return {int_list, {}, {}, {}, LVal, {}};
  }
  // Case 4: a group
  if (full_str.substr(0,1) == "G") {
    size_t n_char = full_str.size();
    if (full_str.substr(0 , 7) != "Group([") {
      std::cerr << "Parsing error for the group 1\n";
      throw TerminalException{1};
    }
    if (full_str.substr(n_char-2 , 2) != "])") {
      std::cerr << "Parsing error for the group 2\n";
      throw TerminalException{1};
    }
    std::vector<std::string_view> LStr = ParseStringByComma(full_str.substr(1,n_char-2));
    std::vector<DataGAP> LVal;
    for (auto & estr : LStr) {
      LVal.push_back(ParseGAPString(estr));
    }
    return {int_group, {}, {}, {}, LVal, {}};
  }
  // Case 6: the record
  if (full_str.substr(0,1) == "r") {
    size_t n_char = full_str.size();
    if (full_str.substr(0 , 4) != "rec(") {
      std::cerr << "Parsing error for the record 1\n";
      throw TerminalException{1};
    }
    if (full_str.substr(n_char-1 , 2) != ")") {
      std::cerr << "Parsing error for the record 2\n";
      throw TerminalException{1};
    }
    std::vector<std::string_view> LStr = ParseStringByComma(full_str.substr(1,n_char-2));
    std::vector<std::pair<std::string, DataGAP>> LVal;
    for (auto & estr : LStr) {
      size_t pos = estr.find(":=");
      std::string name = estr.substr(0, pos);
      std::string_view sstr = estr.substr(pos+2, estr.size() - 2 - pos);
      DataGAP eEnt = ParseGAPString(sstr);
      LVal.push_back({name, eEnt});
    }
    return {int_record, {}, {}, {}, {}, LVal};
  }
  // Case 3: the permutation case
  if (full_str.substr(0,1) == "(") {
    Telt g = ParsePermutation(full_str);
    return {int_permutation, {}, {}, g, {}, {}};
  }
  // Case 1: the element
  T scalar = ParseScalar<T>(full_str);
  return {int_scalar, scalar, {}, {}, {}, {}};
}



template<typename T, typename Telt>
DataGAP ParseGAPString(std::string const& eFile)
{
  std::ifstream is(eFile);
  std::string line;
  std::string full_str;
  auto append_content=[&](std::string const& ustr) {
    size_t n_char = ustr.size();
    for (size_t i_char=0; i_char<n_char; i_char++) {
      std::string echar = ustr.substr(i_char, 1);
      if (echar != " ")
        full_str += echar;
    }
  };
  size_t iLine=0;
  while (std::getline(is, line)) {
    if (iLine == 0) {
      if (line.size() < 8) {
        std::cerr << "The first line should begin by \"return ....\"";
        throw TerminalException{1};
      }
      if (line.substr(0, 7) != "return") {
        std::cerr << "The first line should begin by \"return ....\"";
        throw TerminalException{1};
      }
      append_content(line.substr(7, line.size() - 7));
    } else {
      append_content(line);
    }
    iLine++;
  }
  size_t n_char = full_str.size();
  if (full_str.substr(n_char-1, 1) != ";") {
    std::cerr << "The lčast character is not a semicolon ;. Wrong input\n";
    throw TerminalException{1};
  }
  std::string_view full_view = full_str.substr(0, n_char-1);
  return ParseGAPString(full_view);
}



template<typename T, typename Telt>
T ConvertGAPread_ScalarT(DataGAP<T,Telt> const& data)
{
  if (data.Nature != int_scalar) {
    std::cerr << "It should be a scalar for effective conversion to scalar\n";
    throw TerminalException{1};
  }
  return data.scalar;
}


template<typename T, typename Telt>
MyMatrix<T> ConvertGAPread_MyVectorT(DataGAP<T,Telt> const& data)
{
  if (data.Nature != int_list) {
    std::cerr << "It should be a list for effective conversion to MyVector\n";
    throw TerminalException{1};
  }
  std::vector<T> ListVal;
  for (auto & eScal : data.ListEnt)
    ListVal.push_back(ConvertGAPread_ScalarT(eScal));
  size_t len = ListVal.size();
  MyVector<T> V(len);
  for (size_t i=0; i<len; i++)
    V(i) = ListVal[i];
  return V;
}



template<typename T, typename Telt>
MyMatrix<T> ConvertGAPread_MyMatrixT(DataGAP<T,Telt> const& data)
{
  if (data.Nature != int_list) {
    std::cerr << "It should be a list for effective conversion to MyMatrix\n";
    throw TerminalException{1};
  }
  std::vector<MyVector<T>> ListV;
  for (auto & eVect : data.ListEnt)
    ListV.push_back(ConvertGAPread_MyVectorT(eVect));
  return MatrixFromVectorFamily(ListV);
}



template<typename T, typename Telt>
Telt ConvertGAPread_Permutation(DataGAP<T,Telt> const& data)
{
  if (data.Nature != int_permutation) {
    std::cerr << "It should be a permutation for effective conversion to permutation\n";
    throw TerminalException{1};
  }
  return data.permutation;
}




template<typename T, typename Tgroup>
Tgroup ConvertGAPread_PermutationGroup(DataGAP<T,Tgroup::Telt> const& data, int const& n)
{
  if (data.Nature != int_group) {
    std::cerr << "It should be a group for effective conversion to MyMatrix\n";
    throw TerminalException{1};
  }
  std::vector<Telt> ListGen;
  for (auto & eEnt : data.ListEnt) {
    Telt g1 = ConvertGAPread_Permutation(eEnt);
    Telt g2 = ExtendPermutation(g1, n);
    ListGen.push_back(g2);
  }
  return Tgroup(ListGen, n);
}



template<typename T, typename Telt>
Face ConvertGAPread_Face(DataGAP<T,Telt> const& data, int const& n)
{
  if (data.Nature != int_list) {
    std::cerr << "It should be a list for effective conversion to face\n";
    throw TerminalException{1};
  }
  Face f(n);
  for (auto & eEnt : data.ListEnt) {
    T val_T = ConvertGAPread_ScalarT(eEnt);
    int val_i = UniversalTypeConversion<int,T>(val_T);
    f[val_i - 1] = 1;
  }
  return f;
}



template<typename T, typename Telt>
std::vector<Face> ConvertGAPread_ListFace(DataGAP<T,Telt> const& data, int const& n)
{
  if (data.Nature != int_list) {
    std::cerr << "It should be a list for effective conversion to face\n";
    throw TerminalException{1};
  }
  std::vector<Face> ListFace;
  for (auto & eEnt : data.ListEnt)
    ListFace.push_back(ConvertGAPread_Face(eEnt, n));
  return ListFace;
}










}
#endif

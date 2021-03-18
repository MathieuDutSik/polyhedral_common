#ifndef INCLUDE_POLY_RECURSIVE_DUAL_DESC_H
#define INCLUDE_POLY_RECURSIVE_DUAL_DESC_H



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
    if (echar == "(") {
      LevelParenthesis++;
    }
    if (echar == ")") {
      LevelParenthesis--;
    }
    if (echar == "[") {
      LevelBracket++;
    }
    if (echar == "]") {
      LevelBracket--;
    }
    if (LevelParenthesis == 0 && LevelBracket == 0 && echar == ",") {
      size_t pos_end = i_char;
      insert(pos_start, pos_end);
    }
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
    return {2, {}, str, {}, {}, {}};
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
    return {5, {}, {}, {}, LVal, {}};
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
    return {4, {}, {}, {}, LVal, {}};
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
    return {6, {}, {}, {}, {}, LVal};
  }
  // Case 3: the permutation case
  if (full_str.substr(0,1) == "(") {
    Telt g = ParsePermutation(full_str);
    return {3, {}, {}, g, {}, {}};
  }
  // Case 1: the element
}







#endif

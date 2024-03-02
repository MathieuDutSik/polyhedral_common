// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_DUALDESC_BASIC_DATAFILE_H_
#define SRC_DUALDESC_BASIC_DATAFILE_H_

#include "Basic_file.h"
#include "Boost_bitset.h"
#include <algorithm>
#include <string>
#include <vector>

#ifdef DEBUG
#define DEBUG_BASIC_DATAFILE
#endif

//
// Storing a single number to file
//

struct FileNumber {
private:
  std::FILE *fp;
  std::string file;

public:
  FileNumber(const FileNumber &) = delete;
  FileNumber &operator=(const FileNumber &) = delete;
  FileNumber(FileNumber &&) = delete;
  FileNumber() = delete;

  FileNumber(std::string const &file, bool overwrite) : file(file) {
    if (overwrite) {
      fp = std::fopen(file.data(), "w+");
    } else {
      if (!IsExistingFile(file)) {
        std::cerr << "FileNumber: The file " << file
                  << " should not be missing\n";
        throw TerminalException{1};
      }
      fp = std::fopen(file.data(), "r+");
    }
  }

  ~FileNumber() {
    //    std::cerr << "Clean destruction of FileNumber associated to file=" <<
    //    file << "\n";
    std::fclose(fp);
  }

  size_t getval() {
    size_t i_byte = 0;
    std::fseek(fp, i_byte, SEEK_SET);
    size_t val;
    size_t ret = std::fread(&val, sizeof(size_t), 1, fp);
    if (ret != 1) {
      std::cerr << "FileNumber: Error in getval()\n";
      std::cerr << "val=" << val << " ret=" << ret << " file=" << file << "\n";
      throw TerminalException{1};
    }
    return val;
  }

  void setval(size_t const &val) {
    size_t i_byte = 0;
    std::fseek(fp, i_byte, SEEK_SET);
    size_t n_write = std::fwrite(&val, sizeof(size_t), 1, fp);
    if (n_write != 1) {
      std::cerr << "Error in setval(..), wrong number of written bytes\n";
      throw TerminalException{1};
    }
  }
};

size_t FileNumber_Read(std::string const& FileNb) {
  bool overwrite = false;
  FileNumber fn(FileNb, overwrite);
  return fn.getval();
}

void FileNumber_Write(std::string const& FileNb, size_t const& val) {
  bool overwrite = true;
  FileNumber fn(FileNb, overwrite);
  fn.setval(val);
}

//
// Storing a sequence of bits that can be extended
//

struct FileBool {
private:
  std::FILE *fp;
  size_t n_ent;
  std::string file;

public:
  FileBool(const FileBool &) = delete;
  FileBool &operator=(const FileBool &) = delete;
  FileBool(FileBool &&) = delete;
  FileBool() = delete;

  FileBool(std::string const &file) : file(file) {
    fp = std::fopen(file.data(), "w+");
    n_ent = 0;
  }

  FileBool(std::string const &file, size_t const &_n_ent) : file(file) {
    if (!IsExistingFile(file)) {
      std::cerr << "FileBool: The file " << file << " should not be missing\n";
      throw TerminalException{1};
    }
    fp = std::fopen(file.data(), "r+");
    n_ent = _n_ent;
  }

  ~FileBool() {
    //    std::cerr << "Clean destruction of FileBool associated to file=" <<
    //    file << "\n";
    std::fclose(fp);
  }

  // bool set/get
  bool getbit(size_t const &pos) {
    size_t i_byte = pos / 8;
    size_t i_pos = pos % 8;
    std::fseek(fp, i_byte, SEEK_SET);
    uint8_t val;
    size_t ret = std::fread(&val, sizeof(uint8_t), 1, fp);
    if (ret != 1) {
      std::cerr << "FileBool: Error in getbit(...)\n";
      std::cerr << "pos=" << pos << " i_byte=" << i_byte << "\n";
      std::cerr << "ret=" << ret << " file=" << file << "\n";
      throw TerminalException{1};
    }
    return val >> (i_pos & 0x07) & 1;
  }

  void setbit(size_t const &pos, bool const &val) {
    size_t curr_n_byte = (n_ent + 7) / 8;
    size_t needed_n_byte = (pos + 1 + 7) / 8;
    std::fseek(fp, curr_n_byte, SEEK_SET);
    uint8_t val_u8 = 0;
    for (size_t u = curr_n_byte; u < needed_n_byte; u++) {
      size_t n_write = std::fwrite(&val_u8, sizeof(uint8_t), 1, fp);
      if (n_write != 1) {
        std::cerr << "Error in setbit(..), wrong number of written bytes\n";
        throw TerminalException{1};
      }
    }
    // Now doing the assignment.
    size_t i_byte = pos / 8;
    size_t i_pos = pos % 8;
    std::fseek(fp, i_byte, SEEK_SET);
    size_t ret = std::fread(&val_u8, sizeof(uint8_t), 1, fp);
    if (ret != 1) {
      std::cerr << "FileBool: Error in setbit(...)\n";
      std::cerr << "pos=" << pos << " val=" << val << "\n";
      std::cerr << "curr_n_byte=" << curr_n_byte
                << " needed_n_byte=" << needed_n_byte << "\n";
      std::cerr << "i_byte=" << i_byte << " i_pos=" << i_pos << "\n";
      std::cerr << "ret=" << ret << " file=" << file << "\n";
      throw TerminalException{1};
    }
    std::fseek(fp, i_byte, SEEK_SET);
    val_u8 ^= static_cast<uint8_t>(-static_cast<uint8_t>(val) ^ val_u8) &
              kBitmask[i_pos];
    size_t n_write = std::fwrite(&val_u8, sizeof(uint8_t), 1, fp);
    if (n_write != 1) {
      std::cerr
          << "Error in setbit, number of written elements different from 1\n";
      throw TerminalException{1};
    }
    n_ent = std::max(n_ent, pos + 1);
  }

  void direct_write(std::vector<uint8_t> const &V) {
    size_t start_byte = 0;
    std::fseek(fp, start_byte, SEEK_SET);
    size_t len_w = V.size();
    if (len_w > 0) {
      size_t n_write = std::fwrite(V.data(), sizeof(uint8_t), len_w, fp);
      if (n_write != len_w) {
        std::cerr << "Error in direct_write n_write=" << n_write
                  << " len_w=" << len_w << "\n";
        throw TerminalException{1};
      }
    }
  }

  size_t getpopcnt() {
    size_t curr_n_byte = (n_ent + 7) / 8;
    std::fseek(fp, 0, SEEK_SET);
    std::vector<uint8_t> buffer(curr_n_byte);
    size_t ret = fread(buffer.data(), curr_n_byte * sizeof(uint8_t), 1, fp);
    if (ret != 1) {
      std::cerr << "FileBool: Error in getpopcnt(...)\n";
      throw TerminalException{1};
    }
    // popcnt
    size_t curr_n_ints = n_ent / 32;
    size_t sm = 0;
    uint32_t *ptr = (uint32_t *)(buffer.data());
    for (size_t i = 0; i < curr_n_ints; i++) {
      sm += __builtin_popcount(*(ptr + i));
    }
    // count remainder
    for (size_t i = 4 * curr_n_ints; i < curr_n_byte; i++) {
      uint8_t b = buffer[i];
      for (size_t j = 0; j < 8; j++) {
        if (8 * i + j < n_ent) {
          sm += (b & (uint8_t(1) << j)) ? 1 : 0;
        }
      }
    }
    return sm;
  }
};

// The use of std::vector<uint8_t> is quite not great.
// Could we represent that by the Face type? (from dynamic_byteset)
// The question is about the extendibility of the Face type.
// Clearly, this is a technical debt to pay later one way or the other.
std::vector<uint8_t> FileBool_FullRead(std::string const& eFile, size_t const& n_orbit) {
  FileBool fb(eFile, n_orbit);
  std::vector<uint8_t> l_status(n_orbit);
#ifdef DEBUG_BASIC_DATAFILE
  int sum_status = 0;
#endif
  for (size_t i=0; i<n_orbit; i++) {
    bool test = fb.getbit(i);
#ifdef DEBUG_BASIC_DATAFILE
    sum_status += static_cast<int>(test);
#endif
    uint8_t test_i = static_cast<uint8_t>(test);
    l_status[i] = test_i;
  }
#ifdef DEBUG_BASIC_DATAFILE
  os << "BASIC_DATAFILE: reading database sum_status=" << sum_status << "\n";
#endif
  return l_status;
}

void FileBool_FullWrite(std::string const& eFile, std::vector<uint8_t> const& l_status) {
  FileBool fb(eFile);
  size_t n_obj = l_status.size();
  for (size_t i=0; i<n_obj; i++) {
    bool test_done = static_cast<bool>(l_status[i]);
    fb.setbit(i, test_done);
  }
}




//
// Storing a sequence of faces that can be extended
//

struct FileFace {
private:
  std::FILE *fp;
  std::string file;
  size_t n_face;
  size_t siz;
  size_t ReadSize;
  std::vector<uint8_t> ReadBuffer;
  size_t ZeroSize;
  std::vector<uint8_t> ZeroBuffer;
  void SetBuffer() {
    ReadSize = (siz + 14) / 8;
    ReadBuffer = std::vector<uint8_t>(ReadSize);
    ZeroSize = (siz + 7) / 8;
    ZeroBuffer = std::vector<uint8_t>(ZeroSize, 0);
  }

public:
  FileFace(const FileFace &) = delete;
  FileFace &operator=(const FileFace &) = delete;
  FileFace(FileFace &&) = delete;
  FileFace() = delete;

  FileFace(std::string const &file, size_t const &_siz) : file(file) {
    fp = std::fopen(file.data(), "w+");
    siz = _siz;
    n_face = 0;
    SetBuffer();
  }

  FileFace(std::string const &file, size_t const &_siz, size_t const &_n_face)
      : file(file) {
    if (!IsExistingFile(file)) {
      std::cerr << "FileFace: The file " << file << " should not be missing\n";
      throw TerminalException{1};
    }
    fp = std::fopen(file.data(), "r+");
    siz = _siz;
    n_face = _n_face;
    SetBuffer();
  }

  ~FileFace() {
    //    std::cerr << "Clean destruction of FileFace associated to file=" <<
    //    file << "\n";
    std::fclose(fp);
  }

  size_t get_size() const { return siz; }

  // face set/get
  Face getface(size_t const &pos) {
    size_t start_byte = (pos * siz) / 8;
    size_t start_pos = (pos * siz) % 8;
    std::fseek(fp, start_byte, SEEK_SET);
    size_t end_byte = ((pos + 1) * siz + 7) / 8;
    size_t len = end_byte - start_byte;
    size_t ret = std::fread(ReadBuffer.data(), sizeof(uint8_t), len, fp);
    if (ret != len) {
      std::cerr << "getface: Number of read element different from count. pos="
                << pos << "\n";
      std::cerr << "Please correct. ret=" << ret << " len=" << len
                << " file=" << file << "\n";
      throw TerminalException{1};
    }
    //
    Face f(siz);
    size_t i_byte = 0;
    size_t i_pos = start_pos;
    for (size_t i = 0; i < siz; i++) {
      bool val = ReadBuffer[i_byte] >> (i_pos & 0x07) & 1;
      f[i] = val;
      i_pos++;
      if (i_pos == 8) {
        i_pos = 0;
        i_byte++;
      }
    }
    return f;
  }

  void setface(size_t const &pos, Face const &f) {
    // It is for the first empty byte that we may need to write.
    size_t curr_n_byte = (n_face * siz + 7) / 8;
    // It is the number of the byte just after the last one we need.
    size_t needed_n_byte = ((pos + 1) * siz + 7) / 8;
    std::fseek(fp, curr_n_byte, SEEK_SET);
    // Adding the zero entries.
    size_t len = 0;
    if (needed_n_byte > curr_n_byte)
      len = needed_n_byte - curr_n_byte;
    // We are in the standard case of appending some entries.
    if (len <= ZeroSize && len > 0) {
      size_t n_write = std::fwrite(ZeroBuffer.data(), sizeof(uint8_t), len, fp);
      if (n_write != len) {
        std::cerr << "Error in setface(..), wrong number of written bytes\n";
        throw TerminalException{1};
      }
    } else {
      // This case is rarer. It is for big jumps.
      uint8_t val_u8 = 0;
      for (size_t u = curr_n_byte; u < needed_n_byte; u++) {
        size_t n_write = std::fwrite(&val_u8, sizeof(uint8_t), 1, fp);
        if (n_write != 1) {
          std::cerr << "Error in setface(..), wrong number of written bytes\n";
          throw TerminalException{1};
        }
      }
    }
    // Now doing the assignment.
    size_t start_byte = (pos * siz) / 8;
    size_t start_pos = (pos * siz) % 8;
    size_t end_byte = ((pos + 1) * siz + 7) / 8;
    size_t len_rw = end_byte - start_byte;
    std::fseek(fp, start_byte, SEEK_SET);
    size_t ret = std::fread(ReadBuffer.data(), sizeof(uint8_t), len_rw, fp);
    if (ret != len_rw) {
      std::cerr << "setface: Number of read element different from count. pos="
                << pos << "\n";
      std::cerr << "Please correct. ret=" << ret << " len_rw=" << len_rw
                << " file=" << file << "\n";
      throw TerminalException{1};
    }
    // Assigning the new values
    size_t i_byte = 0;
    size_t i_pos = start_pos;
    for (size_t i = 0; i < siz; i++) {
      bool val = f[i];
      ReadBuffer[i_byte] ^= static_cast<uint8_t>(-static_cast<uint8_t>(val) ^
                                                 ReadBuffer[i_byte]) &
                            kBitmask[i_pos];
      i_pos++;
      if (i_pos == 8) {
        i_pos = 0;
        i_byte++;
      }
    }
    // Writing it out
    std::fseek(fp, start_byte, SEEK_SET);
    if (len_rw > 0) {
      size_t n_write =
          std::fwrite(ReadBuffer.data(), sizeof(uint8_t), len_rw, fp);
      if (n_write != len_rw) {
        std::cerr << "Error in setface(..), wrong number of written bytes\n";
        throw TerminalException{1};
      }
    }
    n_face = std::max(n_face, pos + 1);
  }

  void direct_write(std::vector<uint8_t> const &V) {
    size_t start_byte = 0;
    std::fseek(fp, start_byte, SEEK_SET);
    size_t len_w = V.size();
    if (len_w > 0) {
      size_t n_write = std::fwrite(V.data(), sizeof(uint8_t), len_w, fp);
      if (n_write != len_w) {
        std::cerr << "Error in direct_write n_write=" << n_write
                  << " len_w=" << len_w << "\n";
        throw TerminalException{1};
      }
    }
  }
};

//
// Storing a sequence of data that can be extended
//

size_t read_size_t(std::FILE* fp, size_t pos) {
  size_t ret_val;
  size_t* ptr1 = &ret_val;
  uint8_t* ptr2 = reinterpret_cast<uint8_t*>(ptr1);
  size_t pos_eff = pos * sizeof(size_t);
  std::fseek(fp, pos_eff, SEEK_SET);
  size_t n_read = std::fread(ptr2, sizeof(size_t), 1, fp);
  if (n_read != 1) {
    std::cerr << "n_read=" << n_read << "\n";
    std::cerr << "Failed to read the correct number of entries\n";
    throw TerminalException{1};
  }
  return ret_val;
}

void write_size_t(std::FILE* fp, size_t pos, size_t val) {
  size_t* ptr1 = &val;
  uint8_t* ptr2 = reinterpret_cast<uint8_t*>(ptr1);
  size_t pos_eff = pos * sizeof(size_t);
  std::fseek(fp, pos_eff, SEEK_SET);
  size_t n_write = std::fwrite(ptr2, sizeof(size_t), 1, fp);
  if (n_write != 1) {
    std::cerr << "n_write=" << n_write << "\n";
    std::cerr << "Failed to read the correct number of entries\n";
    throw TerminalException{1};
  }
}



template<typename T>
struct FileData {
private:
  size_t n_ent;
  size_t shift;
  std::FILE *fp_number;
  std::FILE *fp_data;
  std::string file;
  struct IteratorData {
    size_t n_ent;
    size_t pos;
    size_t shift;
    std::FILE *fp_number;
    std::FILE *fp_data;
    T read_state() {
#ifdef DEBUG_BASIC_DATAFILE
      std::cerr << "BASIC_DATAFILE: read_state, pos=" << pos << " n_ent=" << n_ent << "\n";
#endif
      size_t len = read_size_t(fp_number, pos+2);
#ifdef DEBUG_BASIC_DATAFILE
      std::cerr << "BASIC_DATAFILE: read_state, len=" << len << "\n";
#endif
      std::vector<char> V(len);
      std::fseek(fp_data, shift, SEEK_SET);
      size_t n_read = std::fread(V.data(), sizeof(uint8_t), len, fp_data);
#ifdef DEBUG_BASIC_DATAFILE
      std::cerr << "BASIC_DATAFILE: read_state, fread done\n";
#endif
      if (n_read != len) {
        std::cerr << "n_read=" << n_read << " len=" << len << "\n";
        std::cerr << "Failed to read the correct number of entries\n";
        throw TerminalException{1};
      }
      std::string str(V.data(), len);
#ifdef DEBUG_BASIC_DATAFILE
      std::cerr << "BASIC_DATAFILE: read_state, str built\n";
#endif
      std::istringstream iss(str);
#ifdef DEBUG_BASIC_DATAFILE
      std::cerr << "BASIC_DATAFILE: read_state, iss built\n";
#endif
      T val;
      boost::archive::text_iarchive ia(iss);
#ifdef DEBUG_BASIC_DATAFILE
      std::cerr << "BASIC_DATAFILE: read_state, ia built\n";
#endif
      ia >> val;
#ifdef DEBUG_BASIC_DATAFILE
      std::cerr << "BASIC_DATAFILE: read_state, returning\n";
#endif
      return val;
    }
    void single_increase() {
      size_t len = read_size_t(fp_number, pos+2);
      pos++;
      shift += len;
    }
    T operator*() { return read_state(); }
    IteratorData &operator++() {
      single_increase();
      return *this;
    }
    IteratorData operator++(int) {
      IteratorData tmp = *this;
      single_increase();
      return tmp;
    }
    bool operator!=(IteratorData const &iter) {
      return iter.pos != pos;
    }
    bool operator==(IteratorData const &iter) {
      return iter.pos == pos;
    }
  };

public:
  FileData(const FileFace &) = delete;
  FileData &operator=(const FileFace &) = delete;
  FileData(FileData &&) = delete;
  FileData() = delete;

  FileData(std::string const &file, bool overwrite) : file(file) {
    std::string file_number = file + ".nb";
    std::string file_data = file + ".data";
#ifdef DEBUG_BASIC_DATAFILE
    std::cerr << "BASIC_DATAFILE: FileData overwrite=" << overwrite << "\n";
#endif
    if (overwrite) {
      fp_number = std::fopen(file_number.data(), "w+");
      fp_data = std::fopen(file_data.data(), "w+");
      n_ent = 0;
      shift = 0;
      write_size_t(fp_number, 0, n_ent);
      write_size_t(fp_number, 1, shift);
    } else {
#ifdef DEBUG_BASIC_DATAFILE
      std::cerr << "BASIC_DATAFILE: FileData constructor ow=F step 1\n";
#endif
      if (!IsExistingFile(file_number) || !IsExistingFile(file_data)) {
        std::cerr << "FileData: The file " << file << " should not be missing\n";
        throw TerminalException{1};
      }
#ifdef DEBUG_BASIC_DATAFILE
      std::cerr << "BASIC_DATAFILE: FileData constructor ow=F step 2\n";
#endif
      fp_number = std::fopen(file_number.data(), "r+");
      fp_data = std::fopen(file_data.data(), "r+");
#ifdef DEBUG_BASIC_DATAFILE
      std::cerr << "BASIC_DATAFILE: FileData constructor ow=F step 3\n";
#endif
      // Reading the number of elements and the total shift
      n_ent = read_size_t(fp_number, 0);
      shift = read_size_t(fp_number, 1);
#ifdef DEBUG_BASIC_DATAFILE
      std::cerr << "BASIC_DATAFILE: FileData constructor ow=F step 4\n";
#endif
    }
  }

  ~FileData() {
    std::fclose(fp_number);
    std::fclose(fp_data);
  }

  void push_back(T const& val) {
    std::ostringstream ofs;
    boost::archive::text_oarchive oa(ofs);
    oa << val;
    std::string strO = ofs.str();
    size_t len = strO.size();
    write_size_t(fp_number, n_ent+2, len);
    std::fseek(fp_data, shift, SEEK_SET);
    shift += len;
    n_ent++;
    write_size_t(fp_number, 0, n_ent);
    write_size_t(fp_number, 1, shift);
    char* ptr = strO.data();
    size_t n_write = std::fwrite(ptr, sizeof(uint8_t), len, fp_data);
    if (n_write != len) {
      std::cerr << "n_write=" << n_write << " len=" << len << "\n";
      std::cerr << "Failed to read the correct number of entries\n";
      throw TerminalException{1};
    }
  }

  std::vector<size_t> read_all_sizes() {
    std::vector<size_t> l_sizes(n_ent);
    for (size_t i=0; i<n_ent; i++) {
      size_t len = read_size_t(fp_number, i + 2);
      l_sizes[i] = len;
    }
    return l_sizes;
  }

  // The iterator business
  using iterator = IteratorData;
  using const_iterator = IteratorData;
  const_iterator cbegin() const { return {n_ent, 0, 0, fp_number, fp_data}; }
  const_iterator cend() const { return {n_ent, n_ent, 0, fp_number, fp_data}; }
  const_iterator begin() const { return {n_ent, 0, 0, fp_number, fp_data}; }
  const_iterator end() const { return {n_ent, n_ent, 0, fp_number, fp_data}; }
};

template<typename T>
std::vector<T> FileData_FullRead(std::string const& FileDatabase) {
  bool overwrite = false;
  FileData<T> fdata(FileDatabase, overwrite);
  using Iterator = typename FileData<T>::iterator;
  Iterator iter = fdata.begin();
#ifdef DEBUG_BASIC_DATAFILE
  os << "BASIC_DATAFILE: reading database We have iter\n";
#endif
  std::vector<T> l_obj;
  while (iter != fdata.end()) {
#ifdef DEBUG_BASIC_DATAFILE
    os << "BASIC_DATAFILE: reading database Before *iter\n";
#endif
    T val = *iter;
#ifdef DEBUG_BASIC_DATAFILE
    os << "BASIC_DATAFILE: reading database After *iter\n";
#endif
    l_obj.emplace_back(std::move(val));
    iter++;
  }
  return l_obj;
}

template<typename T>
void FileData_FullWrite(std::string const& FileDatabase, std::vector<T> const& l_obj) {
  bool overwrite = true;
  FileData<T> fdata(FileDatabase, overwrite);
  for (auto & val : l_obj) {
    fdata.push_back(val);
  }
}


// clang-format off
#endif  //  SRC_DUALDESC_BASIC_DATAFILE_H_
// clang-format on

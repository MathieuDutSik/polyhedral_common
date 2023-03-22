// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_DUALDESC_BASIC_DATAFILE_H_
#define SRC_DUALDESC_BASIC_DATAFILE_H_

#include "Basic_file.h"
#include "Boost_bitset.h"
#include <algorithm>
#include <string>
#include <vector>

struct FileNumber {
private:
  std::FILE *fp;
  std::string file;

public:
  FileNumber(const FileNumber &) = delete;
  FileNumber &operator=(const FileNumber &) = delete;
  FileNumber(FileNumber &&) = delete;
  FileNumber() = delete;

  FileNumber(std::string const &file, bool is_new) : file(file) {
    if (is_new) {
      if (IsExistingFile(file)) {
        std::cerr << "FileNumber: The file " << file << " should be missing\n";
        throw TerminalException{1};
      }
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
    std::cerr << "Clean destruction of FileNumber associated to file=" << file << "\n";
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
    if (IsExistingFile(file)) {
      std::cerr << "FileBool: The file " << file << " should be missing\n";
      throw TerminalException{1};
    }
    fp = std::fopen(file.data(), "w+");
    n_ent = 0;
  }

  FileBool(std::string const &file, size_t const &_n_ent) {
    if (!IsExistingFile(file)) {
      std::cerr << "FileBool: The file " << file << " should not be missing\n";
      throw TerminalException{1};
    }
    fp = std::fopen(file.data(), "r+");
    n_ent = _n_ent;
  }

  ~FileBool() {
    std::cerr << "Clean destruction of FileBool associated to file=" << file << "\n";
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
};

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
    if (IsExistingFile(file)) {
      std::cerr << "FileFace: The file " << file << " should be missing\n";
      throw TerminalException{1};
    }
    fp = std::fopen(file.data(), "w+");
    siz = _siz;
    n_face = 0;
    SetBuffer();
  }

  FileFace(std::string const &file, size_t const &_siz, size_t const &_n_face) {
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
    std::cerr << "Clean destruction of FileFace associated to file=" << file << "\n";
    std::fclose(fp);
 }

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
    if (len <= ZeroSize) {
      if (len > 0) {
        size_t n_write =
            std::fwrite(ZeroBuffer.data(), sizeof(uint8_t), len, fp);
        if (n_write != len) {
          std::cerr << "Error in setface(..), wrong number of written bytes\n";
          throw TerminalException{1};
        }
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
    size_t n_write =
        std::fwrite(ReadBuffer.data(), sizeof(uint8_t), len_rw, fp);
    if (n_write != len_rw) {
      std::cerr << "Error in setface(..), wrong number of written bytes\n";
      throw TerminalException{1};
    }
    n_face = std::max(n_face, pos + 1);
  }
};

// clang-format off
#endif  //  SRC_DUALDESC_BASIC_DATAFILE_H_
// clang-format on

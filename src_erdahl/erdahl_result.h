// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_ERDAHL_ERDAHL_RESULT_H_
#define SRC_ERDAHL_ERDAHL_RESULT_H_

// clang-format off
#include "Temp_common.h"
#include <optional>
#include <utility>
#include <variant>
// clang-format on

/*
  The outcome of a computation that either succeeds with a value of type
  Tok or fails with an error of type Terr, in the manner of the Rust type
  Result<T, E>.

  The failures of the Erdahl computations are not programming errors: a
  quadratic function that takes a negative value on Z^n is an ordinary
  answer, and the lattice point where it is negative is exactly what the
  caller needs in order to continue (add a constraint, retry). So the error
  carries data and is returned, not thrown.
 */
template <typename Tok, typename Terr> struct ErdahlResult {
  std::variant<Tok, Terr> res;

  static ErdahlResult ok(Tok val) {
    return ErdahlResult{std::variant<Tok, Terr>(std::in_place_index<0>,
                                                std::move(val))};
  }
  static ErdahlResult err(Terr val) {
    return ErdahlResult{std::variant<Tok, Terr>(std::in_place_index<1>,
                                                std::move(val))};
  }

  bool is_ok() const { return res.index() == 0; }
  bool is_err() const { return res.index() == 1; }

  Tok const &get_ok() const {
    if (!is_ok()) {
      std::cerr << "ERDAHL: get_ok called on an error result\n";
      throw TerminalException{1};
    }
    return std::get<0>(res);
  }
  Terr const &get_err() const {
    if (!is_err()) {
      std::cerr << "ERDAHL: get_err called on an ok result\n";
      throw TerminalException{1};
    }
    return std::get<1>(res);
  }

  std::optional<Tok> get_result() const {
    if (is_ok()) {
      return std::get<0>(res);
    }
    return {};
  }
  std::optional<Terr> get_error() const {
    if (is_err()) {
      return std::get<1>(res);
    }
    return {};
  }
};

// clang-format off
#endif  // SRC_ERDAHL_ERDAHL_RESULT_H_
// clang-format on

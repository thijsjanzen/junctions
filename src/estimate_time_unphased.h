// Copyright 2018 - 2024 Thijs Janzen
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
#pragma once

#include <cmath>
#include <iostream>
#include <vector>

//
// from: https://rosettacode.org/wiki/Matrix-exponentiation_operator#C.2B.2B
//

template<int MSize = 7, class T = double >
class SqMx {
  typedef T Ax[MSize][MSize];
  typedef SqMx<MSize, T> Mx;

 private:
  Ax a;
  SqMx() { }

 public:
  explicit SqMx(const Ax &_a) {   // constructor with pre-defined array
    for (int r = 0; r < MSize; r++)
      for (int c = 0; c < MSize; c++)
        a[r][c] = _a[r][c];
  }

  static Mx identity() {
    Mx m;
    for (int r = 0; r < MSize; r++)
      for (int c = 0; c < MSize; c++)
        m.a[r][c] = (r == c ? 1 : 0);
    return m;
  }

  friend std::ostream &operator<<(std::ostream& os, const Mx &p) {
    for (int i = 0; i < MSize; i++) {
      for (int j = 0; j < MSize; j++)
        os << p.a[i][j] << ',';
      os << std::endl;
    }
    return os;
  }

  Mx operator*(const Mx &b) {
    Mx d;
    for (int r = 0; r < MSize; r++) {
      for (int c = 0; c < MSize; c++) {
        d.a[r][c] = 0;
        for (int k = 0; k < MSize; k++)
          d.a[r][c] += a[r][k] * b.a[k][c];
      }
    }
    return d;
  }

  Mx operator^(int n) {
    if (n < 0)
      throw "Negative exponent not implemented";

    Mx d = identity();
    for (Mx sq = *this; n > 0; sq = sq * sq, n /= 2) {
      if (n % 2 != 0) {
        d = d * sq;
      }
    }
    return d;
  }

  std::vector< double > operator()() {
    std::vector< double> output(MSize);
    for (int i = 0; i < MSize; ++i) {
      output[i] = static_cast<double>(a[0][i]);
    }
    return output;
  }
};

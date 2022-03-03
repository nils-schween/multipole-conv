// // Alpha matrix
// class AlphaMatrix {
//   // class invariant: roughly creates the correctly sized matrix and the
//   // introduces a relation between index sets and matrix entries
//  public:
//   AlphaMatrix(unsigned int l)
//       : degree{l},
//         matrix_size{(l + 2) * (l + 1)/2},
//         matrix_elements(matrix_size, std::vector<double>(matrix_size)),
//         index_map_lms(matrix_size),
//         index_map_pqr(matrix_size) {
//     create_index_maps();
//     compute_matrix_elements();
//   }

//   double& operator()(unsigned int i, unsigned j) {
//     return matrix_elements[i][j];
//   }

//   double& operator()(const std::array<unsigned int, 3>& lms,
//                      const std::array<unsigned int, 3>& pqr) {
//     unsigned int i = find_matrix_index(true, lms);
//     unsigned int j = find_matrix_index(false, pqr);
//     return matrix_elements[i][j];
//   }

//   void print();
//   void print_index_maps();
//   void print_tensor_components();

//  private:
//   void create_index_maps();
//   unsigned int find_matrix_index(
//       bool lms,  // lms: true, pqr: false
//       const std::array<unsigned int, 3>& index_tuple);
//   double expansion_coefficient(
//       const std::array<unsigned int, 3>& lms, unsigned int j, unsigned int k,
//       const std::array<unsigned int, 3>& trinomial_indices);
//   void compute_matrix_elements();
//   void print_matrix_element(unsigned int i, unsigned int j);

//   // data members
//   const unsigned int degree;
//   const unsigned int matrix_size;  // quadratic matrix n = m
//   std::vector<std::vector<double>> matrix_elements;
//   std::vector<std::array<unsigned int, 3>> index_map_lms;
//   std::vector<std::array<unsigned int, 3>> index_map_pqr;
// };

// void AlphaMatrix::create_index_maps() {
//   // lms index map
//   for (unsigned int l = (degree % 2 == 0 ? 0 : 1), i = 0; l <= degree; l += 2) {
//     for (unsigned int m = 0; m <= l; ++m) {
//       if (m == 0) {
//         index_map_lms[i] = {l, m, 0};
//         ++i;
//       } else {
//         index_map_lms[i] = {l, m, 0};
//         index_map_lms[i + 1] = {l, m, 1};
//         i += 2;
//       }
//     }
//   }
//   // pqr index map
//   for (unsigned int p = 0, i = 0; p <= degree; ++p) {
//     for (unsigned int q = 0; q <= degree - p; ++q, ++i) {
//       unsigned int r = degree - p - q;
//       index_map_pqr[i] = {p, q, r};
//     }
//   }
// }

// unsigned int AlphaMatrix::find_matrix_index(
//     bool lms, const std::array<unsigned int, 3>& index_tuple) {
//   if (lms) {
//     auto iterator =
//         std::find(index_map_lms.begin(), index_map_lms.end(), index_tuple);
//     if (iterator != index_map_lms.end()) {
//       return std::distance(index_map_lms.begin(), iterator);
//     } else {
//       assert(false && "Invalid l,m,s indices");
//       return 0;  // silence compiler warning: non-void not returning
//     }
//   } else {
//     auto iterator =
//         std::find(index_map_pqr.begin(), index_map_pqr.end(), index_tuple);
//     if (iterator != index_map_pqr.end())
//       return std::distance(index_map_pqr.begin(), iterator);
//     else {
//       assert(false && "Invalid p,q,r indices");
//       return 0;
//     }
//   }
// }

// void AlphaMatrix::print_index_maps() {
//   std::cout << "	Index maps:\n";
//   std::cout << "	i <-> (l, m, s)\n\n";
//   for (unsigned int i = 0; i < index_map_lms.size(); ++i) {
//     std::cout << "	" << i << " <-> (" << index_map_lms[i][0] << ", "
//               << index_map_lms[i][1] << ", " << index_map_lms[i][2] << ")\n";
//   }
//   std::cout << "\n";
//   std::cout << "	i <-> (p, q, r)\n\n";
//   for (unsigned int i = 0; i < index_map_pqr.size(); ++i) {
//     std::cout << "	" << i << " <-> (" << index_map_pqr[i][0] << ", "
//               << index_map_pqr[i][1] << ", " << index_map_pqr[i][2] << ")\n";
//   }
// }

// double AlphaMatrix::expansion_coefficient(
//     const std::array<unsigned int, 3>& lms, unsigned int j, unsigned int k,
//     const std::array<unsigned int, 3>& trinomial_indices) {
//   auto [l, m, s] = lms;
//   assert(l <= degree && m <= l && k <= m && j >= m && j <= l);
//   if (s)  // sin(m\phi)
//     return ALP_coefficient(l, m, j) * sin_mphi_coefficient(m, k) *
//            trinomial((degree - j)/2, trinomial_indices);
//   // cos(m\phi)
//   return ALP_coefficient(l, m, j) * cos_mphi_coefficient(m, k) *
//          trinomial((degree - j)/2, trinomial_indices);
// }

// void AlphaMatrix::compute_matrix_elements() {
//   // NOTE: It the programm is too slow, this loop could be the reason. The
//   // access operator operator()(lms, pqr) is calling std::find twice.
//   // Maybe a clever way to go through the indices would make this uncessary.
//   for (const auto& lms : index_map_lms) {
//     auto [l, m, s] = lms;
//     // 2 * ((l-m)/2) may look strange, but it uses integer division and takes
//     // care of operator precedence (* before /) to avoid calling std::floor. An
//     // equivalent expression is 2 * std::floor((l-m)/2)
//     for (unsigned int j = l - 2 * ((l - m) / 2); j <= l; j += 2) {
//       unsigned int sum_limit = (degree - j) / 2;
//       for (unsigned int k = 0; k <= m; ++k) {
//         for (unsigned int i_1 = 0; i_1 <= sum_limit; ++i_1) {
//           for (unsigned int i_2 = 0; i_2 <= sum_limit - i_1; ++i_2) {
//             unsigned int i_3 = sum_limit - i_1 - i_2;
//             std::array<unsigned int, 3> pqr{k + 2 * i_1, m - k + 2 * i_2,
//                                             j - m + 2 * i_3};
//             this->operator()(lms, pqr) +=
//                 expansion_coefficient(lms, j, k, {i_1, i_2, i_3});
//           }
//         }
//       }
//     }
//   }
// }

// void AlphaMatrix::print_matrix_element(unsigned int i, unsigned int j) {
//   const std::array<unsigned int, 3>& lms = index_map_lms[i];
//   const std::array<unsigned int, 3>& pqr = index_map_pqr[j];
//   std::cout << "alpha^" << lms[0] << lms[1] << lms[2] << "_" << pqr[0] << pqr[1]
//             << pqr[2] << " = " << matrix_elements[i][j];
// }

// void AlphaMatrix::print() {
//   for (unsigned int i = 0; i < matrix_size; ++i) {
//     std::cout << "	";
//     for (unsigned int j = 0; j < matrix_size; ++j) {
//       if (matrix_elements[i][j] != 0.) {
//         print_matrix_element(i, j);
//         std::cout << "    ";
//       }
//     }
//     std::cout << "\n";
//   }
// }

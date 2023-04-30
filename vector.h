//
// Created by Vladislav on 27/04/2023.
//

#ifndef MATRIX_VECTOR_H
#define MATRIX_VECTOR_H
#include "matrix.h"
namespace Algebra {
    template<typename value>
    class Vector : public Matrix<value> {
    public:
        Vector(const std::vector<value> &arr){
            Matrix<value>::n = arr.size();
            Matrix<value>::M.resize(Matrix<value>::n, std::vector<value> (Matrix<value>::m));
            for (size_t i = 0; i < Matrix<value>::n; ++i) {
                Matrix<value>::M[i][0] = arr[i];
            }
        }
        size_t Size() const {
            return Matrix<value>::get_column();
        }
        friend Matrix<value> Pack(const std::vector<Vector<value> > &arr) {
            for (size_t i = 1; i < arr.size(); ++i) {
                if (arr[i - 1].Size() != arr[i].Size()) {
                    Matrix<value>::fail("Diff sizes of vectors in Pack");
                }
            }
        }
    };
}
#endif //MATRIX_VECTOR_H

#ifndef GERSCHGORIN_HPP
#define GERSCHGORIN_HPP

#include "eigen.hpp"
#include <vector>
#include <utility> // std::pair

// GerschgorinSolver inherits EigenAnalysis (which inherits Matrix).
//
// Gerschgorin Circle Theorem:
//   For each row i, define:
//     centre  c_i = A(i,i)          (diagonal element)
//     radius  r_i = sum_{j != i} |A(i,j)|   (sum of absolute off-diagonal elements)
//
//   Then every eigenvalue of A lies in at least one Gerschgorin disc
//     D_i = { z in C : |z - c_i| <= r_i }
//   i.e. the eigenvalue interval on the real line is  [c_i - r_i, c_i + r_i].
//
// This class:
//   - computes all discs,
//   - unions them to obtain the overall bounding interval,
//   - prints per-disc info and the union interval.
class GerschgorinSolver : public EigenAnalysis
{
public:
    // A disc: centre, radius, and the derived real interval [lo, hi]
    struct Disc
    {
        int    row;
        double centre;
        double radius;
        double lo;   // centre - radius
        double hi;   // centre + radius
    };

private:
    std::vector<Disc> discs;   // one per row
    double unionLo;            // leftmost  boundary of union
    double unionHi;            // rightmost boundary of union
    bool   computed;

public:
    // Construct from an n x n matrix
    explicit GerschgorinSolver(const Matrix &A);

    // Compute all Gerschgorin discs
    void compute() override;

    // Print per-disc results and the union interval
    void printResults() const override;

    // Accessors
    const std::vector<Disc> &getDiscs()  const;
    double getUnionLo() const;
    double getUnionHi() const;
};

#endif

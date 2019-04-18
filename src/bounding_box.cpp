#include "Rcpp.h"

#include <deque>
#include <algorithm>
#include <stdexcept>
#include <string>

// [[Rcpp::export(rng=false)]]
Rcpp::List bounding_box(Rcpp::IntegerVector runs, Rcpp::StringVector values,
    Rcpp::IntegerVector index1, Rcpp::StringVector ref_chr1, Rcpp::IntegerVector ref_start1, Rcpp::IntegerVector ref_end1,
    Rcpp::IntegerVector index2, Rcpp::StringVector ref_chr2, Rcpp::IntegerVector ref_start2, Rcpp::IntegerVector ref_end2,
    Rcpp::LogicalVector reflect)
{
    const size_t ngroups=runs.size();
    if (runs.size()!=values.size()) {
        throw std::runtime_error("'runs' and 'values' should have the same length");
    }
    if (index1.size()!=index2.size()) {
        throw std::runtime_error("'index1' and 'index2' should have the same length");
    }
    if (reflect.size()!=1) {
        throw std::runtime_error("'reflect' should be a logical scalar");
    }
    const bool do_reflect=reflect[0];

    // Setting up output vectors.
    Rcpp::StringVector cout_left(ngroups), cout_right(ngroups);
    Rcpp::IntegerVector sout_left(ngroups), eout_left(ngroups), sout_right(ngroups), eout_right(ngroups);

    auto i1it=index1.begin();
    auto i2it=index2.begin();
    std::deque<int> left_starts, left_ends, right_starts, right_ends;

    for (size_t r=0; r<ngroups; ++r) {
        const int currun=runs[r];
        if (currun==0) { throw std::runtime_error("empty group"); }

        Rcpp::String left_chr=ref_chr1[*i1it], right_chr=ref_chr2[*i2it];
        if (left_chr > right_chr) {
            std::swap(left_chr, right_chr);
        }
        const bool intra=(left_chr==right_chr);

        for (int i=1; i<currun; ++i) {
            ++i1it;
            ++i2it;

            Rcpp::String chr1=ref_chr1[*i1it], chr2=ref_chr2[*i2it];
            const int start1=ref_start1[*i1it];
            const int start2=ref_start2[*i2it];
            const int end1=ref_end1[*i1it];
            const int end2=ref_end2[*i2it];

            if (chr1==left_chr && chr2==right_chr) {
                // Reflecting around the diagonal to ensure left anchor < right anchor.
                if (intra && do_reflect && (start1 > start2 || (start1==start2 && end1 > end2))) {
                    left_starts.push_back(start2);
                    left_ends.push_back(end2);
                    right_starts.push_back(start1);
                    right_ends.push_back(end1);
                } else {
                    left_starts.push_back(start1);
                    left_ends.push_back(end1);
                    right_starts.push_back(start2);
                    right_ends.push_back(end2);
                }
            } else if (chr1==right_chr && chr2==left_chr) {
                left_starts.push_back(start2);
                left_ends.push_back(end2);
                right_starts.push_back(start1);
                right_ends.push_back(end1);
            } else {
                throw std::runtime_error(std::string("multiple chromosomes for group '")
                    + Rcpp::as<std::string>(values[r]) + "'");
            }
        }

        // Figuring out max and min in the starts and ends.
        cout_left[r]=left_chr;
        cout_right[r]=right_chr;
        sout_left[r]=*std::min_element(left_starts.begin(), left_starts.end());
        eout_left[r]=*std::max_element(left_ends.begin(), left_ends.end());
        sout_right[r]=*std::min_element(right_starts.begin(), right_starts.end());
        eout_right[r]=*std::max_element(right_ends.begin(), right_ends.end());

        left_starts.clear();
        left_ends.clear();
        right_starts.clear();
        right_ends.clear();
    }
    
    return Rcpp::List::create(
        Rcpp::List::create(cout_left, sout_left, eout_left),
        Rcpp::List::create(cout_right, sout_right, eout_right)
    );
}
